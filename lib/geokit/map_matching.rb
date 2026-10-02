module Geokit
  # Offline map matching (map snapping) for timestamped trajectories.
  #
  # Implements a Hidden Markov Model matcher in the spirit of
  # Newson & Krumm (2009): GPS position error is the emission model and
  # the discrepancy between network travel distance and observed
  # displacement is the transition model. The most likely sequence of
  # road segments is decoded with the Viterbi algorithm.
  #
  # Everything is offline and caller-supplied: observations, candidate
  # road segments per observation, and a network distance oracle. No map
  # API is contacted and no road data is baked in.
  #
  # Example:
  #   matcher = Geokit::MapMatching::Matcher.new(
  #     network: my_network, # responds to route_distance(from, to)
  #     sigma: 10.0,         # GPS position error std-dev, in +units+
  #     beta: 40.0,          # transition tolerance scale, in +units+
  #     max_time_gap: 300,   # seconds; longer gaps split the trajectory
  #     units: :meters
  #   )
  #   result = matcher.match(observations, candidates)
  module MapMatching
    # A single timestamped GPS observation.
    # +time+ is anything responding to +to_f+ (Time, DateTime, Numeric);
    # +point+ is any Geokit::LatLng-normalizable value.
    Observation = Struct.new(:time, :point)

    # A road segment. +id+ is caller-chosen, +polyline+ is an ordered
    # array of Geokit::LatLng-normalizable points, +one_way+ marks a
    # directed segment. Directionality is enforced by the network oracle,
    # not by the matcher itself.
    Segment = Struct.new(:id, :polyline, :one_way)

    # A road segment projected against one observation. Instances are
    # built by the matcher and handed to the network oracle, so the
    # oracle can locate the projection along the segment via +offset+
    # (distance from the polyline start, in the matcher's units) or
    # +fraction+ (0..1 along the polyline).
    class Candidate
      attr_reader :segment, :point, :distance, :offset, :fraction

      def initialize(segment:, point:, distance:, offset:, fraction:)
        @segment  = segment
        @point    = point     # projected position on the polyline (LatLng)
        @distance = distance  # perpendicular distance observation->projection
        @offset   = offset    # distance along the polyline to the projection
        @fraction = fraction  # offset / polyline length (0 for point segments)
      end
    end

    # One observation matched to one candidate.
    class MatchedPoint
      attr_reader :observation, :candidate

      def initialize(observation, candidate)
        @observation = observation
        @candidate = candidate
      end

      def segment
        candidate.segment
      end

      # Projected (snapped) position of the observation on the segment.
      def point
        candidate.point
      end

      # Perpendicular distance between observation and projection.
      def distance
        candidate.distance
      end
    end

    # A maximal run of consecutively matched observations connected by
    # reachable transitions. A MatchResult holds one or more of these;
    # it never contains a path that drives through walls.
    class MatchedSegment
      attr_reader :points, :score

      def initialize(points, score)
        @points = points
        @score = score
      end

      def segments
        points.map(&:segment)
      end
    end

    # Why two neighbouring matched pieces are not connected.
    # +kind+ is one of :no_candidates (the observation had no candidate
    # segments and was dropped), :time_gap (observations further apart
    # than max_time_gap), or :unreachable (no network route connects the
    # best candidates of the two observations, e.g. a one-way street
    # driven against its direction). +observation+ is the observation
    # that caused the break (for :time_gap and :unreachable, the first
    # observation after the break).
    Break = Struct.new(:kind, :observation)

    # The outcome of a match run.
    class MatchResult
      attr_reader :segments, :breaks

      def initialize(segments, breaks)
        @segments = segments
        @breaks = breaks
      end

      # All matched points in chronological order.
      def points
        segments.flat_map(&:points)
      end

      def matched?
        !segments.empty?
      end
    end

    class Error < StandardError; end

    # Raised when two observations share the same timestamp: the
    # transition model is undefined for zero-duration steps, so this is
    # rejected instead of silently producing nonsense.
    class DuplicateTimestampError < Error; end

    class Matcher
      attr_reader :sigma, :beta, :max_time_gap, :units

      # network: object responding to
      #   route_distance(from_candidate, to_candidate) -> Numeric or nil
      # (a callable is also accepted). Distances are in +units+; nil or
      # Float::INFINITY means not reachable.
      def initialize(network:, sigma:, beta:, max_time_gap: nil, units: :meters)
        raise ArgumentError, 'sigma must be positive' unless sigma.to_f.positive?
        raise ArgumentError, 'beta must be positive' unless beta.to_f.positive?

        @router = self.class.normalize_router(network)
        @sigma = sigma.to_f
        @beta = beta.to_f
        @max_time_gap = max_time_gap
        @units = units
        # Touch the unit tables so an unsupported unit fails fast.
        Geokit::LatLng.units_per_latitude_degree(units) ||
          raise(ArgumentError, "#{units} is an unsupported unit of length.")
      end

      def self.normalize_router(network)
        return network.method(:route_distance) if network.respond_to?(:route_distance)
        return network if network.respond_to?(:call)

        raise ArgumentError, 'network must respond to route_distance(from, to) or to call'
      end

      # observations: array of Observation (must have distinct times).
      # candidates: array parallel to +observations+; element i is the
      # array of candidate Segments for observation i (may be empty).
      def match(observations, candidates)
        unless observations.length == candidates.length
          raise ArgumentError, 'observations and candidates must have the same length'
        end

        entries = observations.zip(candidates).map do |obs, segs|
          [normalize_observation(obs), segs]
        end
        entries.sort_by! { |obs, _| observation_time(obs) }
        check_duplicate_times!(entries)

        breaks = []
        windows = split_windows(entries, breaks)
        segments = windows.flat_map { |window| match_window(window, breaks) }
        breaks.sort_by! { |brk| observation_time(brk.observation) }
        MatchResult.new(segments, breaks)
      end

      private

      def normalize_observation(obs)
        Observation.new(obs.time, Geokit::LatLng.normalize(obs.point))
      end

      def observation_time(obs)
        time = obs.time
        time = time.to_f if time.respond_to?(:to_f)
        unless time.is_a?(Numeric)
          raise ArgumentError, "observation time #{obs.time.inspect} is not numeric"
        end

        time.to_f
      end

      def check_duplicate_times!(entries)
        entries.each_cons(2) do |(prev, _), (cur, _)|
          next unless observation_time(prev) == observation_time(cur)

          raise DuplicateTimestampError,
                "duplicate observation timestamp #{cur.time.inspect}"
        end
      end

      # Groups observations into windows separated by unmatchable
      # observations (no candidates) and by long signal gaps.
      def split_windows(entries, breaks)
        windows = []
        current = []
        previous_time = nil
        entries.each do |obs, segs|
          time = observation_time(obs)
          if segs.empty?
            breaks << Break.new(:no_candidates, obs)
            windows << current unless current.empty?
            current = []
          else
            if @max_time_gap && previous_time &&
               time - previous_time > @max_time_gap
              breaks << Break.new(:time_gap, obs)
              windows << current unless current.empty?
              current = []
            end
            current << [obs, project_candidates(obs, segs)]
          end
          previous_time = time
        end
        windows << current unless current.empty?
        windows
      end

      # Runs Viterbi over one window, then cuts the decoded path wherever
      # the best transition was unreachable, yielding wall-free segments.
      def match_window(window, breaks)
        path, unreachable = viterbi(window)
        segments = []
        points = []
        score = 0.0
        path.each_with_index do |(candidate, emission, transition), i|
          if i.positive? && unreachable[i]
            segments << MatchedSegment.new(points, score)
            points = []
            score = 0.0
            breaks << Break.new(:unreachable, window[i].first)
          end
          points << MatchedPoint.new(window[i].first, candidate)
          score += emission
          score += transition if transition
        end
        segments << MatchedSegment.new(points, score) unless points.empty?
        segments
      end

      # Returns [path, unreachable_flags]. Each path element is
      # [candidate, emission_log_score, transition_log_score_or_nil].
      def viterbi(window)
        first_candidates = window.first[1]
        dp = first_candidates.map { |c| emission_score(c.distance) }
        back = []

        window.each_cons(2) do |(prev_obs, prev_candidates), (obs, candidates)|
          displacement = Geokit::LatLng.distance_between(
            prev_obs.point, obs.point, units: @units
          )
          step_back = []
          new_dp = candidates.map do |candidate|
            emission = emission_score(candidate.distance)
            best = nil
            prev_candidates.each_with_index do |prev_candidate, i|
              transition = transition_score(prev_candidate, candidate, displacement)
              total = dp[i] + transition
              best = [total, i, transition] if best.nil? || total > best[0]
            end
            step_back << best[1]
            best[0] + emission
          end
          back << step_back
          dp = new_dp
        end

        backtrace(window, back, dp)
      end

      def backtrace(window, back, dp)
        index = dp.each_with_index.max_by { |score, _| score }[1]
        path = Array.new(window.length)
        unreachable = Array.new(window.length, false)
        (window.length - 1).downto(0) do |i|
          candidate = window[i][1][index]
          emission = emission_score(candidate.distance)
          transition = nil
          if i.positive?
            prev_index = back[i - 1][index]
            displacement = Geokit::LatLng.distance_between(
              window[i - 1][0].point, window[i][0].point, units: @units
            )
            transition = transition_score(
              window[i - 1][1][prev_index], candidate, displacement
            )
            unreachable[i] = transition == -Float::INFINITY
            index = prev_index
          end
          path[i] = [candidate, emission, transition]
        end
        [path, unreachable]
      end

      # Emission model: zero-mean Gaussian position error.
      def emission_score(distance)
        -Math.log(@sigma * Math.sqrt(2.0 * Math::PI)) -
          0.5 * (distance / @sigma)**2
      end

      # Transition model: exponential decay on the absolute difference
      # between network route distance and observed displacement.
      # Unreachable pairs score -Infinity.
      def transition_score(from_candidate, to_candidate, displacement)
        route = @router.call(from_candidate, to_candidate)
        return -Float::INFINITY if route.nil? || route == Float::INFINITY
        unless route.is_a?(Numeric)
          raise ArgumentError, "network returned non-numeric distance #{route.inspect}"
        end

        -Math.log(@beta) - (route - displacement).abs / @beta
      end

      # Projects the observation onto each candidate segment polyline.
      def project_candidates(obs, segments)
        segments.map do |segment|
          point, distance, offset, fraction = project(obs.point, segment)
          Candidate.new(
            segment: segment, point: point,
            distance: distance, offset: offset, fraction: fraction
          )
        end
      end

      # Point-to-polyline projection in a local tangent plane around the
      # observation, using Geokit's own units-per-degree conversions.
      # Returns [projected LatLng, distance, offset, fraction].
      def project(point, segment)
        polyline = normalized_polyline(segment)
        per_lat = Geokit::LatLng.units_per_latitude_degree(@units)
        per_lng = Geokit::LatLng.units_per_longitude_degree(point.lat, @units)
        xy = polyline.map do |p|
          [(p.lng - point.lng) * per_lng, (p.lat - point.lat) * per_lat]
        end

        best = nil # [dist_sq, along, px, py]
        cumulative = 0.0
        xy.each_cons(2) do |(ax, ay), (bx, by)|
          dx = bx - ax
          dy = by - ay
          length_sq = dx * dx + dy * dy
          t = length_sq.zero? ? 0.0 : -(ax * dx + ay * dy) / length_sq
          t = 0.0 if t.negative?
          t = 1.0 if t > 1.0
          px = ax + t * dx
          py = ay + t * dy
          dist_sq = px * px + py * py
          along = cumulative + Math.sqrt(length_sq) * t
          best = [dist_sq, along, px, py] if best.nil? || dist_sq < best[0]
          cumulative += Math.sqrt(length_sq)
        end
        if best.nil? # single-point polyline
          best = [xy[0][0]**2 + xy[0][1]**2, 0.0, xy[0][0], xy[0][1]]
        end

        length = segment_length(segment, polyline)
        projected = Geokit::LatLng.new(
          point.lat + best[3] / per_lat,
          point.lng + best[2] / per_lng
        )
        fraction = length.zero? ? 0.0 : best[1] / length
        [projected, Math.sqrt(best[0]), best[1], fraction]
      end

      def normalized_polyline(segment)
        @polylines ||= {}
        @polylines[segment.object_id] ||=
          Array(segment.polyline).map { |p| Geokit::LatLng.normalize(p) }
      end

      def segment_length(segment, polyline)
        @segment_lengths ||= {}
        @segment_lengths[segment.object_id] ||= polyline.each_cons(2).sum do |a, b|
          Geokit::LatLng.distance_between(a, b, units: @units)
        end
      end
    end
  end
end
