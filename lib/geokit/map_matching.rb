module Geokit
  # Offline map matching (snapping GPS trajectories to a road network).
  #
  # Implements a Hidden Markov Model solved with the Viterbi algorithm,
  # following Newson & Krumm, "Hidden Markov Map Matching Through Noise
  # and Sparseness" (2009):
  #
  # * Emission model: GPS position error. The great-circle distance
  #   between an observation and its projection on a candidate road
  #   segment is modeled as a zero-mean Gaussian with caller-provided
  #   +sigma+ (meters).
  # * Transition model: the absolute difference between the road-network
  #   driving distance and the observed straight-line displacement is
  #   modeled as an exponential distribution with caller-provided +beta+
  #   (meters).
  #
  # Everything is computed locally: no geocoding, no map APIs and no
  # built-in road data. The caller provides:
  #
  # * timestamped observations (Geokit::MapMatching::Observation),
  # * candidate road segments per observation
  #   (Geokit::MapMatching::Candidate, e.g. via Edge#project),
  # * a router object answering road-network distances. It must respond
  #   to +route_distance(from_candidate, to_candidate)+ (or +call+) and
  #   return meters, or +nil+ when the transition is not drivable (for
  #   example against a one-way street or across disconnected parts of
  #   the network).
  #
  # Segmentation semantics: the matcher never fabricates a route through
  # the network. The trajectory is split into separate matched segments
  # whenever
  #
  # * an observation has no candidates (it is left unmatched),
  # * timestamps are not strictly increasing (duplicate or out-of-order
  #   observations start a new segment),
  # * the gap between two observations exceeds the optional
  #   +max_time_gap+,
  # * no candidate pair between two consecutive observations is
  #   reachable according to the router.
  module MapMatching
    # Local mean Earth radius (meters) used for polyline projection.
    EARTH_RADIUS_M = 6_371_000.0

    # Log probability of observing a GPS position +distance_m+ meters
    # away from the true position (Gaussian emission model).
    def self.emission_log_probability(distance_m, sigma)
      -Math.log(sigma * Math.sqrt(2.0 * Math::PI)) -
        0.5 * (distance_m / sigma)**2
    end

    # Log probability of the road-network distance +route_m+ given the
    # observed straight-line displacement +linear_m+ (exponential
    # transition model).
    def self.transition_log_probability(route_m, linear_m, beta)
      -Math.log(beta) - (route_m - linear_m).abs / beta
    end

    # Great-circle distance between two points, in meters.
    def self.distance_m(point_a, point_b)
      Geokit::LatLng.distance_between(point_a, point_b, units: :kms) * 1000.0
    end

    # A single timestamped GPS observation. +time+ may be a Time or any
    # comparable numeric (seconds); it is only used for ordering and for
    # gap comparison against Matcher's +max_time_gap+.
    class Observation
      attr_reader :point, :time

      def initialize(point, time)
        @point = Geokit::LatLng.normalize(point)
        @time = time
      end
    end

    # A road segment, described by an id and a polyline geometry
    # (array of lat/lng compatible values, at least two points).
    # Direction of travel and one-way restrictions are not interpreted
    # here; they belong to the router that answers route distances.
    class Edge
      attr_reader :id, :geometry

      def initialize(id, geometry)
        @id = id
        @geometry = geometry.map { |point| Geokit::LatLng.normalize(point) }
        if @geometry.size < 2
          raise ArgumentError, 'edge geometry needs at least two points'
        end
      end

      # Length of the polyline in meters.
      def length_m
        @length_m ||= @geometry.each_cons(2).sum do |from, to|
          MapMatching.distance_m(from, to)
        end
      end

      # Projects +point+ onto the polyline and returns a Candidate with
      # the projected position, its distance to the observation (meters)
      # and the fraction (0..1) of the edge traveled up to the
      # projection. Uses a local equirectangular approximation, accurate
      # for typical road-segment lengths.
      def project(point)
        point = Geokit::LatLng.normalize(point)
        lat0 = point.lat * Math::PI / 180.0
        lng0 = point.lng * Math::PI / 180.0
        cos0 = Math.cos(lat0)
        xy = @geometry.map do |vertex|
          [EARTH_RADIUS_M * cos0 * (vertex.lng * Math::PI / 180.0 - lng0),
           EARTH_RADIUS_M * (vertex.lat * Math::PI / 180.0 - lat0)]
        end

        best = nil
        traveled = 0.0
        xy.each_cons(2) do |a, b|
          dx = b[0] - a[0]
          dy = b[1] - a[1]
          seg_len_sq = dx * dx + dy * dy
          seg_len = Math.sqrt(seg_len_sq)
          t = seg_len_sq.zero? ? 0.0 : -(a[0] * dx + a[1] * dy) / seg_len_sq
          t = [[t, 0.0].max, 1.0].min
          proj_x = a[0] + t * dx
          proj_y = a[1] + t * dy
          dist = Math.hypot(proj_x, proj_y)
          if best.nil? || dist < best[:distance]
            best = {distance: dist, x: proj_x, y: proj_y,
                    along: traveled + t * seg_len}
          end
          traveled += seg_len
        end

        lat = point.lat + best[:y] / EARTH_RADIUS_M * 180.0 / Math::PI
        lng = point.lng + best[:x] / (EARTH_RADIUS_M * cos0) * 180.0 / Math::PI
        Candidate.new(edge: self, point: Geokit::LatLng.new(lat, lng),
                      distance: best[:distance],
                      fraction: best[:along] / traveled)
      end
    end

    # A projection of one observation onto one Edge.
    class Candidate
      # Projected position on the edge (Geokit::LatLng).
      attr_reader :edge, :point
      # Great-circle distance between observation and projection, meters.
      attr_reader :distance
      # Position along the edge, 0.0 (start) to 1.0 (end).
      attr_reader :fraction

      def initialize(edge:, point:, distance:, fraction:)
        @edge = edge
        @point = point
        @distance = distance
        @fraction = fraction
      end
    end

    # The HMM/Viterbi map matcher. All model parameters are supplied
    # explicitly by the caller:
    #
    # * +router+: object responding to
    #   +route_distance(from_candidate, to_candidate)+ (or +call+),
    #   returning driving distance in meters or +nil+ when unreachable.
    # * +sigma+: GPS position error standard deviation, meters (> 0).
    # * +beta+: transition model scale, meters (> 0).
    # * +max_time_gap+: optional; consecutive observations further apart
    #   than this (same unit as observation times, seconds for Time) are
    #   never connected and start a new segment.
    class Matcher
      attr_reader :router, :sigma, :beta, :max_time_gap

      def initialize(router:, sigma:, beta:, max_time_gap: nil)
        unless sigma.is_a?(Numeric) && sigma.positive?
          raise ArgumentError, 'sigma must be a positive number (meters)'
        end
        unless beta.is_a?(Numeric) && beta.positive?
          raise ArgumentError, 'beta must be a positive number (meters)'
        end
        if !max_time_gap.nil? &&
           !(max_time_gap.is_a?(Numeric) && max_time_gap.positive?)
          raise ArgumentError, 'max_time_gap must be nil or a positive number'
        end
        unless router.respond_to?(:route_distance) || router.respond_to?(:call)
          raise ArgumentError,
                'router must respond to route_distance(from, to) or call'
        end

        @router = router
        @sigma = sigma.to_f
        @beta = beta.to_f
        @max_time_gap = max_time_gap
      end

      # Matches a trajectory.
      #
      # +observations+ is an array of Observation, [point, time] pairs or
      # {point:, time:} hashes. +candidate_lists+ is a same-sized array
      # whose entries are arrays of Candidate (empty or nil when an
      # observation has no candidates).
      #
      # Returns a Result. The result holds one or more MatchedPath
      # segments; observations that could not be matched have a nil entry
      # in Result#matches.
      def match(observations, candidate_lists)
        obs = observations.map { |o| normalize_observation(o) }
        lists = candidate_lists.map { |list| list || [] }
        if obs.size != lists.size
          raise ArgumentError,
                'observations and candidate lists must have the same size'
        end

        segments = []
        aligned = Array.new(obs.size)
        columns = []

        obs.each_with_index do |observation, i|
          candidates = lists[i]
          if candidates.empty?
            finalize_run(columns, segments, aligned, obs, lists)
            next
          end

          emissions = candidates.map do |candidate|
            MapMatching.emission_log_probability(candidate.distance, sigma)
          end

          if columns.empty? || break_before?(obs, columns.last[:obs_index], i)
            finalize_run(columns, segments, aligned, obs, lists)
            columns << {obs_index: i, logps: emissions,
                        backpointers: Array.new(candidates.size)}
            next
          end

          column = extend_run(obs, lists, columns, i, candidates, emissions)
          if column
            columns << column
          else
            # No reachable transition: close the current segment and
            # start a new one instead of jumping through the network.
            finalize_run(columns, segments, aligned, obs, lists)
            columns << {obs_index: i, logps: emissions,
                        backpointers: Array.new(candidates.size)}
          end
        end
        finalize_run(columns, segments, aligned, obs, lists)

        Result.new(segments: segments, matches: aligned)
      end

      private

      def normalize_observation(value)
        case value
        when Observation then value
        when Array then Observation.new(value[0], value[1])
        when Hash
          point = value[:point] || value['point']
          time = value.key?(:time) ? value[:time] : value['time']
          Observation.new(point, time)
        else
          raise ArgumentError, "cannot interpret #{value.inspect} as an observation"
        end
      end

      def route_distance(from_candidate, to_candidate)
        if router.respond_to?(:route_distance)
          router.route_distance(from_candidate, to_candidate)
        else
          router.call(from_candidate, to_candidate)
        end
      end

      # True when observation +i+ must not be connected to the previous
      # candidate-bearing observation +j+: non-increasing (duplicate or
      # out-of-order) timestamps, or a gap beyond max_time_gap.
      def break_before?(obs, j, i)
        gap = obs[i].time - obs[j].time
        return true if gap <= 0
        !max_time_gap.nil? && gap > max_time_gap
      end

      # One Viterbi step between the last column and observation +i+.
      # Returns the new column, or nil when no candidate pair is
      # reachable.
      def extend_run(obs, lists, columns, i, candidates, emissions)
        prev = columns.last
        prev_index = prev[:obs_index]
        linear = MapMatching.distance_m(obs[prev_index].point, obs[i].point)
        logps = Array.new(candidates.size, -Float::INFINITY)
        backpointers = Array.new(candidates.size)
        reachable = false

        candidates.each_with_index do |candidate_b, bi|
          best = -Float::INFINITY
          best_ai = nil
          lists[prev_index].each_with_index do |candidate_a, ai|
            prev_logp = prev[:logps][ai]
            next if prev_logp == -Float::INFINITY
            route_m = route_distance(candidate_a, candidate_b)
            next if route_m.nil?
            logp = prev_logp +
                   MapMatching.transition_log_probability(route_m, linear, beta)
            if logp > best
              best = logp
              best_ai = ai
            end
          end
          unless best_ai.nil?
            reachable = true
            logps[bi] = best + emissions[bi]
            backpointers[bi] = best_ai
          end
        end

        return nil unless reachable

        {obs_index: i, logps: logps, backpointers: backpointers}
      end

      # Backtracks the best path through the accumulated Viterbi columns
      # and appends a MatchedPath to +segments+.
      def finalize_run(columns, segments, aligned, obs, lists)
        return if columns.empty?

        last = columns.last
        best_index = last[:logps].each_with_index.max_by { |logp, _| logp }[1]
        score = last[:logps][best_index]
        matches = []
        index = best_index
        columns.reverse_each do |column|
          obs_index = column[:obs_index]
          match = Match.new(index: obs_index, observation: obs[obs_index],
                            candidate: lists[obs_index][index])
          aligned[obs_index] = match
          matches.unshift(match)
          index = column[:backpointers][index]
        end
        segments << MatchedPath.new(matches: matches, score: score)
        columns.clear
      end
    end

    # The result of Matcher#match.
    class Result
      # Array of MatchedPath, in chronological order.
      attr_reader :segments

      def initialize(segments:, matches:)
        @segments = segments
        @matches = matches
      end

      # Array aligned with the input observations: a Match for every
      # matched observation, nil for observations without candidates.
      attr_reader :matches

      # Sum of the segment log-probability scores.
      def score
        segments.sum(&:score)
      end

      def success?
        segments.any?
      end
    end

    # One connected piece of a matched trajectory: the most likely
    # sequence of road segments for a run of consecutive observations.
    class MatchedPath
      # Array of Match, in chronological order.
      attr_reader :matches
      # Total log probability (emission + transition) of this path.
      attr_reader :score

      def initialize(matches:, score:)
        @matches = matches
        @score = score
      end

      # The matched road segment sequence (one Edge per observation).
      def edges
        matches.map(&:edge)
      end

      def observations
        matches.map(&:observation)
      end

      def size
        matches.size
      end
    end

    # The match of a single observation: its winning candidate and
    # therefore its projected position on the matched road network.
    class Match
      # Position of the observation in the input array.
      attr_reader :index
      attr_reader :observation, :candidate

      def initialize(index:, observation:, candidate:)
        @index = index
        @observation = observation
        @candidate = candidate
      end

      # Projected (snapped) position on the matched road segment.
      def point
        candidate.point
      end

      def edge
        candidate.edge
      end
    end
  end
end
