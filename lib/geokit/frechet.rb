module Geokit
  # Computes the *continuous* Frechet distance between two polyline
  # trajectories.
  #
  # The Frechet distance is commonly illustrated as the shortest leash
  # that allows a person and a dog to walk along the two curves from
  # start to finish. Either one may stop and wait, but neither may walk
  # backwards. Unlike a discrete (vertex-to-vertex dynamic programming)
  # distance, the continuous variant considers every point of both
  # polylines, so subdividing a segment into smaller pieces never
  # changes the result, and unlike the Hausdorff distance it respects
  # the order and direction in which the points are traversed.
  #
  # The library works on the sphere, the Frechet distance does not:
  # it needs a metric plane. The caller therefore supplies a local
  # projection that maps each Geokit::LatLng to planar [x, y]
  # coordinates (meters). Geokit deliberately does not pick a map
  # extent or projection for you; use equirectangular_projection (or
  # any callable of your own) with an origin that suits your data.
  #
  #   projection = Geokit::Frechet.equirectangular_projection(
  #     Geokit::LatLng.new(37.0, -122.0)
  #   )
  #   result = Geokit::Frechet.distance(track_a, track_b,
  #                                     projection: projection)
  #   result.distance     # best estimate, in meters
  #   result.precision    # width of the certified [lower, upper] bracket
  #   result.witness      # monotone pairing of arc-length positions
  #
  # The solver is the classical Alt & Godau free-space-diagram
  # construction combined with a binary search on the leash length. It
  # returns certified lower and upper bounds instead of pretending that
  # a fixed sampling density is the continuous answer.
  module Frechet
    DEFAULT_TOLERANCE = 1e-6 # meters
    MAX_ITERATIONS = 100

    # Outcome of a Frechet computation. +lower_bound+ and +upper_bound+
    # bracket the true continuous Frechet distance; +witness+ is a
    # monotone chain of [s, t] pairs (arc-length positions in meters
    # along the first and second trajectory) starting at [0, 0] and
    # ending at [length_a, length_b]. Linear interpolation between
    # consecutive witness pairs yields a continuous, monotone matching
    # whose leash never exceeds +upper_bound+.
    class Result
      attr_reader :lower_bound, :upper_bound, :witness

      def initialize(lower_bound:, upper_bound:, witness:)
        @lower_bound = lower_bound
        @upper_bound = upper_bound
        @witness = witness
      end

      # Best estimate of the Frechet distance (midpoint of the bracket).
      def distance
        (lower_bound + upper_bound) / 2.0
      end

      # Achieved precision: upper_bound - lower_bound.
      def precision
        upper_bound - lower_bound
      end
    end

    # Computes the continuous Frechet distance between two trajectories.
    #
    # +path_a+ and +path_b+ are non-empty enumerables of points (typically
    # Geokit::LatLng). +projection+ is mandatory and must respond to
    # +call+, mapping a point to planar [x, y] coordinates in meters.
    # +tolerance+ is the absolute precision (in meters) to which the
    # distance is certified.
    #
    # Raises ArgumentError for empty trajectories, a missing/invalid
    # projection, or a non-positive tolerance.
    def self.distance(path_a, path_b, projection:, tolerance: DEFAULT_TOLERANCE,
                      max_iterations: MAX_ITERATIONS)
      unless projection.respond_to?(:call)
        raise ArgumentError, 'projection must respond to #call'
      end
      raise ArgumentError, 'tolerance must be positive' unless tolerance.to_f > 0

      a = project_path(path_a, projection, 'path_a')
      b = project_path(path_b, projection, 'path_b')
      Solver.new(a, b).solve(tolerance: tolerance.to_f, max_iterations: max_iterations)
    end

    # Builds a simple equirectangular local projection around an
    # explicitly chosen origin. Suitable for small areas (city scale);
    # applications that need something fancier should pass their own
    # callable to .distance instead.
    def self.equirectangular_projection(origin)
      origin = Geokit::LatLng.normalize(origin)
      cos_lat = Math.cos(origin.lat * Math::PI / 180.0)
      meters_per_degree =
        Geokit::Mappable::ClassMethods::EARTH_RADIUS_IN_METERS * Math::PI / 180.0
      lambda do |point|
        ll = Geokit::LatLng.normalize(point)
        [
          (ll.lng - origin.lng) * meters_per_degree * cos_lat,
          (ll.lat - origin.lat) * meters_per_degree
        ]
      end
    end

    def self.project_path(path, projection, name)
      raise ArgumentError, "#{name} must not be empty" if path.nil? || path.empty?

      path.map do |point|
        projected = projection.call(point)
        unless projected.respond_to?(:to_a) && projected.to_a.size == 2
          raise ArgumentError, "projection must return [x, y], got #{projected.inspect}"
        end

        x, y = projected.to_a
        [x.to_f, y.to_f]
      end
    end
    private_class_method :project_path

    # Alt & Godau decision procedure + parametric binary search.
    class Solver #:nodoc: all
      INFINITY = Float::INFINITY

      # Relative/absolute slack used inside the decision procedure so
      # that floating-point noise at exactly-tight leashes cannot break
      # the reachability chain. The reported upper bound always includes
      # this slack, so the certified bracket stays honest.
      SLACK = 1e-9

      # Tolerance (in vertex-index parameter units) used when
      # backtracking the witness, to absorb floating-point noise at
      # cell boundaries.
      WITNESS_TOL = 1e-9

      def initialize(a, b)
        @a = a
        @b = b
        @lengths_a = segment_lengths(a)
        @lengths_b = segment_lengths(b)
        @cum_a = cumulative(@lengths_a)
        @cum_b = cumulative(@lengths_b)
      end

      def solve(tolerance:, max_iterations:)
        return degenerate_result if @a.size == 1 || @b.size == 1

        lo = [distance(@a.first, @b.first), distance(@a.last, @b.last)].max
        hi = max_pairwise_distance
        return exact_result(lo) if hi <= lo

        iterations = 0
        while with_slack(hi) - lo > tolerance && iterations < max_iterations
          mid = (lo + hi) / 2.0
          if feasible?(mid)
            hi = mid
          else
            lo = mid
          end
          iterations += 1
        end
        Result.new(lower_bound: lo, upper_bound: with_slack(hi),
                   witness: build_witness(hi))
      end

      private

      def segment_lengths(points)
        points.each_cons(2).map { |p, q| distance(p, q) }
      end

      def cumulative(lengths)
        cum = [0.0]
        lengths.each { |length| cum << cum[-1] + length }
        cum
      end

      def distance(p, q)
        Math.hypot(p[0] - q[0], p[1] - q[1])
      end

      def max_pairwise_distance
        hi = 0.0
        @a.each do |p|
          @b.each do |q|
            d = distance(p, q)
            hi = d if d > hi
          end
        end
        hi
      end

      # One of the trajectories is a single point: the answer is exact
      # (the point waits while the other curve is traversed), no search.
      def degenerate_result
        if @a.size == 1 && @b.size == 1
          d = distance(@a[0], @b[0])
          Result.new(lower_bound: d, upper_bound: d, witness: [[0.0, 0.0]])
        elsif @a.size == 1
          d = @b.map { |q| distance(@a[0], q) }.max
          Result.new(lower_bound: d, upper_bound: d,
                     witness: [[0.0, 0.0], [0.0, @cum_b[-1]]])
        else
          d = @a.map { |p| distance(p, @b[0]) }.max
          Result.new(lower_bound: d, upper_bound: d,
                     witness: [[0.0, 0.0], [@cum_a[-1], 0.0]])
        end
      end

      def exact_result(value)
        Result.new(lower_bound: value, upper_bound: value,
                   witness: build_witness(value))
      end

      def with_slack(eps)
        eps * (1.0 + SLACK) + 1e-12
      end

      # --- Free-space diagram decision procedure ---------------------

      # Interval of u in [0, 1] where the point travelling from p to q
      # stays within +eps+ of the fixed point. Returns [lo, hi] or nil.
      def free_interval(fixed, p, q, eps)
        dx = q[0] - p[0]
        dy = q[1] - p[1]
        fx = p[0] - fixed[0]
        fy = p[1] - fixed[1]
        a = dx * dx + dy * dy
        if a.zero?
          return fx * fx + fy * fy <= eps * eps ? [0.0, 1.0] : nil
        end

        b = 2.0 * (dx * fx + dy * fy)
        c = fx * fx + fy * fy - eps * eps
        disc = b * b - 4.0 * a * c
        scale = b * b + 4.0 * a * (c.abs + eps * eps + 1.0)
        return nil if disc < -1e-9 * scale

        root = Math.sqrt([disc, 0.0].max)
        lo = (-b - root) / (2.0 * a)
        hi = (-b + root) / (2.0 * a)
        lo = 0.0 if lo < 1e-12
        hi = 1.0 if hi > 1.0 - 1e-12
        lo = 0.0 if lo < 0.0
        hi = 1.0 if hi > 1.0
        lo <= hi ? [lo, hi] : nil
      end

      # Intersection of +interval+ with [low, 1]; nil when empty.
      def clip_from(interval, low)
        return nil unless interval

        lo = interval[0] > low ? interval[0] : low
        lo <= interval[1] ? [lo, interval[1]] : nil
      end

      def interval_low(interval)
        interval ? interval[0] : INFINITY
      end

      # Fills @lr/@br (reachable left/bottom edges) and @rr/@tr
      # (reachable right/top edges) for every cell and reports whether
      # the top-right corner of the diagram is reachable, i.e. whether
      # the Frechet distance is at most +eps+.
      def feasible?(eps)
        eps = with_slack(eps)
        m = @a.size - 1
        n = @b.size - 1
        @lr = Array.new(m) { Array.new(n) }
        @br = Array.new(m) { Array.new(n) }
        @rr = Array.new(m) { Array.new(n) }
        @tr = Array.new(m) { Array.new(n) }

        m.times do |i|
          n.times do |j|
            left_free = free_interval(@a[i], @b[j], @b[j + 1], eps)
            bottom_free = free_interval(@b[j], @a[i], @a[i + 1], eps)

            left_in = i.zero? ? (j.zero? ? left_free : nil) : @rr[i - 1][j]
            bottom_in = j.zero? ? (i.zero? ? bottom_free : nil) : @tr[i][j - 1]

            left_low = interval_low(left_in)
            left_low = 0.0 if bottom_in && bottom_in[0] <= 0.0
            bottom_low = interval_low(bottom_in)
            bottom_low = 0.0 if left_in && left_in[0] <= 0.0

            left_reach = clip_from(left_free, left_low)
            bottom_reach = clip_from(bottom_free, bottom_low)

            right_free = free_interval(@a[i + 1], @b[j], @b[j + 1], eps)
            top_free = free_interval(@b[j + 1], @a[i], @a[i + 1], eps)

            right_low = interval_low(left_reach)
            right_low = 0.0 if bottom_reach
            top_low = interval_low(bottom_reach)
            top_low = 0.0 if left_reach

            @lr[i][j] = left_reach
            @br[i][j] = bottom_reach
            @rr[i][j] = clip_from(right_free, right_low)
            @tr[i][j] = clip_from(top_free, top_low)
          end
        end

        last_right = @rr[m - 1][n - 1]
        last_top = @tr[m - 1][n - 1]
        (last_right && last_right[1] >= 1.0) || (last_top && last_top[1] >= 1.0) ? true : false
      end

      # --- Witness reconstruction ------------------------------------

      # Backtracks a monotone path through the reachable free space at
      # +eps+ (which must be feasible). Returns a monotone chain of
      # [s, t] arc-length pairs in meters.
      def build_witness(eps)
        feasible?(eps)
        m = @a.size - 1
        n = @b.size - 1
        chain = []
        gi = m.to_f
        gj = n.to_f
        loop do
          chain << [gi, gj]
          break if gi.zero? && gj.zero?

          gi, gj = backtrack_step(gi, gj, m, n)
        end
        chain.reverse!
        to_arc_length(chain)
      end

      # One monotone step towards the origin, from the global
      # vertex-index parameter point (gi, gj).
      def backtrack_step(gi, gj, m, n)
        candidate_cells(gi, gj, m, n) do |ci, cj, x, y|
          left = @lr[ci][cj]
          if left && left[0] <= y + WITNESS_TOL && (x > WITNESS_TOL || left[0] < y - WITNESS_TOL)
            return [ci.to_f, cj + [left[0], y].min]
          end

          bottom = @br[ci][cj]
          if bottom && bottom[0] <= x + WITNESS_TOL && (y > WITNESS_TOL || bottom[0] < x - WITNESS_TOL)
            return [ci + [bottom[0], x].min, cj.to_f]
          end

          corner_reachable = (left && left[0] <= 0.0) || (bottom && bottom[0] <= 0.0)
          return [ci.to_f, cj.to_f] if corner_reachable && (x > 0.0 || y > 0.0)
        end
        raise 'could not reconstruct a monotone witness'
      end

      # Yields each diagram cell containing the global point (gi, gj),
      # as [cell_i, cell_j, local_x, local_y].
      def candidate_cells(gi, gj, m, n)
        is = gi == gi.floor && gi > 0 && gi <= m ? [gi.floor - 1, gi.floor] : [gi.floor]
        js = gj == gj.floor && gj > 0 && gj <= n ? [gj.floor - 1, gj.floor] : [gj.floor]
        is.each do |ci|
          next if ci < 0 || ci >= m

          js.each do |cj|
            next if cj < 0 || cj >= n

            yield ci, cj, gi - ci, gj - cj
          end
        end
      end

      def to_arc_length(chain)
        chain.map do |s, t|
          [arc_length(@cum_a, @lengths_a, s), arc_length(@cum_b, @lengths_b, t)]
        end.each_with_object([]) do |pair, out|
          out << pair unless out.last == pair
        end
      end

      def arc_length(cum, lengths, parameter)
        i = parameter.floor
        return cum[-1] if i >= lengths.size

        cum[i] + (parameter - i) * lengths[i]
      end
    end
  end
end
