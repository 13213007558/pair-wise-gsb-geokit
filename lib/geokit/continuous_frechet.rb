module Geokit
  # Computes the continuous Fréchet distance between two projected polylines
  # using a free-space decision procedure and bisection. The application
  # supplies a pointwise local projection so both trajectories are represented
  # in the same planar, metric coordinate system. Distances are returned in
  # the projection's units (typically meters).
  #
  # The result includes numeric lower and upper bounds and a monotonic
  # continuous witness. Each witness pair contains arc lengths in the two
  # projected trajectories; either coordinate may stay fixed, but neither may
  # decrease.
  class ContinuousFrechet
    EPSILON = 1e-12
    MACHINE_EPSILON = Float::EPSILON

    Result = Struct.new(:distance, :lower_bound, :upper_bound, :absolute_error,
                        :exact, :iterations, :witness)

    Witness = Struct.new(:path_a, :path_b, :arc_length_a, :arc_length_b,
                         :pairs) do
      def monotonic?
        pairs.each_cons(2).all? do |first, second|
          second[0] >= first[0] && second[1] >= first[1]
        end
      end
    end

    PlaneCurve = Struct.new(:points, :lengths, :total_length)

    Reach = Struct.new(:border, :row, :column, :lo, :hi, :previous,
                       :entry_s, :entry_t)
    Entrance = Struct.new(:reach, :s, :t, :local)

    def self.distance_between(path_a, path_b, projection = nil, options = {},
                              &block_projection)
      projection ||= block_projection
      points_a = Array(path_a).map { |point| LatLng.normalize(point) }
      points_b = Array(path_b).map { |point| LatLng.normalize(point) }

      return empty_result(points_a, points_b) if
        points_a.empty? || points_b.empty?

      unless projection.respond_to?(:call)
        raise ArgumentError,
              'A callable local metric projection is required.'
      end

      curve_a = build_curve(points_a, projection)
      curve_b = build_curve(points_b, projection)
      tolerance = options.fetch(:absolute_tolerance, 1e-6).to_f
      relative_tolerance = options.fetch(:relative_tolerance, 1e-12).to_f
      iterations_limit = options.fetch(:iterations, 64).to_i
      validate_options(tolerance, relative_tolerance, iterations_limit)

      if curve_a.total_length.zero? || curve_b.total_length.zero?
        return constant_curve_result(points_a, points_b, curve_a, curve_b)
      end

      exact_witness = decide(curve_a, curve_b, 0.0, 0.0)
      if exact_witness
        witness = build_witness(points_a, points_b, curve_a, curve_b,
                                exact_witness)
        return Result.new(0.0, 0.0, 0.0, 0.0, true, 0, witness)
      end

      upper = initial_upper_bound(curve_a, curve_b)
      best_witness = decide(curve_a, curve_b, upper, numeric_slack(curve_a,
                                                                   curve_b))
      until best_witness
        upper *= 2.0
        best_witness = decide(curve_a, curve_b, upper,
                              numeric_slack(curve_a, curve_b))
      end

      lower = 0.0
      iterations = 0
      slack = numeric_slack(curve_a, curve_b)
      iterations_limit.times do
        iterations += 1
        gap = upper - lower
        target_gap = [
          0.0,
          tolerance + relative_tolerance * upper.abs - 2.0 * slack,
        ].max
        break if gap <= target_gap

        middle = lower + gap / 2.0
        candidate_witness = decide(curve_a, curve_b, middle,
                                   slack)
        if candidate_witness
          upper = middle
          best_witness = candidate_witness
        else
          lower = middle
        end
      end

      lower_bound = [0.0, lower - slack].max
      upper_bound = upper + slack
      witness = build_witness(points_a, points_b, curve_a, curve_b,
                              best_witness)
      Result.new(upper, lower_bound, upper_bound,
                 upper_bound - lower_bound, false, iterations, witness)
    end

    def self.empty_result(points_a, points_b)
      if points_a.empty? && points_b.empty?
        witness = Witness.new(points_a, points_b, 0.0, 0.0, [])
        Result.new(0.0, 0.0, 0.0, 0.0, true, 0, witness)
      else
        length_a = points_a.empty? ? 0.0 : nil
        length_b = points_b.empty? ? 0.0 : nil
        witness = Witness.new(points_a, points_b, length_a, length_b, [])
        Result.new(Float::INFINITY, Float::INFINITY, Float::INFINITY, 0.0,
                   true, 0, witness)
      end
    end
    private_class_method :empty_result

    def self.validate_options(tolerance, relative_tolerance, iterations)
      unless tolerance >= 0 && relative_tolerance >= 0 && iterations > 0
        raise ArgumentError, 'Invalid Fréchet tolerance or iteration options.'
      end
    end
    private_class_method :validate_options

    def self.build_curve(points, projection)
      projected = points.map { |point| project_point(point, projection) }
      lengths = [0.0]
      total_length = 0.0

      projected.each_cons(2) do |first, second|
        total_length += distance(first, second)
        lengths << total_length
      end

      PlaneCurve.new(projected, lengths, total_length)
    end
    private_class_method :build_curve

    def self.project_point(point, projection)
      projected = projection.call(point)
      unless projected.is_a?(Array) && projected.size == 2 &&
             projected.all? { |coordinate| coordinate.is_a?(Numeric) &&
                                          coordinate.finite? }
        raise ArgumentError,
              'Projection must return finite [x, y] metric coordinates.'
      end
      projected.map(&:to_f)
    end
    private_class_method :project_point

    def self.initial_upper_bound(curve_a, curve_b)
      upper = distance(curve_a.points.first, curve_b.points.first)
      curve_a.points.each do |point_a|
        curve_b.points.each do |point_b|
          upper = [upper, distance(point_a, point_b)].max
        end
      end
      [upper, EPSILON].max
    end
    private_class_method :initial_upper_bound

    def self.numeric_slack(curve_a, curve_b)
      scale = [1.0, curve_a.total_length, curve_b.total_length].max
      MACHINE_EPSILON * scale
    end
    private_class_method :numeric_slack

    def self.distance(first, second)
      Math.hypot(first[0] - second[0], first[1] - second[1])
    end
    private_class_method :distance

    def self.constant_curve_result(points_a, points_b, curve_a, curve_b)
      maximum_distance = 0.0
      curve_a.points.each do |point_a|
        curve_b.points.each do |point_b|
          maximum_distance =
            [maximum_distance, distance(point_a, point_b)].max
        end
      end

      pairs = if curve_a.total_length.zero? && curve_b.total_length.zero?
                [[0.0, 0.0]]
              elsif curve_a.total_length.zero?
                vertex_pairs(false, curve_b.points.size)
              else
                vertex_pairs(true, curve_a.points.size)
              end
      witness = build_witness(points_a, points_b, curve_a, curve_b, pairs)
      Result.new(maximum_distance, maximum_distance, maximum_distance, 0.0,
                 true, 0, witness)
    end
    private_class_method :constant_curve_result

    def self.vertex_pairs(constant_first, count)
      parameters = (0...count).map(&:to_f)
      if constant_first
        parameters.map { |parameter| [parameter, 0.0] }
      else
        parameters.map { |parameter| [0.0, parameter] }
      end
    end
    private_class_method :vertex_pairs

    def self.decide(curve_a, curve_b, radius, slack)
      effective_radius = radius + slack
      rows = curve_a.points.size - 1
      columns = curve_b.points.size - 1
      origin = Reach.new(:origin, -1, -1, 0.0, 0.0, nil, 0.0, 0.0)
      left_in = Array.new(rows + 1) { Array.new(columns) }
      bottom_in = Array.new(rows) { Array.new(columns + 1) }

      rows.times do |row|
        columns.times do |column|
          left_reach = if row.zero? && column.zero?
                         origin
                       elsif row.zero?
                         nil
                       else
                         left_in[row][column]
                       end
          bottom_reach = if row.zero? && column.zero?
                           origin
                         elsif column.zero?
                           nil
                         else
                           bottom_in[row][column]
                         end

          left_interval = free_vertical(curve_a, curve_b, row, column,
                                        effective_radius)
          bottom_interval = free_horizontal(curve_a, curve_b, row, column,
                                            effective_radius)
          right_interval = free_vertical(curve_a, curve_b, row + 1, column,
                                         effective_radius)
          top_interval = free_horizontal(curve_a, curve_b, row, column,
                                         effective_radius, true)

          left = incoming_from_vertical(left_reach, left_interval, row,
                                        column)
          bottom = incoming_from_horizontal(bottom_reach, bottom_interval,
                                            row, column)

          right = nil
          top = nil

          if bottom && right_interval
            right = make_reach(:right, row, column, right_interval[0],
                               right_interval[1], bottom[:reach],
                               bottom[:s], bottom[:t])
          elsif left && right_interval
            right_lo = [right_interval[0], left[:local]].max
            if right_lo <= right_interval[1] + EPSILON
              right = make_reach(:right, row, column, right_lo,
                                 right_interval[1], left[:reach],
                                 left[:s], left[:t])
            end
          end

          if left && top_interval
            top = make_reach(:top, row, column, top_interval[0],
                             top_interval[1], left[:reach],
                             left[:s], left[:t])
          end
          if top.nil? && bottom && top_interval
            top_lo = [top_interval[0], bottom[:local]].max
            if top_lo <= top_interval[1] + EPSILON
              top = make_reach(:top, row, column, top_lo,
                               top_interval[1], bottom[:reach],
                               bottom[:s], bottom[:t])
            end
          end

          left_in[row + 1][column] = right
          bottom_in[row][column + 1] = top
        end
      end

      right = left_in[rows][columns - 1]
      top = bottom_in[rows - 1][columns]
      return reconstruct(right, rows.to_f, columns.to_f) if
        right && right.hi >= 1.0 - EPSILON
      return reconstruct(top, rows.to_f, columns.to_f) if
        top && top.hi >= 1.0 - EPSILON

      nil
    end
    private_class_method :decide

    def self.incoming_from_vertical(reach, free_interval, row, column)
      incoming(reach, free_interval, row, column, true)
    end
    private_class_method :incoming_from_vertical

    def self.incoming_from_horizontal(reach, free_interval, row, column)
      incoming(reach, free_interval, row, column, false)
    end
    private_class_method :incoming_from_horizontal

    def self.incoming(reach, free_interval, row, column, vertical)
      return nil unless reach && free_interval

      local = [reach.lo, free_interval[0]].max
      high = [reach.hi, free_interval[1]].min
      return nil if local > high + EPSILON

      local = soft_clamp(local)
      s = vertical ? row.to_f : row + local
      t = vertical ? column + local : column.to_f
      Entrance.new(reach, s, t, local)
    end
    private_class_method :incoming

    def self.make_reach(border, row, column, lo, hi, previous, entry_s,
                        entry_t)
      Reach.new(border, row, column, soft_clamp(lo), soft_clamp(hi), previous,
                entry_s.to_f, entry_t.to_f)
    end
    private_class_method :make_reach

    def self.free_vertical(curve_a, curve_b, row, column, radius)
      free_interval(curve_a.points[row], curve_b.points[column],
                    curve_b.points[column + 1], radius)
    end
    private_class_method :free_vertical

    def self.free_horizontal(curve_a, curve_b, row, column, radius,
                             top = false)
      fixed_index = column + (top ? 1 : 0)
      free_interval(curve_b.points[fixed_index], curve_a.points[row],
                    curve_a.points[row + 1], radius)
    end
    private_class_method :free_horizontal

    def self.free_interval(fixed_point, segment_start, segment_end, radius)
      dx = segment_end[0] - segment_start[0]
      dy = segment_end[1] - segment_start[1]
      delta_x = fixed_point[0] - segment_start[0]
      delta_y = fixed_point[1] - segment_start[1]
      a = dx * dx + dy * dy
      b = -2.0 * (delta_x * dx + delta_y * dy)
      c = delta_x * delta_x + delta_y * delta_y - radius * radius

      a = 0.0 if a.abs <= EPSILON

      return [0.0, 1.0] if a.zero? && c <= EPSILON
      return nil if a.zero?

      discriminant = b * b - 4.0 * a * c
      distance_scale = [1.0, Math.sqrt(a), radius.abs].max
      squared_tolerance = if radius.zero?
                            0.0
                          else
                            MACHINE_EPSILON * distance_scale**2
                          end
      return nil if discriminant < -squared_tolerance

      discriminant = 0.0 if discriminant.abs <= squared_tolerance
      root = Math.sqrt(discriminant)
      low = (-b - root) / (2.0 * a)
      high = (-b + root) / (2.0 * a)
      low = [0.0, low].max
      high = [1.0, high].min
      high >= low - EPSILON ? [soft_clamp(low), soft_clamp(high)] : nil
    end
    private_class_method :free_interval

    def self.soft_clamp(value)
      return 0.0 if value <= EPSILON
      return 1.0 if value >= 1.0 - EPSILON
      value
    end
    private_class_method :soft_clamp

    def self.reconstruct(reach, final_s, final_t)
      pairs = []
      current = reach
      s = final_s
      t = final_t

      loop do
        pairs << [s, t]
        break if current.border == :origin

        unless [s, t] == [current.entry_s, current.entry_t]
          pairs << [current.entry_s, current.entry_t]
        end
        s = current.entry_s
        t = current.entry_t
        current = current.previous
      end

      pairs.reverse
    end
    private_class_method :reconstruct

    def self.build_witness(points_a, points_b, curve_a, curve_b,
                           parameter_pairs)
      arc_pairs = parameter_pairs.map do |s, t|
        [arc_position(curve_a, s), arc_position(curve_b, t)]
      end
      Witness.new(points_a, points_b, curve_a.total_length,
                  curve_b.total_length, arc_pairs)
    end
    private_class_method :build_witness

    def self.arc_position(curve, parameter)
      index = [[parameter.floor, 0].max, curve.lengths.size - 1].min
      index = curve.lengths.size - 1 if index == curve.lengths.size - 1
      start_length = curve.lengths[index]
      segment_length = if index < curve.lengths.size - 1
                         curve.lengths[index + 1] - start_length
                       else
                         0.0
                       end
      local = soft_clamp(parameter - index)
      start_length + segment_length * local
    end
    private_class_method :arc_position
  end
end
