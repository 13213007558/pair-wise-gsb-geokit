module Geokit
  # Immutable result of a minimum enclosing spherical cap calculation.
  class SphericalCapResult
    attr_reader :center, :angular_radius, :radius, :units,
                :support_ids, :support_points, :unique
    alias_method :unique?, :unique

    def initialize(center:, angular_radius:, radius:, units:, support_ids:,
                   support_points:, unique:)
      @center = center
      @angular_radius = angular_radius
      @radius = radius
      @units = units
      @support_ids = support_ids
      @support_points = support_points
      @unique = unique
      freeze
    end

    def radius_degrees
      LatLng.rad2deg(angular_radius)
    end

    # Returns true when the optimal center is unique.
    def center_unique?
      unique?
    end
  end

  # Calculates the smallest spherical cap containing a set of points.
  class SphericalCap
    include SphericalCapVectors
    include SphericalCapHull

    EPS = 1e-12
    HEMISPHERE = Math::PI / 2.0

    def self.minimum(points, options = {})
      new(points, options).send(:solve)
    end

    def self.minimum_enclosing(points, options = {})
      minimum(points, options)
    end

    private_class_method :new
    private

    def initialize(points, options)
      @units = LatLng.send(:get_units!, options)
      @input = Array(points)
      raise ArgumentError, 'expected at least one point' if @input.empty?

      labels = options[:ids] || options[:labels]
      if labels && labels.length != @input.length
        raise ArgumentError, 'ids must contain one entry per point'
      end

      @records = @input.each_with_index.map do |point, index|
        location = LatLng.normalize(point)
        validate_coordinates!(location)
        identity = location_identity(location)
        { id: labels ? labels[index] : (identity.nil? ? index : identity),
          point: point,
          location: location,
          vector: to_vector(location) }
      end
    end

    def solve
      unique_points = dedupe(@records)
      best = initial_candidate(unique_points)
      enumerate_pairs(unique_points) { |candidate| best = choose(best, candidate) }
      enumerate_triangle_planes(unique_points) do |plane|
        best = choose(best, candidate(plane))
        best = choose(best, candidate(scale(plane, -1.0)))
      end

      vector = best[:vector]
      minimum_dot = @records.map { |record| dot(vector, record[:vector]) }.min
      radius_angle = acos_clamped(minimum_dot)
      support = @records.select do |record|
        dot(vector, record[:vector]) <= minimum_dot + 1e-10
      end

      SphericalCapResult.new(
        center: to_lat_lng(vector),
        angular_radius: radius_angle,
        radius: radius_angle * LatLng.units_sphere_multiplier(@units),
        units: @units,
        support_ids: support.map { |record| record[:id] },
        support_points: support.map { |record| record[:point] },
        unique: unique_solution?(best, vector, minimum_dot, radius_angle)
      )
    end

    def initial_candidate(unique_points)
      if unique_points.length == 2 && antipodal?(unique_points[0], unique_points[1])
        return candidate(perpendicular(canonical_direction(unique_points[0][:vector])))
      end
      unique_points.map { |record| candidate(record[:vector]) }.reduce { |best, challenger|
        choose(best, challenger)
      }
    end

    def enumerate_pairs(points)
      points.combination(2) do |first, second|
        next if antipodal?(first, second)
        total = add(first[:vector], second[:vector])
        yield candidate(normalize(total))
        yield candidate(normalize(scale(total, -1.0)))
      end
    end

    def candidate(vector)
      minimum_dot = @records.map { |record| dot(vector, record[:vector]) }.min
      { vector: vector, minimum_dot: minimum_dot, alternate: false }
    end

    def choose(best, challenger)
      if challenger[:minimum_dot] > best[:minimum_dot] + EPS
        challenger
      elsif (challenger[:minimum_dot] - best[:minimum_dot]).abs <= EPS
        tied_candidate(best, challenger)
      else
        best
      end
    end

    def tied_candidate(best, challenger)
      return best if same_vector?(best[:vector], challenger[:vector])
      tied = best.dup
      tied[:alternate] = true
      tied[:vector] = challenger[:vector] if
        lex_smaller?(challenger[:vector], best[:vector])
      tied
    end

    def unique_solution?(best, center, minimum_dot, angular_radius)
      return true if minimum_dot > EPS
      return hemisphere_unique?(center, minimum_dot) if
        (angular_radius - HEMISPHERE).abs <= EPS
      !best[:alternate]
    end

    def hemisphere_unique?(center, minimum_dot)
      return false if minimum_dot > EPS || antipodal_center_feasible?(center)

      active = @records.map { |record| record[:vector] }
                       .select { |point| dot(center, point) <= 1e-10 }
      east = normalize(cross(center, [0.0, 0.0, 1.0])) rescue nil
      east = normalize(cross(center, [0.0, 1.0, 0.0])) if east.nil?
      north = normalize(cross(east, center))
      angles = active.map do |point|
        tangent = subtract(point, scale(center, dot(center, point)))
        Math.atan2(dot(north, tangent), dot(east, tangent))
      end.sort
      gaps = angles.each_cons(2).map { |left, right| right - left }
      gaps << angles.first + (2.0 * Math::PI) - angles.last
      gaps.max < Math::PI - 1e-9
    end

    def antipodal_center_feasible?(center)
      @records.all? { |record| dot(scale(center, -1.0), record[:vector]) >= -1e-10 }
    end

    def dedupe(records)
      records.each_with_object([]) do |record, unique|
        duplicate = unique.any? do |existing|
          norm(subtract(record[:vector], existing[:vector])) < EPS
        end
        unique << record unless duplicate
      end
    end

    def antipodal?(first, second)
      (dot(first[:vector], second[:vector]) + 1.0).abs <= EPS
    end

    def perpendicular(vector)
      axis = [0.0, 0.0, 1.0]
      return [1.0, 0.0, 0.0] if dot(vector, axis).abs > 1.0 - EPS
      normalize(cross(vector, axis))
    end

    def canonical_direction(vector)
      opposite = scale(vector, -1.0)
      lex_smaller?(opposite, vector) ? opposite : vector
    end

    def location_identity(location)
      return location.identity if location.respond_to?(:identity)
      return location.id if location.respond_to?(:id)
      nil
    end

    def validate_coordinates!(location)
      return if finite?(location.lat) && finite?(location.lng) &&
                location.lat.between?(-90.0, 90.0)
      raise ArgumentError, "invalid location: #{location}"
    end

    def finite?(value)
      value.is_a?(Numeric) && value.to_f.finite?
    end
  end
end
