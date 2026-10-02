module Geokit
  class SphericalCap
    RADIUS_EPSILON = 1e-10
    VECTOR_EPSILON = 1e-12
    CENTER_EPSILON = 1e-8
    ANGLE_EPSILON = 1e-9

    Vector = Struct.new(:x, :y, :z)
    Group = Struct.new(:vector, :indices)
    Candidate = Struct.new(:center)

    attr_reader :center, :angular_radius, :radius, :units,
                :support_indices, :support_points,
                :boundary_indices, :boundary_points,
                :nonunique_reason

    def initialize(center, angular_radius, units, points, support_indices,
                   boundary_indices, center_unique, nonunique_reason,
                   center_vector: nil)
      @center = center
      @center_vector = center_vector || self.class.send(:vector, center)
      @angular_radius = angular_radius
      @units = units
      @radius = angular_radius * LatLng.units_sphere_multiplier(units)
      @points = points
      @support_indices = support_indices
      @support_points = support_indices.map { |index| points[index] }
      @boundary_indices = boundary_indices
      @boundary_points = boundary_indices.map { |index| points[index] }
      @center_unique = center_unique
      @nonunique_reason = nonunique_reason
    end

    def angular_radius_degrees
      LatLng.rad2deg(angular_radius)
    end

    def radius_in_units(units)
      units = LatLng.get_units!(units: units)
      angular_radius * LatLng.units_sphere_multiplier(units)
    end

    def contains?(point)
      angular_distance(point) <= angular_radius + ANGLE_EPSILON
    end

    def angular_distance(point)
      point = LatLng.normalize(point)
      self.class.send(:angular_distance, self.class.send(:vector, point), @center_vector)
    end

    def unique_center?
      @center_unique
    end
    alias unique? unique_center?

    def nonunique_center?
      !unique_center?
    end
    alias nonunique? nonunique_center?

    def self.for_points(*points, units: Geokit.default_units)
      if points.one? && points.first.is_a?(Array) && !coordinate_pair?(points.first)
        points = points.first
      end

      raise ArgumentError, 'at least one point is required' if points.empty?

      units = LatLng.get_units!(units: units)
      points = points.map { |point| LatLng.normalize(point) }
      unless points.all?(&:valid?)
        raise ArgumentError, 'all points must have valid latitude and longitude'
      end

      groups = points.each_with_index.map do |point, index|
        Group.new(vector(point), [index])
      end
      groups = merge_groups(groups)
      groups.sort_by! { |group| [group.vector.x, group.vector.y, group.vector.z] }

      candidates = build_candidates(groups)
      evaluated = candidates.map do |candidate|
        minimum_dot = groups.map { |group| dot(candidate.center, group.vector) }.min
        [Math.acos(clamp(minimum_dot)), candidate]
      end
      best_radius = evaluated.map(&:first).min
      best_evaluated = evaluated.select do |radius, _|
        (radius - best_radius).abs <= RADIUS_EPSILON
      end

      center_groups = distinct_centers(best_evaluated)
      center_group = center_groups.min_by do |group|
        representative = group[:representative]
        [representative.x, representative.y, representative.z]
      end
      center_vector = center_group[:representative]

      boundary_indices = points.each_with_index.map do |point, index|
        distance = angular_distance(vector(point), center_vector)
        index if (distance - best_radius).abs <= ANGLE_EPSILON
      end.compact
      boundary_groups = groups.select do |group|
        group.indices.any? { |index| boundary_indices.include?(index) }
      end

      support_groups = support_groups(
        boundary_groups, center_vector, best_radius
      )
      support_indices = support_groups.map { |group| group.indices.min }.sort
      center = point_from_vector(center_vector)
      center_unique, reason = center_uniqueness(
        best_radius, center_vector, groups, boundary_groups, center_groups
      )

      new(
        center, best_radius, units, points, support_indices,
        boundary_indices, center_unique, reason,
        center_vector: center_vector
      )
    end

    def self.coordinate_pair?(value)
      value.is_a?(Array) && value.size == 2 && value.first.is_a?(Numeric)
    end
    private_class_method :coordinate_pair?

    def self.merge_groups(groups)
      groups.each_with_object([]) do |group, merged|
        existing = merged.find do |candidate|
          angular_distance(group.vector, candidate.vector) <= ANGLE_EPSILON
        end
        if existing
          existing.indices.concat(group.indices).sort!
        else
          merged << group
        end
      end
    end
    private_class_method :merge_groups

    def self.build_candidates(groups)
      candidates = groups.map { |group| Candidate.new(group.vector) }

      groups.combination(2).each do |first, second|
        add_pair_candidates(candidates, first, second)
      end

      groups.combination(3).each do |first, second, third|
        normal = normalize(cross(subtract(second.vector, first.vector),
                                 subtract(third.vector, first.vector)))
        add_center(candidates, normal, [first, second, third]) if normal
      end

      candidates
    end
    private_class_method :build_candidates

    def self.add_pair_candidates(candidates, first, second)
      pair = [first, second]

      if dot(first.vector, second.vector) <= -1 + VECTOR_EPSILON
        add_center(candidates, perpendicular(first.vector), pair)
      else
        midpoint = normalize(add(first.vector, second.vector))
        add_center(candidates, midpoint, pair) if midpoint
      end

      normal = normalize(cross(first.vector, second.vector))
      add_center(candidates, normal, pair) if normal
    end
    private_class_method :add_pair_candidates

    def self.add_center(candidates, center, _support_groups)
      candidates << Candidate.new(center)
      candidates << Candidate.new(scale(center, -1))
    end
    private_class_method :add_center

    def self.distinct_centers(evaluated_candidates)
      evaluated_candidates.each_with_object([]) do |(_, candidate), groups|
        group = groups.find do |existing|
          angular_distance(existing[:representative], candidate.center) <= CENTER_EPSILON
        end

        if group
          group[:centers] << candidate.center
        else
          groups << {
            representative: candidate.center,
            centers: [candidate.center]
          }
        end
      end.map do |group|
        group[:representative] = group[:centers].min_by do |center|
          [center.x, center.y, center.z]
        end
        group
      end
    end
    private_class_method :distinct_centers

    def self.support_groups(boundary_groups, center, radius)
      if radius <= ANGLE_EPSILON
        return boundary_groups.sort_by { |group| group.indices.min }.first(1)
      end

      boundary_groups.combination(2).each do |first, second|
        pair_support = pair_support?(
          first.vector, second.vector, center, radius
        )
        return [first, second] if pair_support
      end

      boundary_groups.combination(3).each do |first, second, third|
        normal = normalize(cross(subtract(second.vector, first.vector),
                                 subtract(third.vector, first.vector)))
        if normal &&
           (angular_distance(normal, center) <= CENTER_EPSILON ||
            angular_distance(scale(normal, -1), center) <= CENTER_EPSILON)
          return [first, second, third]
        end
      end

      boundary_groups.sort_by { |group| group.indices.min }.first(1)
    end
    private_class_method :support_groups

    def self.pair_support?(first, second, center, radius)
      if dot(first, second) <= -1 + VECTOR_EPSILON
        return near_hemisphere?(radius) &&
               dot(first, center).abs <= CENTER_EPSILON
      end

      midpoint = normalize(add(first, second))
      return true if midpoint && angular_distance(midpoint, center) <= CENTER_EPSILON

      return false unless near_hemisphere?(radius)

      normal = normalize(cross(first, second))
      normal && angular_distance(normal, center) <= CENTER_EPSILON
    end
    private_class_method :pair_support?

    def self.center_uniqueness(radius, center, groups, boundary_groups, center_groups)
      if (radius - Math::PI / 2).abs <= RADIUS_EPSILON
        rank = vector_rank(boundary_groups.map(&:vector))
        has_interior_point = groups.any? { |group|
          dot(center, group.vector) > CENTER_EPSILON
        }
        return [true, nil] if rank == 2 && has_interior_point

        reason = if rank == 2
                   :opposite_hemisphere_centers
                 else
                   :hemisphere_center_continuum
                 end
        return [false, reason]
      end

      return [true, nil] if radius < Math::PI / 2

      [center_groups.one?, center_groups.one? ? nil : :multiple_spherical_caps]
    end
    private_class_method :center_uniqueness

    def self.near_hemisphere?(radius)
      (radius - Math::PI / 2).abs <= RADIUS_EPSILON
    end
    private_class_method :near_hemisphere?

    def self.vector_rank(vectors)
      first = vectors.find { |vector| norm(vector) > VECTOR_EPSILON }
      return 0 unless first

      second = vectors.find do |vector|
        norm(cross(first, vector)) > CENTER_EPSILON
      end
      return 1 unless second

      normal = cross(first, second)
      third = vectors.find { |vector| dot(normal, vector).abs > CENTER_EPSILON }
      third ? 3 : 2
    end
    private_class_method :vector_rank

    def self.perpendicular(vector)
      axes = [
        Vector.new(1, 0, 0), Vector.new(0, 1, 0), Vector.new(0, 0, 1)
      ]
      axis = axes.max_by { |candidate| (dot(vector, candidate).abs - 1).abs }
      normalize(cross(vector, axis))
    end
    private_class_method :perpendicular

    def self.vector(point)
      latitude = LatLng.deg2rad(point.lat)
      longitude = LatLng.deg2rad(point.lng)
      cosine_latitude = Math.cos(latitude)
      Vector.new(
        cosine_latitude * Math.cos(longitude),
        cosine_latitude * Math.sin(longitude),
        Math.sin(latitude)
      )
    end
    private_class_method :vector

    def self.point_from_vector(vector)
      latitude = Math.asin(clamp(vector.z))
      longitude = Math.atan2(vector.y, vector.x)
      if (latitude.abs - Math::PI / 2).abs <= VECTOR_EPSILON
        longitude = 0
      end
      LatLng.new(LatLng.rad2deg(latitude), LatLng.rad2deg(longitude))
    end
    private_class_method :point_from_vector

    def self.add(first, second)
      Vector.new(first.x + second.x, first.y + second.y, first.z + second.z)
    end
    private_class_method :add

    def self.subtract(first, second)
      Vector.new(first.x - second.x, first.y - second.y, first.z - second.z)
    end
    private_class_method :subtract

    def self.scale(vector, factor)
      Vector.new(vector.x * factor, vector.y * factor, vector.z * factor)
    end
    private_class_method :scale

    def self.dot(first, second)
      first.x * second.x + first.y * second.y + first.z * second.z
    end
    private_class_method :dot

    def self.cross(first, second)
      Vector.new(
        first.y * second.z - first.z * second.y,
        first.z * second.x - first.x * second.z,
        first.x * second.y - first.y * second.x
      )
    end
    private_class_method :cross

    def self.norm(vector)
      Math.sqrt(dot(vector, vector))
    end
    private_class_method :norm

    def self.normalize(vector)
      magnitude = norm(vector)
      return if magnitude <= VECTOR_EPSILON
      scale(vector, 1 / magnitude)
    end
    private_class_method :normalize

    def self.angular_distance(first, second)
      Math.acos(clamp(dot(first, second)))
    end
    private_class_method :angular_distance

    def self.clamp(value)
      [[value, -1.0].max, 1.0].min
    end
    private_class_method :clamp
  end
end
