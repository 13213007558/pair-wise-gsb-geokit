module Geokit
  # Geodesic polygon on an ideal sphere. Rings are arrays of LatLng-like
  # points and every edge is the shorter great-circle arc.
  #
  # Outer rings are oriented so that the covered region is on the left when
  # the ring is viewed from outside the sphere. Holes use the opposite
  # boundary orientation; their algebraic contribution is therefore
  # subtracted. Reversing an outer ring selects the complementary spherical
  # region, which may be larger than one hemisphere. Rings should be simple;
  # outer rings should be disjoint and holes should lie within an outer ring.
  class SphericalPolygon
    RAD = Math::PI / 180.0
    DEFAULT_RADIUS = 6_376_772.71
    CENTER_EPSILON = 1e-12

    attr_reader :outer_rings, :holes, :radius

    def initialize(outer_rings, holes = [], radius: DEFAULT_RADIUS)
      rings = wrap_outer_rings(outer_rings)
      @outer_rings = rings.map { |ring| normalize_ring(ring) }
      @holes = wrap_rings(holes).map { |ring| normalize_ring(ring) }
      @radius = radius.to_f

      if @outer_rings.empty?
        raise ArgumentError, 'SphericalPolygon requires an outer ring'
      end
      raise ArgumentError, 'radius must be positive' if @radius <= 0

      calculate_geometry
    end

    # Solid angle in steradians. Outer rings add area; holes subtract it.
    def solid_angle
      @solid_angle
    end

    # Surface area in square radius units. The default radius is in meters.
    def area(radius: @radius)
      solid_angle * radius.to_f**2
    end
    alias_method :spherical_area, :area

    # Unnormalized first surface moment, integral of position over the area.
    # Components are in the same Cartesian xyz coordinate system as unit
    # vectors derived from longitude and latitude and have cubic radius units.
    def first_moment(radius: @radius)
      @unit_first_moment.map { |component| component * radius.to_f**3 }
    end
    alias_method :spherical_first_moment, :first_moment

    # Direction of the first surface moment, converted back to lat/lng. This
    # is the radial projection of the surface center of mass and is nil when
    # the moment is too close to zero to identify a unique direction.
    def center(epsilon: CENTER_EPSILON)
      return nil if @unit_first_moment_norm <= epsilon

      direction = @unit_first_moment.map do |component|
        component / @unit_first_moment_norm
      end
      vector_to_latlng(direction)
    end
    alias_method :spherical_center, :center
    alias_method :spherical_centroid, :center
    alias_method :centroid, :center

    private

    def wrap_outer_rings(outer_rings)
      wrap_rings(outer_rings)
    end

    def wrap_rings(rings)
      value = Array(rings)
      return [] if value.empty?
      return value if value.all? { |ring| ring_collection?(ring) }

      [value]
    end

    def ring_collection?(value)
      value.respond_to?(:to_ary) && !coordinate_pair?(value)
    end

    def coordinate_pair?(value)
      value.respond_to?(:to_ary) &&
        value.size == 2 &&
        value.all? { |component| component.is_a?(Numeric) }
    end

    def normalize_ring(points)
      ring = Array(points).map { |point| Geokit::LatLng.normalize(point) }
      raise ArgumentError, 'a ring requires at least three points' if ring.size < 3

      ring = ring.chunk_while(&:==).map(&:first)
      ring.pop if ring.last == ring.first && ring.size > 1
      ring.shift if ring.first == ring.last

      if ring.size < 3
        raise ArgumentError, 'a ring requires at least three distinct points'
      end

      ring << ring.first
    end

    def calculate_geometry
      outer_area = 0.0
      moment = [0.0, 0.0, 0.0]

      @outer_rings.each do |ring|
        ring_area, ring_moment = ring_geometry(ring)
        outer_area += ring_area
        moment = add_vectors(moment, ring_moment)
      end

      @holes.each do |ring|
        left_area, ring_moment = ring_geometry(ring)
        outer_area += left_area - 4.0 * Math::PI
        moment = add_vectors(moment, ring_moment)
      end

      @solid_angle = outer_area
      @unit_first_moment = moment
      @unit_first_moment_norm = vector_norm(moment)
    end

    def ring_geometry(closed_ring)
      ring = closed_ring[0...-1]
      vertices = ring.map { |point| point_to_vector(point) }
      edge_normals = []
      edge_angles = []

      vertices.each_with_index do |start, index|
        finish = vertices[(index + 1) % vertices.size]
        edge_normals << edge_normal(start, finish)
        edge_angles << central_angle(start, finish)
      end

      total_turn = vertices.each_with_index.sum do |vertex, index|
        previous_normal = edge_normals[(index - 1) % vertices.size]
        normal = edge_normals[index]
        signed_angle(previous_normal, normal, vertex)
      end

      left_area = (2.0 * Math::PI - total_turn) % (4.0 * Math::PI)
      moment = edge_normals.each_with_index.each_with_object([0.0, 0.0, 0.0]) do |(normal, index), sum|
        contribution = normal.map { |component| 0.5 * edge_angles[index] * component }
        sum.replace(add_vectors(sum, contribution))
      end

      [left_area, moment]
    end

    def point_to_vector(point)
      latitude = point.lat * RAD
      longitude = point.lng * RAD
      cos_latitude = Math.cos(latitude)

      [
        cos_latitude * Math.cos(longitude),
        cos_latitude * Math.sin(longitude),
        Math.sin(latitude)
      ]
    end

    def vector_to_latlng(vector)
      x, y, z = vector
      latitude = Math.asin(clamp(z / vector_norm(vector), -1.0, 1.0))
      longitude = if Math.cos(latitude).abs < CENTER_EPSILON
                    0.0
                  else
                    Math.atan2(y, x)
                  end
      Geokit::LatLng.new(latitude / RAD, longitude / RAD)
    end

    def edge_normal(start, finish)
      normal = cross_product(start, finish)
      norm = vector_norm(normal)
      return normal.map { |component| component / norm } unless norm < 1e-15

      if dot_product(start, finish) < -1.0 + 1e-12
        raise ArgumentError,
          'antipodal adjacent points do not define a unique great-circle arc'
      end

      raise ArgumentError, 'a ring contains a zero-length edge'
    end

    def central_angle(start, finish)
      cosine = clamp(dot_product(start, finish), -1.0, 1.0)
      angle = Math.acos(cosine)
      return angle unless angle >= Math::PI - 1e-12

      raise ArgumentError,
        'spherical polygon edges must be shorter than one half great circle'
    end

    def signed_angle(from, to, axis)
      Math.atan2(
        dot_product(axis, cross_product(from, to)),
        dot_product(from, to)
      )
    end

    def cross_product(left, right)
      [
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0]
      ]
    end

    def dot_product(left, right)
      left[0] * right[0] + left[1] * right[1] + left[2] * right[2]
    end

    def add_vectors(left, right)
      [left[0] + right[0], left[1] + right[1], left[2] + right[2]]
    end

    def vector_norm(vector)
      Math.sqrt(dot_product(vector, vector))
    end

    def clamp(value, minimum, maximum)
      [[minimum, value].max, maximum].min
    end
  end
end
