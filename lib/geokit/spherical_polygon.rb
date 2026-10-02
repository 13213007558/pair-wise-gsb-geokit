module Geokit
  # Spherical geodesic polygons made from one or more outer rings and holes.
  # Rings are arrays of Geokit::LatLng-compatible points. Each pair of
  # adjacent vertices is connected by the shorter great-circle arc, so an
  # antipodal pair is invalid because its arc is not unique.
  #
  # Ring orientation follows the boundary, not the visual winding in a
  # longitude/latitude chart: the represented region must be on the left of
  # every directed edge when viewed from outside the sphere. A ring smaller
  # than a hemisphere therefore has counter-clockwise outer winding. A small
  # hole uses the opposite winding. Reversing a small outer ring represents
  # the complementary region larger than a hemisphere.
  class SphericalPolygon
    UNIT_RADIUS = 1.0
    private_constant :UNIT_RADIUS

    # Minimal three-dimensional vector used for unit-sphere moments.
    class Vector
      attr_reader :x, :y, :z

      def self.[](x, y, z); new(x, y, z); end

      def initialize(x, y, z)
        @x = x.to_f
        @y = y.to_f
        @z = z.to_f
      end

      def +(other)
        self.class.new(x + other.x, y + other.y, z + other.z)
      end

      def -(other); self.class.new(x - other.x, y - other.y, z - other.z); end

      def -@; self.class.new(-x, -y, -z); end

      def *(scale); self.class.new(x * scale, y * scale, z * scale); end

      def /(scale); self.class.new(x / scale, y / scale, z / scale); end

      def dot(other)
        x * other.x + y * other.y + z * other.z
      end

      def cross(other)
        self.class.new(y * other.z - z * other.y,
                       z * other.x - x * other.z,
                       x * other.y - y * other.x)
      end

      def norm
        Math.sqrt(dot(self))
      end

      def normalize
        self / norm
      end
      def to_a
        [x, y, z]
      end

      def ==(other)
        other.is_a?(self.class) && to_a == other.to_a
      end

      alias_method :eql?, :==

      def inspect
        "#<#{self.class.name} #{x.inspect}, #{y.inspect}, #{z.inspect}>"
      end
    end

    attr_reader :outer_rings, :holes

    def initialize(outer_rings = [], holes: [], whole: false)
      if whole && (!Array(outer_rings).empty? || !Array(holes).empty?)
        raise ArgumentError, 'the whole sphere cannot have boundary rings'
      end

      @whole = whole
      @outer_rings = wrap_rings(outer_rings).map { |ring| normalize_ring(ring) }
      @holes = wrap_rings(holes).map { |ring| normalize_ring(ring) }
      calculate_geometry
    end

    def self.whole_sphere
      new([], whole: true)
    end

    # Solid angle in steradians on the unit sphere.
    def solid_angle
      @solid_angle
    end

    # Spherical surface area. Without +radius+ or +units+, the result is the
    # solid angle in steradians. With units, it is square distance units.
    def area(radius: nil, units: nil)
      @solid_angle * sphere_radius(radius: radius, units: units)**2
    end

    # Unnormalized first moment integral over a sphere of uniform surface
    # density. The default result uses the unit sphere. With a radius or units,
    # positions and area are scaled, giving distance units cubed.
    def first_moment(radius: nil, units: nil)
      @first_moment * sphere_radius(radius: radius, units: units)**3
    end

    # The spherical center is the normalized first-moment direction. It is
    # independent of the sphere radius. When opposing surface balances the
    # moment relative to the rings' scale, no unique center exists and nil is
    # returned.
    def center(epsilon: 1e-12)
      scale = [@moment_scale, @solid_angle.abs].max
      return nil if @first_moment.norm <= epsilon * scale

      direction = @first_moment.normalize
      latitude = Geokit::LatLng.rad2deg(Math.asin(clamp(direction.z, -1.0, 1.0)))
      longitude = Geokit::LatLng.rad2deg(Math.atan2(direction.y, direction.x))
      Geokit::LatLng.new(latitude, longitude)
    end

    alias_method :spherical_centroid, :center

    private

    def wrap_rings(value)
      rings = Array(value)
      return rings if rings.empty?
      point_like?(rings.first) ? [rings] : rings
    end

    def point_like?(value)
      return true unless value.is_a?(Array)
      return false if value.empty? || value.first.is_a?(Array)

      value.length == 2 && value.all? { |component| component.is_a?(Numeric) }
    end

    def normalize_ring(points)
      ring = Array(points).map { |point| Geokit::LatLng.normalize(point) }
      raise ArgumentError, 'a ring requires at least three points' if ring.length < 3

      ring.pop if ring.length > 1 && ring.last == ring.first
      ring = ring.each_with_object([]) do |point, unique|
        unique << point unless unique.last == point
      end
      raise ArgumentError, 'a ring requires at least three unique points' if ring.length < 3

      ring
    end

    def calculate_geometry
      @solid_angle = @whole ? 4.0 * Math::PI : 0.0
      @first_moment = Vector.new(0.0, 0.0, 0.0)
      @moment_scale = @solid_angle

      @outer_rings.each do |ring|
        signed_angle, moment, path_length = ring_geometry(ring)
        @moment_scale += path_length
        raise ArgumentError, 'a degenerate ring has zero solid angle' if signed_angle.zero?
        @solid_angle += signed_angle.positive? ? signed_angle : signed_angle + 4.0 * Math::PI
        @first_moment += moment
      end

      @holes.each do |ring|
        signed_angle, moment, path_length = ring_geometry(ring)
        @moment_scale += path_length
        raise ArgumentError, 'a degenerate ring has zero solid angle' if signed_angle.zero?
        @solid_angle += signed_angle.negative? ? signed_angle : signed_angle - 4.0 * Math::PI
        @first_moment += moment
      end
    end

    def ring_geometry(ring)
      vertices = ring.map { |point| unit_vector(point) }
      moment = Vector.new(0.0, 0.0, 0.0)
      path_length = 0.0

      vertices.each_with_index do |start, index|
        finish = vertices[(index + 1) % vertices.length]
        edge = edge_moment(start, finish)
        moment += edge
        path_length += edge.norm * 2.0
      end

      [signed_solid_angle(vertices), moment, path_length]
    end

    def edge_moment(start, finish)
      normal = start.cross(finish)
      cross_norm = normal.norm
      dot_product = clamp(start.dot(finish), -1.0, 1.0)
      angle = Math.atan2(cross_norm, dot_product)

      if angle.zero?
        Vector.new(0.0, 0.0, 0.0)
      else
        if dot_product.negative? && cross_norm <= 1e-15
          raise ArgumentError,
            'antipodal ring vertices do not define a unique great-circle arc'
        end
        normal * (angle / (2.0 * cross_norm))
      end
    end

    def signed_solid_angle(vertices)
      reference = reference_vector(vertices)
      angle = vertices.each_with_index.sum do |start, index|
        finish = vertices[(index + 1) % vertices.length]
        numerator = reference.dot(start.cross(finish))
        denominator = 1.0 + reference.dot(start) + start.dot(finish) +
                      finish.dot(reference)
        2.0 * Math.atan2(numerator, denominator)
      end

      two_pi = 2.0 * Math::PI
      four_pi = 4.0 * Math::PI
      wrapped = angle - four_pi * ((angle + two_pi) / four_pi).floor
      if (wrapped.abs - two_pi).abs <= 1e-12
        two_pi * (angle >= 0.0 ? 1.0 : -1.0)
      else
        wrapped
      end
    end

    def reference_vector(vertices)
      edges = vertices.each_with_index.map do |start, index|
        [start, vertices[(index + 1) % vertices.length]]
      end

      1000.times do |index|
        longitude = index * 2.399963229728653
        z = 1.0 - 2.0 * ((index * 0.6180339887498949 + 0.5) % 1.0)
        horizontal = Math.sqrt([1.0 - z * z, 0.0].max)
        candidate = Vector.new(
          horizontal * Math.cos(longitude),
          horizontal * Math.sin(longitude),
          z
        )

        next if vertices.any? { |vertex| candidate.dot(vertex).abs > 0.999999999999 }
        next if edges.any? do |start, finish|
          edge_normal = start.cross(finish)
          edge_norm = edge_normal.norm
          edge_norm > 1e-15 && (candidate.dot(edge_normal) / edge_norm).abs < 1e-10
        end
        return candidate
      end

      raise ArgumentError, 'could not choose a stable reference direction for the ring'
    end

    def unit_vector(point)
      latitude = Geokit::LatLng.deg2rad(point.lat)
      longitude = Geokit::LatLng.deg2rad(point.lng)
      cosine_latitude = Math.cos(latitude)
      Vector.new(
        cosine_latitude * Math.cos(longitude),
        cosine_latitude * Math.sin(longitude),
        Math.sin(latitude)
      )
    end

    def sphere_radius(radius:, units:)
      return radius.to_f if radius
      return UNIT_RADIUS unless units
      Geokit::LatLng.units_sphere_multiplier(units).to_f
    end

    def clamp(value, minimum, maximum)
      return minimum if value < minimum
      return maximum if value > maximum
      value
    end
  end
end
