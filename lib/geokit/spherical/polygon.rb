module Geokit
  module Spherical
    # A polygon on the sphere whose edges are great-circle short arcs
    # (strictly shorter than half a great circle).
    #
    # Ring orientation convention: the polygon interior is on the LEFT of
    # travel, i.e. outer rings run counterclockwise and holes run
    # clockwise when seen from outside the sphere. Input rings are
    # re-oriented automatically when they enclose less than a hemisphere;
    # each input ring must therefore enclose less than a hemisphere.
    #
    # Rings may be given closed (last point equal to the first) or open;
    # they are stored open. Rotating the start vertex of a ring never
    # changes the geometry.
    class Polygon
      # Arrays of Geokit::LatLng, stored open (no duplicated closing point).
      attr_reader :outer, :holes

      # +outer+ and +holes+ are arrays of Geokit::LatLng or [lat, lng]
      # pairs. +validate:+ is internal: boolean results are constructed
      # already valid and correctly oriented.
      def initialize(outer, holes = [], validate: true)
        if validate
          @outer = normalize_input_ring(outer)
          @holes = holes.map { |ring| normalize_input_ring(ring) }
          validate_and_orient!
        else
          @outer = outer
          @holes = holes
        end
      end

      # All rings as arrays of unit vectors, outer ring first.
      def vector_rings
        @vector_rings ||= ([outer] + holes).map do |ring|
          ring.map { |ll| Vector.from_lat_lng(ll.lat, ll.lng) }
        end
      end

      # Spherical area in steradians (unit sphere).
      def area
        signed = vector_rings.sum { |ring| Core.signed_ring_area(ring) }
        signed %= 4.0 * Math::PI
        signed += 4.0 * Math::PI if signed < 0
        signed
      end

      # True when +point+ (Geokit::LatLng or [lat, lng]) is inside the
      # polygon. Points on the boundary count as contained.
      def contains?(point)
        p = Vector.from_lat_lng(*self.class.coordinates_of(point))
        rings = vector_rings
        return true if Core.on_boundary?(rings, p)
        Core.polygon_winding_number(rings, p) > 0.5
      end

      def union(other)
        Boolean.union(self, other)
      end

      def intersection(other)
        Boolean.intersection(self, other)
      end

      def difference(other)
        Boolean.difference(self, other)
      end

      # Extracts [lat, lng] floats from a Geokit::LatLng or a pair.
      def self.coordinates_of(point)
        if point.respond_to?(:lat) && point.respond_to?(:lng)
          [point.lat.to_f, point.lng.to_f]
        elsif point.is_a?(Array) && point.size == 2
          [point[0].to_f, point[1].to_f]
        else
          raise ArgumentError, "#{point.inspect} cannot be read as lat/lng"
        end
      end

      private

      def normalize_input_ring(ring)
        points = ring.map { |pt| Geokit::LatLng.new(*self.class.coordinates_of(pt)) }
        points.pop while points.size > 1 && points.first == points.last
        deduped = points.each_with_object([]) do |pt, acc|
          acc << pt if acc.empty? || acc.last != pt
        end
        if deduped.size < 3
          raise InvalidPolygonError, "a ring needs at least 3 distinct vertices, got #{deduped.size}"
        end
        deduped
      end

      def validate_and_orient!
        rings = vector_rings
        rings.each { |ring| validate_edges!(ring) }
        validate_no_crossings!(rings)
        orient_rings!
        validate_nesting!(rings)
        @vector_rings = nil # recompute after possible re-orientation
        vector_rings
      end

      def validate_edges!(ring)
        ring.each_index do |i|
          angle = Vector.angle_between(ring[i], ring[(i + 1) % ring.size])
          if angle < Core::ANG_TOL
            raise InvalidPolygonError, "degenerate edge between vertices #{i} and #{i + 1}"
          end
          if (Math::PI - angle).abs < Core::ANG_TOL
            raise AntipodalEdgeError,
              "edge between vertices #{i} and #{i + 1} is (nearly) half a great circle; " +              'its short arc is ambiguous - split it with an intermediate vertex'
          end
        end
      end

      def validate_no_crossings!(rings)
        edges = []
        rings.each_with_index do |ring, ri|
          ring.each_index do |i|
            edges << [ri, i, ring.size, ring[i], ring[(i + 1) % ring.size]]
          end
        end
        edges.combination(2) { |e1, e2| check_edge_pair(e1, e2) }
      end

      def check_edge_pair(e1, e2)
        ri1, i1, n1, a, b = e1
        ri2, i2, n2, c, d = e2
        same_ring = ri1 == ri2
        adjacent = same_ring && ((i1 + 1) % n1 == i2 || (i2 + 1) % n2 == i1)
        type, points = Core.arc_intersections(a, b, c, d)
        return if points.empty?

        problem =
          if type == :coincident
            dedupe_vectors(points).size >= 2 || (same_ring && !adjacent)
          else
            proper = points.any? do |p|
              strictly_interior?(a, b, p) && strictly_interior?(c, d, p)
            end
            proper || (same_ring && !adjacent)
          end
        return unless problem

        if same_ring
          raise SelfIntersectionError,
            "ring #{ri1} crosses or touches itself between edges #{i1} and #{i2}"
        else
          raise RingNestingError, "rings #{ri1} and #{ri2} cross or overlap"
        end
      end

      def strictly_interior?(a, b, p)
        Core.on_arc?(a, b, p) &&
          Vector.angle_between(a, p) > Core::ANG_TOL &&
          Vector.angle_between(p, b) > Core::ANG_TOL
      end

      def dedupe_vectors(points)
        points.each_with_object([]) do |p, acc|
          acc << p unless acc.any? { |q| Vector.angle_between(p, q) < Core::ANG_TOL }
        end
      end

      # Enforce the orientation convention (interior on the left). Rings
      # enclosing more than a hemisphere cannot be told apart from their
      # complement, so input rings must enclose less than a hemisphere.
      def orient_rings!
        rings = vector_rings
        rings.each_with_index do |ring, idx|
          area = Core.signed_ring_area(ring)
          if area.abs < 1e-12
            raise InvalidPolygonError, "ring #{idx} encloses (nearly) zero area"
          end
          want_positive = idx.zero? # outer ring positive, holes negative
          reverse_ring!(idx) if (area.positive? != want_positive)
        end
      end

      def reverse_ring!(idx)
        if idx.zero?
          @outer = @outer.reverse
        else
          @holes[idx - 1] = @holes[idx - 1].reverse
        end
        @vector_rings = nil
      end

      def validate_nesting!(rings)
        outer = rings[0]
        hole_rings = rings[1..] || []
        hole_rings.each_with_index do |hole, hi|
          probe = edge_midpoint(hole)
          if Core.on_boundary?([outer], probe) ||
              Core.ring_winding_number(outer, probe).abs < 0.5
            raise RingNestingError, "hole #{hi} is not strictly inside the outer ring"
          end
          hole_rings.each_with_index do |other, oi|
            next if oi == hi
            if Core.ring_winding_number(other, probe).abs > 0.5
              raise RingNestingError, "hole #{hi} is nested inside hole #{oi}"
            end
          end
        end
      end

      def edge_midpoint(ring)
        Vector.normalize(Vector.add(ring[0], ring[1]))
      end
    end

    # The result of a spherical boolean operation: zero or more polygons,
    # each with one outer ring and any number of holes.
    class MultiPolygon
      attr_reader :polygons

      def initialize(polygons = [])
        @polygons = polygons
      end

      def empty?
        @polygons.empty?
      end

      # Spherical area in steradians (unit sphere).
      def area
        @polygons.sum(&:area)
      end

      def contains?(point)
        @polygons.any? { |polygon| polygon.contains?(point) }
      end

      def union(other)
        Boolean.union(self, other)
      end

      def intersection(other)
        Boolean.intersection(self, other)
      end

      def difference(other)
        Boolean.difference(self, other)
      end
    end
  end
end
