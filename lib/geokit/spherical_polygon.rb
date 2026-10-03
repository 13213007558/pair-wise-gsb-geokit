require 'geokit/spherical/vector_math'
require 'geokit/spherical/ring'
require 'geokit/spherical/boolean_engine'

module Geokit
  # A polygon on the sphere, made of one or more parts; each part has one
  # outer ring (shell) and zero or more holes. Edges are great-circle minor
  # arcs: every edge is strictly shorter than a half circumference, so it is
  # the unique shortest great-circle arc between its endpoints.
  #
  # All computations use local coordinates only (unit vectors in 3D); no
  # planar projection is involved, which makes polygons crossing the
  # antimeridian or surrounding a pole behave exactly like any other.
  #
  # Conventions
  # -----------
  # * Ring orientation does not matter on input: every ring is normalized so
  #   that the polygon interior is on the LEFT of the traversal direction
  #   (shells counter-clockwise seen from outside, holes clockwise). Output
  #   rings follow the same convention.
  # * Each ring interior must be smaller than a hemisphere (2*pi steradians).
  # * Rings are open: the last vertex is not repeated (a repeated closing
  #   vertex is dropped on input).
  # * Rotating the start vertex of an input ring does not change the
  #   resulting geometry.
  #
  # Invalid input is reported, never silently repaired:
  # * self-intersecting or self-touching rings,
  # * (nearly) antipodal consecutive vertices / half-circumference edges,
  #   whose great circle is ambiguous,
  # * degenerate (zero-length or near-zero-length) edges,
  # * holes outside their shell, overlapping shells, rings that cross or
  #   touch each other within one polygon.
  #
  # Usage:
  #   poly = Geokit::SphericalPolygon.new([
  #     { shell: [Geokit::LatLng.new(...), ...], holes: [[...]] }
  #   ])
  #   poly.contains?(Geokit::LatLng.new(0, 180))
  #   poly.area # steradians
  #   poly.union(other) / poly.intersection(other) / poly.difference(other)
  #   (also available as |, & and -)
  class SphericalPolygon
    # Base error for spherical polygon problems.
    class Error < StandardError; end
    # Raised when an input ring or polygon structure is invalid.
    class InvalidRingError < Error; end

    # One connected part of the polygon: +shell+ is an array of
    # Geokit::LatLng, +holes+ is an array of such arrays. Rings are open
    # (no repeated closing point) and oriented with the polygon interior on
    # the left of the traversal direction.
    Part = Struct.new(:shell, :holes)

    # Smallest accepted edge length, in radians.
    MIN_EDGE_ANGLE = 1e-9
    # Edges whose endpoints are this close to antipodal are rejected: a
    # half-circumference edge does not define a unique great circle.
    ANTIPODAL_DOT = -1.0 + 1e-9
    # Rings enclosing less than this many steradians are rejected.
    MIN_AREA = 1e-12
    # Rings whose interior is (nearly) exactly a hemisphere are ambiguous:
    # both sides have the same area, so the interior cannot be decided.
    HEMISPHERE_AREA = 2.0 * Math::PI

    attr_reader :parts

    # +parts+ is an array of part definitions. Each part is either a Hash
    # +{ shell: [...], holes: [[...], ...] }+ or a plain array of points
    # (treated as a shell without holes). A single flat array of points is
    # also accepted and treated as one shell. Points are Geokit::LatLng or
    # any object responding to +lat+/+lng+, or a +[lat, lng]+ pair.
    def initialize(parts = [])
      parts = [{ shell: parts }] if flat_point_list?(parts)
      @part_rings = []
      @rings = []
      @parts = []
      Array(parts).each do |part|
        part = { shell: part } unless part.is_a?(Hash)
        shell = normalize_ring(part.fetch(:shell), :shell)
        holes = part.fetch(:holes, []).map { |h| normalize_ring(h, :hole) }
        @part_rings << [shell, holes]
        @rings << [shell, 1]
        holes.each { |h| @rings << [h, -1] }
        @parts << Part.new(ring_to_lat_lngs(shell), holes.map { |h| ring_to_lat_lngs(h) })
      end
      validate_structure!
    end

    def empty?
      @rings.empty?
    end

    # Area covered by the polygon. With the default +radius+ of 1 the result
    # is in steradians; pass a sphere radius (e.g. in kilometres) to get the
    # corresponding squared unit.
    def area(radius = 1.0)
      total = @rings.sum do |vecs, sign|
        left = Spherical::Ring.left_area(vecs)
        sign.positive? ? left : left - 4.0 * Math::PI
      end
      total * radius * radius
    end

    # True when the point is inside the polygon or on its boundary.
    def contains?(point)
      v = to_vector(point)
      return true if boundary_vector?(v)
      membership_score(v) >= 1
    end

    # True when the point lies on an edge (within a small tolerance).
    def on_boundary?(point)
      boundary_vector?(to_vector(point))
    end

    def union(other)
      boolean(other, :union)
    end

    def intersection(other)
      boolean(other, :intersection)
    end

    def difference(other)
      boolean(other, :difference)
    end

    alias_method :|, :union
    alias_method :&, :intersection
    alias_method :-, :difference

    protected

    # Rings as [vectors, sign] pairs (+1 for shells, -1 for holes), oriented
    # with the polygon interior on the left. Used by the boolean engine.
    attr_reader :rings

    private

    # +1 inside a shell, -1 inside a hole void, summed over all rings: the
    # total is >= 1 exactly when the point is inside the polygon.
    def membership_score(v)
      @rings.sum do |vecs, sign|
        inside = Spherical::Ring.inside_left?(vecs, v) ? 1 : 0
        sign.positive? ? inside : inside - 1
      end
    end

    def boolean(other, op)
      engine = Spherical::BooleanEngine.new(rings, other.send(:rings))
      self.class.send(:from_engine_rings, engine.run(op))
    end

    def flat_point_list?(parts)
      parts.is_a?(Array) && !parts.empty? && parts.all? { |p| point_like?(p) }
    end

    def point_like?(p)
      (p.respond_to?(:lat) && p.respond_to?(:lng)) ||
        (p.is_a?(Array) && p.size == 2)
    end

    def to_vector(point)
      if point.respond_to?(:lat) && point.respond_to?(:lng)
        Spherical::VectorMath.lat_lng_to_vector(point.lat, point.lng)
      else
        Spherical::VectorMath.lat_lng_to_vector(point[0], point[1])
      end
    end

    def boundary_vector?(v)
      @rings.any? { |vecs, | Spherical::Ring.on_boundary?(vecs, v) }
    end

    # Converts to unit vectors, removes duplicate/closing vertices, checks
    # the ring is valid and orients it according to its role: shells get a
    # positive signed area (interior on the left), holes a negative one
    # (polygon interior on the left, void on the right).
    def normalize_ring(points, role)
      vecs = Array(points).map { |p| to_vector(p) }
      vecs = dedupe_vertices(vecs)
      if vecs.size < 3
        raise InvalidRingError,
              "ring has #{vecs.size} distinct vertices; at least 3 are required"
      end
      validate_edges!(vecs)
      validate_simple!(vecs)
      area = Spherical::Ring.signed_area(vecs)
      if area.abs < MIN_AREA
        raise InvalidRingError, 'ring encloses (nearly) zero area'
      end
      if (area.abs - HEMISPHERE_AREA).abs < 1e-9
        raise InvalidRingError,
              'ring encloses (nearly) exactly a hemisphere; the interior is ambiguous'
      end
      vecs.reverse! if (role == :shell && area < 0) || (role == :hole && area > 0)
      vecs
    end

    def dedupe_vertices(vecs)
      vecs = vecs.dup
      vecs.pop while vecs.size > 1 && same_vector?(vecs.last, vecs.first)
      out = [vecs.first]
      vecs[1..-1].each { |v| out << v unless same_vector?(v, out.last) }
      out
    end

    def same_vector?(a, b)
      (a[0] - b[0]).abs < 1e-15 && (a[1] - b[1]).abs < 1e-15 && (a[2] - b[2]).abs < 1e-15
    end

    def validate_edges!(vecs)
      vecs.each_index do |i|
        a = vecs[i]
        b = vecs[(i + 1) % vecs.size]
        d = Spherical::VectorMath.dot(a, b)
        angle = Math.atan2(Spherical::VectorMath.norm(Spherical::VectorMath.cross(a, b)), d)
        if angle < MIN_EDGE_ANGLE
          raise InvalidRingError,
                "degenerate edge #{i}: consecutive vertices are (nearly) identical"
        end
        if d < ANTIPODAL_DOT
          raise InvalidRingError,
                "edge #{i} is (nearly) a half circumference: antipodal endpoints "                 'do not define a unique great-circle arc'
        end
      end
    end

    def validate_simple!(vecs)
      n = vecs.size
      n.times do |i|
        (i + 1...n).each do |j|
          next if j == i + 1 || (i.zero? && j == n - 1)
          if Spherical::Ring.edges_intersect?(vecs[i], vecs[(i + 1) % n], vecs[j], vecs[(j + 1) % n])
            raise InvalidRingError,
                  "ring is self-intersecting or self-touching (edges #{i} and #{j})"
          end
        end
      end
    end

    def validate_structure!
      validate_rings_do_not_touch!
      validate_nesting!
    end

    def validate_rings_do_not_touch!
      @rings.map(&:first).combination(2) do |r1, r2|
        r1.each_index do |i|
          r2.each_index do |j|
            if Spherical::Ring.edges_intersect?(r1[i], r1[(i + 1) % r1.size], r2[j], r2[(j + 1) % r2.size])
              raise InvalidRingError, 'rings of the same polygon must not cross or touch'
            end
          end
        end
      end
    end

    def validate_nesting!
      @part_rings.each_with_index do |(shell, holes), pi|
        holes.each do |hole|
          unless Spherical::Ring.inside_left?(shell, hole[0])
            raise InvalidRingError, 'hole lies outside its shell'
          end
        end
        holes.combination(2) do |h1, h2|
          # A stored hole ring has the void on its right, so a vertex of
          # another hole must be inside its left region (outside the void).
          if !Spherical::Ring.inside_left?(h2, h1[0]) ||
             !Spherical::Ring.inside_left?(h1, h2[0])
            raise InvalidRingError, 'holes of the same polygon overlap'
          end
        end
        @part_rings.each_with_index do |(other_shell, _), pj|
          next if pi == pj
          if Spherical::Ring.inside_left?(other_shell, shell[0])
            raise InvalidRingError, 'shells of the same polygon overlap'
          end
        end
      end
    end

    def ring_to_lat_lngs(vecs)
      vecs.map do |v|
        lat, lng = Spherical::VectorMath.vector_to_lat_lng(v)
        Geokit::LatLng.new(lat, lng)
      end
    end

    class << self
      private

      # Builds a polygon from boolean-engine output: rings of unit vectors,
      # each with the polygon interior on its left. Rings with a positive
      # signed area are shells, the others are holes of the smallest shell
      # containing them.
      # Builds a polygon from boolean-engine output: rings of unit vectors,
      # each with the polygon interior on its left. A ring is a hole when
      # its void (the region on its right) is enclosed by another ring;
      # otherwise it is a shell. This nesting-based classification also
      # handles results that wrap the whole sphere (annulus-like unions),
      # where a shell can legitimately enclose more than a hemisphere.
      def from_engine_rings(rings)
        infos = rings.map do |ring|
          left = Spherical::Ring.left_area(ring)
          { ring: ring, left: left, void: 4.0 * Math::PI - left, container: nil }
        end
        infos.each do |info|
          sample = exterior_sample(info[:ring])
          infos.each do |other|
            next if other.equal?(info)
            next unless other[:left] > info[:void]
            next unless Spherical::Ring.inside_left?(other[:ring], sample)
            info[:container] = other if info[:container].nil? ||
                                        other[:left] < info[:container][:left]
          end
        end
        # Mutual containment happens for annulus-like results (two boundary
        # rings around the whole sphere); the ring with the larger left
        # region becomes the shell.
        infos.each do |info|
          container = info[:container]
          next unless container && container[:container] && container[:container].equal?(info)
          if info[:left] >= container[:left]
            info[:container] = nil
          else
            container[:container] = nil
          end
        end
        shells = infos.reject { |info| info[:container] }.map { |info| info[:ring] }
        holes = infos.select { |info| info[:container] }.map { |info| info[:ring] }
        grouped = shells.map { |s| [s, []] }
        holes.each do |hole|
          sample = exterior_sample(hole)
          best = nil
          grouped.each_with_index do |(shell, _), i|
            next unless Spherical::Ring.inside_left?(shell, sample)
            best = i if best.nil? ||
                        Spherical::Ring.left_area(shell) < Spherical::Ring.left_area(grouped[best][0])
          end
          raise Error, 'internal error: hole without a containing shell' if best.nil?
          grouped[best][1] << hole
        end
        poly = allocate
        poly.instance_variable_set(:@part_rings, grouped)
        poly.instance_variable_set(:@rings,
                                   shells.map { |s| [s, 1] } + holes.map { |h| [h, -1] })
        poly.instance_variable_set(:@parts, grouped.map do |shell, hs|
          Part.new(vec_ring_to_lat_lngs(shell), hs.map { |h| vec_ring_to_lat_lngs(h) })
        end)
        poly
      end

      # A point just outside the polygon next to the ring's first edge (the
      # polygon interior is on the left of every output ring, so the
      # exterior is on the right).
      def exterior_sample(ring)
        a, b = ring[0], ring[1]
        vm = Spherical::VectorMath
        n = vm.normalize(vm.cross(a, b))
        m = vm.normalize(vm.add(a, b))
        vm.normalize(vm.sub(m, vm.scale(n, 1e-7)))
      end

      def vec_ring_to_lat_lngs(vecs)
        vecs.map do |v|
          lat, lng = Spherical::VectorMath.vector_to_lat_lng(v)
          Geokit::LatLng.new(lat, lng)
        end
      end
    end
  end
end
