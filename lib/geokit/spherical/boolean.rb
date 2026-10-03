module Geokit
  module Spherical
    # Boolean operations (union / intersection / difference) on spherical
    # polygons.
    #
    # Algorithm: an arrangement of great-circle short arcs on the unit
    # sphere.
    #
    # 1. All edge/edge intersections between the two operands are computed
    #    exactly (great-circle crossings, endpoint touches and coincident
    #    arc overlaps) and every edge is split at those points.
    # 2. Each resulting segment is classified against the other operand by
    #    testing its midpoint with a winding number. Coincident segments
    #    (shared arcs) are detected structurally by identical endpoints,
    #    never by guessing with a tolerance.
    # 3. Segments are kept or dropped according to the operation and
    #    stitched into rings: at every vertex the outgoing segment making
    #    the smallest counterclockwise angle with the incoming direction
    #    is followed, which keeps the result interior on the left and
    #    resolves shared vertices and tangencies deterministically.
    # 4. Rings are grouped into polygons (outer rings + holes) by nesting
    #    depth.
    #
    # All kept segments carry the direction in which the result interior
    # lies on their left, so output rings follow the same convention as
    # the input: outer rings counterclockwise, holes clockwise.
    class Boolean
      # Chord/angle tolerance (radians) used to merge intersection points
      # that are mathematically identical but differ by float noise.
      VERTEX_TOL = 1e-9

      Segment = Struct.new(:from, :to, :source)

      def self.union(a, b)
        new(a, b).run(:union)
      end

      def self.intersection(a, b)
        new(a, b).run(:intersection)
      end

      def self.difference(a, b)
        new(a, b).run(:difference)
      end

      def initialize(a, b)
        @rings_a = rings_of(a)
        @rings_b = rings_of(b)
        @pool = []
      end

      def run(op)
        edges_a = @rings_a.flat_map { |ring| ring_edges(ring) }
        edges_b = @rings_b.flat_map { |ring| ring_edges(ring) }
        segments_a = split_all(edges_a, edges_b, :a)
        segments_b = split_all(edges_b, edges_a, :b)
        kept = select_segments(segments_a, segments_b, op)
        build_result(stitch(kept))
      end

      private

      def rings_of(operand)
        case operand
        when MultiPolygon
          operand.polygons.flat_map(&:vector_rings)
        when Polygon
          operand.vector_rings
        else
          raise ArgumentError, "expected a Geokit::Spherical::Polygon or MultiPolygon, got #{operand.class}"
        end
      end

      # Canonicalize a vertex so that mathematically identical points are
      # the same object and can be compared with equal?.
      def intern(v)
        existing = @pool.find { |p| Vector.angle_between(p, v) < VERTEX_TOL }
        return existing if existing
        @pool << v
        v
      end

      def ring_edges(ring)
        ring.each_index.map do |i|
          [intern(ring[i]), intern(ring[(i + 1) % ring.size])]
        end
      end

      def split_all(edges, other_edges, source)
        splits = Array.new(edges.size) { [] }
        edges.each_with_index do |edge, i|
          other_edges.each do |other|
            _type, points = Core.arc_intersections(edge[0], edge[1], other[0], other[1])
            points.each do |p|
              ip = intern(p)
              splits[i] << ip unless ip.equal?(edge[0]) || ip.equal?(edge[1])
            end
          end
        end
        edges.each_with_index.flat_map do |edge, i|
          split_edge(edge[0], edge[1], splits[i].uniq, source)
        end
      end

      def split_edge(a, b, points, source)
        nh = Vector.normalize(Vector.cross(a, b))
        sorted = points.sort_by do |p|
          Math.atan2(Vector.dot(Vector.cross(a, p), nh), Vector.dot(a, p))
        end
        chain = [a] + sorted + [b]
        chain.each_cons(2).filter_map do |u, v|
          Segment.new(u, v, source) unless u.equal?(v)
        end
      end

      def select_segments(segments_a, segments_b, op)
        unmatched_b = segments_b.dup
        kept = []
        segments_a.each do |sa|
          shared = unmatched_b.find { |sb| same_pair?(sa, sb) }
          if shared
            unmatched_b.delete(shared)
            same_direction = sa.from.equal?(shared.from)
            case op
            when :union, :intersection
              kept << sa if same_direction
            when :difference
              kept << sa unless same_direction
            end
          else
            inside = polygon_contains?(@rings_b, midpoint(sa))
            case op
            when :union, :difference
              kept << sa unless inside
            when :intersection
              kept << sa if inside
            end
          end
        end
        unmatched_b.each do |sb|
          inside = polygon_contains?(@rings_a, midpoint(sb))
          case op
          when :union
            kept << sb unless inside
          when :intersection
            kept << sb if inside
          when :difference
            kept << Segment.new(sb.to, sb.from, sb.source) if inside
          end
        end
        kept
      end

      def same_pair?(s1, s2)
        (s1.from.equal?(s2.from) && s1.to.equal?(s2.to)) ||
          (s1.from.equal?(s2.to) && s1.to.equal?(s2.from))
      end

      def midpoint(segment)
        Vector.normalize(Vector.add(segment.from, segment.to))
      end

      def polygon_contains?(rings, p)
        if Core.on_boundary?(rings, p)
          raise Error, 'internal error: classification point lies on the other operand boundary'
        end
        Core.polygon_winding_number(rings, p) > 0.5
      end

      def stitch(segments)
        outgoing = Hash.new { |hash, key| hash[key] = [] }
        segments.each { |s| outgoing[s.from.object_id] << s }
        unused = {}
        segments.each { |s| unused[s.object_id] = true }
        rings = []
        segments.each do |first|
          next unless unused.delete(first.object_id)
          ring = [first.from]
          current = first
          loop do
            ring << current.to
            break if current.to.equal?(first.from)
            candidates = outgoing[current.to.object_id].select { |s| unused[s.object_id] }
            raise Error, 'unbalanced boundary graph' if candidates.empty?
            nxt = candidates.min_by { |s| turn_angle(current, s) }
            unused.delete(nxt.object_id)
            current = nxt
            raise Error, 'boundary ring failed to close' if ring.size > segments.size + 1
          end
          rings << ring
        end
        rings
      end

      # Counterclockwise angle in [0, 2PI) from the arrival direction of
      # +seg_in+ to the departure direction of +seg_out+, measured in the
      # tangent plane at their shared vertex.
      def turn_angle(seg_in, seg_out)
        v = seg_in.to
        t_in = Vector.normalize(Vector.cross(Vector.cross(seg_in.from, seg_in.to), v))
        t_out = Vector.normalize(Vector.cross(Vector.cross(seg_out.from, seg_out.to), v))
        angle = Math.atan2(Vector.dot(Vector.cross(t_in, t_out), v), Vector.dot(t_in, t_out))
        angle += 2.0 * Math::PI if angle.negative?
        angle
      end

      def build_result(rings)
        return MultiPolygon.new([]) if rings.empty?
        ring_vecs = rings.map { |ring| ring[0...-1] }
        count = ring_vecs.size
        contains = Array.new(count) { Array.new(count, false) }
        (0...count).each do |i|
          (0...count).each do |j|
            next if i == j
            probe = edge_midpoint(ring_vecs[j])
            contains[i][j] = Core.ring_winding_number(ring_vecs[i], probe).abs > 0.5
          end
        end
        depth = (0...count).map { |j| (0...count).count { |i| contains[i][j] } }
        enclosed = ring_vecs.map do |rv|
          area = Core.signed_ring_area(rv)
          area.negative? ? area + 4.0 * Math::PI : area
        end
        polygons = (0...count).select { |i| depth[i].even? }.map do |outer_idx|
          hole_indexes = (0...count).select do |h|
            depth[h].odd? && contains[outer_idx][h] &&
              (0...count).none? do |other|
                other != outer_idx && depth[other].even? &&
                  contains[other][h] && enclosed[other] < enclosed[outer_idx]
              end
          end
          Polygon.new(
            to_lat_lng_ring(ring_vecs[outer_idx]),
            hole_indexes.map { |h| to_lat_lng_ring(ring_vecs[h]) },
            validate: false
          )
        end
        MultiPolygon.new(polygons)
      end

      def edge_midpoint(ring)
        Vector.normalize(Vector.add(ring[0], ring[1]))
      end

      def to_lat_lng_ring(ring)
        ring.map { |v| Geokit::LatLng.new(*Vector.to_lat_lng(v)) }
      end
    end
  end
end
