module Geokit
  module Spherical
    # Internal engine for boolean operations (union / intersection /
    # difference) between two spherical polygons.
    #
    # Strategy (arrangement based, no planar projection anywhere):
    #
    # 1. Collect every edge (great-circle minor arc) of both polygons.
    # 2. Compute all intersections between edges of different polygons,
    #    including shared vertices and collinear overlapping arcs, and split
    #    the edges there into atomic segments.
    # 3. Merge duplicate segments (shared arcs) so each geometric segment
    #    appears exactly once.
    # 4. Classify each segment by sampling a point slightly to its left and
    #    right and evaluating the boolean operation on the point-in-polygon
    #    results of both sides. A segment is kept (with a direction) when
    #    the operation result differs across the segment; the kept direction
    #    always has the result interior on its left. Tangencies and shared
    #    arcs need no special cases: they are decided by the two samples.
    # 5. Stitch kept directed segments into rings. At each vertex the walk
    #    continues with the first kept edge met when rotating clockwise
    #    from the reverse of the incoming edge, which pairs edges correctly
    #    even where the result touches itself at a vertex.
    #
    # The output rings all have the polygon interior on their left; rings
    # with positive signed area are shells, rings with negative signed area
    # are holes.
    class BooleanEngine
      include VectorMath

      COLLINEAR_TOL = 1e-9
      SNAP_DOT = 1.0 - 1e-12
      MIN_ANGLE = 1e-9
      SAMPLE_OFFSET = 1e-7

      def initialize(rings_a, rings_b)
        @polys = [rings_a, rings_b]
      end

      # op is :union, :intersection or :difference (a minus b).
      # Returns an array of rings; each ring is an array of unit vectors.
      def run(op)
        edges = collect_edges
        splits = compute_splits(edges)
        segments = build_segments(edges, splits)
        kept = select_segments(segments, op)
        assemble_rings(kept)
      end

      private

      def collect_edges
        edges = []
        @polys.each_with_index do |rings, poly|
          rings.each do |ring, _sign|
            ring.each_index do |i|
              edges << [ring[i], ring[(i + 1) % ring.size], poly]
            end
          end
        end
        edges
      end

      def compute_splits(edges)
        splits = Hash.new { |h, k| h[k] = [] }
        normals = edges.map { |a, b, _| cross(a, b) }
        edges.each_with_index do |(a1, a2, pa), i|
          (i + 1...edges.size).each do |j|
            b1, b2, pb = edges[j]
            next if pa == pb
            n1 = normals[i]
            n2 = normals[j]
            d = cross(n1, n2)
            dn = norm(d)
            if dn <= COLLINEAR_TOL * norm(n1) * norm(n2)
              # Same great circle: overlapping (or touching) arcs. Split each
              # arc at the endpoints of the other that lie on it.
              record_endpoint_split(splits, edges, i, b1)
              record_endpoint_split(splits, edges, i, b2)
              record_endpoint_split(splits, edges, j, a1)
              record_endpoint_split(splits, edges, j, a2)
            else
              p = scale(d, 1.0 / dn)
              [p, scale(p, -1.0)].each do |pt|
                next unless Ring.point_on_arc?(pt, a1, a2) && Ring.point_on_arc?(pt, b1, b2)
                pt = snap_point(pt, [a1, a2, b1, b2])
                add_split(splits, edges, i, pt)
                add_split(splits, edges, j, pt)
                break
              end
            end
          end
        end
        splits
      end

      def record_endpoint_split(splits, edges, idx, endpoint)
        a, b, = edges[idx]
        add_split(splits, edges, idx, endpoint) if Ring.point_on_arc?(endpoint, a, b)
      end

      # Intersection points that (nearly) coincide with an existing vertex
      # are replaced by that vertex so both polygons share identical
      # coordinates for the point.
      def snap_point(pt, vertices)
        vertices.each do |v|
          return v if dot(pt, v) >= SNAP_DOT
        end
        pt
      end

      def add_split(splits, edges, idx, pt)
        a, b, = edges[idx]
        theta = Math.atan2(norm(cross(a, b)), dot(a, b))
        t = Math.atan2(norm(cross(a, pt)), dot(a, pt))
        return if t < MIN_ANGLE || theta - t < MIN_ANGLE
        list = splits[idx]
        return if list.any? { |t0, _| (t0 - t).abs < MIN_ANGLE }
        list << [t, pt]
      end

      def build_segments(edges, splits)
        segments = []
        edges.each_with_index do |(a, b, _poly), i|
          chain = [a]
          splits[i].sort_by(&:first).each { |_t, pt| chain << pt }
          chain << b
          chain.each_cons(2) do |p, q|
            next if Math.atan2(norm(cross(p, q)), dot(p, q)) < MIN_ANGLE
            segments << [p, q]
          end
        end
        dedupe_segments(segments)
      end

      # Shared arcs produce the same atomic segment once per polygon; keep
      # only one copy (direction is irrelevant at this stage).
      def dedupe_segments(segments)
        seen = {}
        segments.each_with_object([]) do |seg, uniq|
          k1 = vertex_key(seg[0])
          k2 = vertex_key(seg[1])
          key = (k1 <=> k2) <= 0 ? [k1, k2] : [k2, k1]
          next if seen[key]
          seen[key] = true
          uniq << seg
        end
      end

      def select_segments(segments, op)
        kept = []
        segments.each do |p, q|
          n = cross(p, q)
          nn = norm(n)
          next if nn < 1e-15
          n = scale(n, 1.0 / nn)
          m = normalize(add(p, q))
          left = normalize(add(m, scale(n, SAMPLE_OFFSET)))
          right = normalize(add(m, scale(n, -SAMPLE_OFFSET)))
          lres = apply_op(op, inside_polygon?(left, 0), inside_polygon?(left, 1))
          rres = apply_op(op, inside_polygon?(right, 0), inside_polygon?(right, 1))
          kept << [p, q] if lres && !rres
          kept << [q, p] if rres && !lres
        end
        kept
      end

      def inside_polygon?(pt, poly)
        @polys[poly].sum do |ring, sign|
          inside = Ring.inside_left?(ring, pt) ? 1 : 0
          sign.positive? ? inside : inside - 1
        end >= 1
      end

      def apply_op(op, a, b)
        case op
        when :union then a || b
        when :intersection then a && b
        when :difference then a && !b
        else raise ArgumentError, "unknown boolean operation: #{op.inspect}"
        end
      end

      def assemble_rings(kept)
        outgoing = Hash.new { |h, k| h[k] = [] }
        kept.each_with_index { |(p, _q), i| outgoing[vertex_key(p)] << i }
        angles = kept.map { |p, q| tangent_angle(p, q) }
        visited = Array.new(kept.size, false)
        rings = []
        kept.each_index do |start|
          next if visited[start]
          ring = []
          current = start
          loop do
            if ring.size > kept.size
              raise SphericalPolygon::Error, 'ring assembly failed: boundary is not closed'
            end
            visited[current] = true
            p, q = kept[current]
            ring << p
            candidates = outgoing[vertex_key(q)]
            if candidates.empty?
              raise SphericalPolygon::Error, 'ring assembly failed: unbalanced boundary'
            end
            nxt = if candidates.size == 1
                    candidates[0]
                  else
                    rev = tangent_angle(q, p)
                    candidates.min_by { |ci| (rev - angles[ci]) % (2.0 * Math::PI) }
                  end
            break if nxt == start
            current = nxt
          end
          rings << ring
        end
        rings
      end

      # Heading of the geodesic from +at+ towards +to+, measured in the
      # tangent plane at +at+ against a deterministic local basis.
      def tangent_angle(at, to)
        dir = sub(to, scale(at, dot(at, to)))
        axis = at[2].abs < 0.9 ? [0.0, 0.0, 1.0] : [1.0, 0.0, 0.0]
        e1 = normalize(cross(at, axis))
        e2 = cross(at, e1)
        Math.atan2(dot(dir, e2), dot(dir, e1))
      end

      def vertex_key(v)
        [v[0].round(9), v[1].round(9), v[2].round(9)]
      end
    end
  end
end
