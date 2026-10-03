module Geokit
  module Spherical
    # Internal predicates and measures for rings of unit vectors.
    #
    # Ring orientation convention used across the library: the interior of
    # the polygon is always on the LEFT of the traversal direction. For a
    # shell this means the bounded region on the left is inside the polygon;
    # for a hole the ring is traversed the other way so that the polygon
    # interior is still on the left (the void is on the right). Every ring
    # interior must be smaller than a hemisphere.
    module Ring
      extend VectorMath

      class << self
        # Signed spherical area of the region on the left of the ring, in
        # steradians, normalized to (-2*pi, 2*pi): positive when the left
        # region is smaller than a hemisphere, negative when it is larger
        # (e.g. for a stored hole ring the result is minus the void area).
        def signed_area(vecs)
          area = left_area(vecs)
          area >= 2.0 * Math::PI ? area - 4.0 * Math::PI : area
        end

        # Area of the region on the left of the ring, in [0, 4*pi).
        def left_area(vecs)
          fan_solid_angle(vecs, vecs[0]) % (4.0 * Math::PI)
        end

        # Sum of the signed solid angles of the fan triangles (+apex+, a, b)
        # over the ring edges (van Oosterom & Strackee formula). The result
        # equals the left-region area modulo 4*pi.
        def fan_solid_angle(vecs, apex)
          sum = 0.0
          vecs.each_index do |i|
            a = vecs[i]
            b = vecs[(i + 1) % vecs.size]
            num = dot(apex, cross(a, b))
            den = 1.0 + dot(apex, a) + dot(a, b) + dot(b, apex)
            sum += 2.0 * Math.atan2(num, den)
          end
          sum
        end

        # True when +point+ lies inside the region on the left of the ring.
        #
        # Spherical ray casting: a reference point is placed just to the
        # RIGHT of a ring edge, which for a simple ring is guaranteed to be
        # outside the left region. The point is inside exactly when the
        # minor arc from +point+ to the reference crosses the ring an odd
        # number of times. When +point+ is very close to the reference, the
        # question is answered by a local side test on the nearest edge
        # instead. This uses only local vector arithmetic (no projection, no
        # antimeridian or pole special cases). Points exactly on the ring
        # are not handled here; see on_boundary?.
        def inside_left?(vecs, point)
          anchors = edge_anchors(vecs)
          anchors.each do |(normal, mid)|
            [0.0, 0.9, -0.9].each do |angle|
              dir = rotate_around(normal, mid, angle)
              q = normalize(sub(mid, scale(dir, 1e-6)))
              if chord(point, q) < 1e-5
                # Point is close to this anchor: decide locally by the side
                # of the nearest edge.
                nearest = nearest_edge(vecs, point)
                return dot(nearest, point) > 0
              end
              result = crossing_parity(vecs, point, q)
              return result unless result.nil?
            end
          end
          raise SphericalPolygon::Error, 'point-in-ring test failed: degenerate geometry'
        end

        # Parity of proper crossings between the minor arc (p, q) and the
        # ring (true/false for odd/even), or nil when the configuration is
        # degenerate (the arc is ill-defined, or an endpoint lies on the
        # ring, or a ring vertex lies on the arc, or a ring edge is nearly
        # collinear with the arc). Crossings are computed as actual
        # intersection points of the supporting great circles; note that the
        # planar "straddle" shortcut is NOT valid on the sphere (the two
        # circles meet twice, possibly on different arcs).
        def crossing_parity(vecs, p, q)
          n_arc = cross(p, q)
          return nil if norm(n_arc) < 1e-12
          n_arc = normalize(n_arc)
          crossings = 0
          vecs.each_index do |i|
            a = vecs[i]
            b = vecs[(i + 1) % vecs.size]
            n_e = cross(a, b)
            return nil if norm(n_e) < 1e-15
            n_e = normalize(n_e)
            d = cross(n_arc, n_e)
            dn = norm(d)
            return nil if dn < 1e-3
            d = scale(d, 1.0 / dn)
            [d, scale(d, -1.0)].each do |cand|
              next unless on_arc_strict?(cand, p, q) && point_on_arc?(cand, a, b)
              return nil if chord(cand, a) < 1e-9 || chord(cand, b) < 1e-9
              crossings += 1
            end
          end
          crossings.odd?
        end

        # True when +point+ lies on one of the ring edges (within tolerance).
        def on_boundary?(vecs, point, tol = 1e-9)
          vecs.each_index do |i|
            a = vecs[i]
            b = vecs[(i + 1) % vecs.size]
            n = cross(a, b)
            nn = norm(n)
            next if nn.zero?
            next unless dot(n, point).abs <= tol * nn
            return true if point_on_arc?(point, a, b)
          end
          false
        end

        # True when point +p+ lies on the minor arc from +a+ to +b+, given
        # that +p+ is (numerically) on the same great circle. Tested via
        # the arc midpoint: +p+ must be at least as close to the midpoint
        # direction as the endpoints are. (The naive "close to both
        # endpoints" test wrongly accepts major-arc points for arcs wider
        # than 120 degrees.)
        def point_on_arc?(p, a, b)
          m = normalize(add(a, b))
          dot(p, m) >= dot(a, m) - 1e-12
        end

        # True when the minor arcs (a1, a2) and (b1, b2) have any point in
        # common: a proper crossing, a shared/touching endpoint, or a
        # collinear overlap. Computed via actual intersection points, not
        # the (on the sphere insufficient) straddle test.
        def edges_intersect?(a1, a2, b1, b2)
          n1 = cross(a1, a2)
          n2 = cross(b1, b2)
          d = cross(n1, n2)
          dn = norm(d)
          if dn <= 1e-9 * norm(n1) * norm(n2)
            return point_on_arc?(a1, b1, b2) || point_on_arc?(a2, b1, b2) ||
                   point_on_arc?(b1, a1, a2) || point_on_arc?(b2, a1, a2)
          end
          p = scale(d, 1.0 / dn)
          q = scale(p, -1.0)
          (point_on_arc?(p, a1, a2) && point_on_arc?(p, b1, b2)) ||
            (point_on_arc?(q, a1, a2) && point_on_arc?(q, b1, b2))
        end

        # Chord (straight-line) distance between two unit vectors.
        def chord(a, b)
          norm(sub(a, b))
        end

        # Rotates vector +v+ around unit axis +axis+ by +angle+ (Rodrigues).
        def rotate_around(v, axis, angle)
          c = Math.cos(angle)
          s = Math.sin(angle)
          add(add(scale(v, c), scale(cross(axis, v), s)),
              scale(axis, dot(axis, v) * (1.0 - c)))
        end

        private

        # Unit normal of the ring edge that +point+ is closest to (by
        # perpendicular distance to the edge's great circle when the
        # projection falls on the arc, otherwise by endpoint distance).
        def nearest_edge(vecs, point)
          best = nil
          best_dist = Float::INFINITY
          vecs.each_index do |i|
            a = vecs[i]
            b = vecs[(i + 1) % vecs.size]
            n = cross(a, b)
            nn = norm(n)
            next if nn < 1e-15
            n = scale(n, 1.0 / nn)
            dist = if point_on_arc?(point, a, b)
                     dot(n, point).abs
                   else
                     [chord(point, a), chord(point, b)].min
                   end
            if dist < best_dist
              best_dist = dist
              best = n
            end
          end
          best
        end

        # [normal, midpoint] pairs for each edge: the unit normal of the
        # edge's plane (pointing left of the traversal) and the normalized
        # edge midpoint.
        def edge_anchors(vecs)
          vecs.each_index.map do |i|
            a = vecs[i]
            b = vecs[(i + 1) % vecs.size]
            n = cross(a, b)
            n = normalize(n)
            [n, normalize(add(a, b))]
          end
        end

        # Strict on-arc test used for crossing points: +cand+ must lie
        # between +p+ and +q+ on their common great circle, with a small
        # tolerance for floating point error in +cand+.
        def on_arc_strict?(cand, p, q)
          m = normalize(add(p, q))
          dot(cand, m) >= dot(p, m) - 1e-13
        end
      end
    end
  end
end
