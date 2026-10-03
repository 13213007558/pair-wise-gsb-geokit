module Geokit
  module Spherical
    # Low level spherical geometry: great-circle short arcs, signed ring
    # areas and winding numbers. Everything works on unit 3D vectors, so no
    # planar projection is ever involved.
    module Core
      # Angular tolerance (radians) used to decide whether a point lies on
      # an arc. ~1e-9 rad is about 6 mm on the Earth surface.
      ANG_TOL = 1e-9

      module_function

      # Is unit vector +p+ on the short great-circle arc from +a+ to +b+?
      def on_arc?(a, b, p, tol = ANG_TOL)
        n = Vector.cross(a, b)
        ln = Vector.norm(n)
        return false if ln < 1e-15
        nh = Vector.scale(n, 1.0 / ln)
        return false if Vector.dot(p, nh).abs > tol
        Vector.dot(Vector.cross(a, p), nh) >= -tol &&
          Vector.dot(Vector.cross(p, b), nh) >= -tol
      end

      # Intersections between the short arcs a->b and c->d.
      #
      # Returns [type, points] where type is:
      # * :cross      - the supporting great circles are distinct; +points+
      #                 holds 0, 1 or 2 unit vectors where the arcs meet
      #                 (proper crossings as well as endpoint touches).
      # * :coincident - both arcs lie on the same great circle; +points+
      #                 holds the endpoints of either arc that lie on the
      #                 other arc (empty when the arcs are disjoint, one
      #                 point when they touch, two or more when they
      #                 overlap over a non-degenerate interval).
      def arc_intersections(a, b, c, d, tol = ANG_TOL)
        n1 = Vector.cross(a, b)
        n2 = Vector.cross(c, d)
        x = Vector.cross(n1, n2)
        s = Vector.norm(x)
        if s <= tol * Vector.norm(n1) * Vector.norm(n2)
          points = []
          points << a if on_arc?(c, d, a, tol)
          points << b if on_arc?(c, d, b, tol)
          points << c if on_arc?(a, b, c, tol)
          points << d if on_arc?(a, b, d, tol)
          [:coincident, points]
        else
          xn = Vector.scale(x, 1.0 / s)
          points = []
          [xn, Vector.negate(xn)].each do |p|
            points << p if on_arc?(a, b, p, tol) && on_arc?(c, d, p, tol)
          end
          [:cross, points]
        end
      end

      # Signed solid angle of the triangle (a, b, c) of unit vectors, in
      # (-2PI, 2PI). Positive when a, b, c are counterclockwise as seen
      # from outside the sphere.
      def tri_solid_angle(a, b, c)
        num = Vector.dot(a, Vector.cross(b, c))
        den = 1.0 + Vector.dot(a, b) + Vector.dot(b, c) + Vector.dot(c, a)
        2.0 * Math.atan2(num, den)
      end

      # Signed spherical area (steradians) enclosed by a ring of unit
      # vectors, with the interior conventionally on the left of travel.
      # Computed by fan triangulation, so it is exact up to float rounding
      # and independent of any map projection. Positive for rings whose
      # (smaller than a hemisphere) interior is on the left.
      def signed_ring_area(ring)
        return 0.0 if ring.size < 3
        a0 = ring[0]
        sum = 0.0
        (1...(ring.size - 1)).each do |i|
          sum += tri_solid_angle(a0, ring[i], ring[i + 1])
        end
        sum
      end

      # Winding number of a ring of unit vectors around unit vector +p+.
      # Roughly +1 or -1 when +p+ is inside the ring, ~0 when outside.
      # +p+ must not lie on the ring itself.
      def ring_winding_number(ring, p)
        sum = 0.0
        ring.each_index do |i|
          a = ring[i]
          b = ring[(i + 1) % ring.size]
          u = Vector.sub(a, Vector.scale(p, Vector.dot(a, p)))
          v = Vector.sub(b, Vector.scale(p, Vector.dot(b, p)))
          sum += Math.atan2(Vector.dot(Vector.cross(u, v), p), Vector.dot(u, v))
        end
        sum / (2.0 * Math::PI)
      end

      # Combined winding number of a set of oriented rings (outer rings
      # counterclockwise, holes clockwise) around +p+. Points inside the
      # polygon yield ~1, points inside a hole or outside yield ~0.
      def polygon_winding_number(rings, p)
        rings.sum { |ring| ring_winding_number(ring, p) }
      end

      # Is +p+ on any edge of the given rings?
      def on_boundary?(rings, p, tol = ANG_TOL)
        rings.any? do |ring|
          ring.each_index.any? do |i|
            on_arc?(ring[i], ring[(i + 1) % ring.size], p, tol)
          end
        end
      end
    end
  end
end
