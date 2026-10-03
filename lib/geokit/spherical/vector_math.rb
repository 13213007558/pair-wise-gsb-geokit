module Geokit
  module Spherical
    # Internal 3D vector helpers. All spherical geometry is computed with
    # local coordinates only: points are unit vectors in R3, edges are
    # great-circle minor arcs (strictly shorter than a half circumference).
    # No planar projection is ever used, so the antimeridian and the poles
    # are not special in any way.
    module VectorMath
      module_function

      def dot(a, b)
        a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
      end

      def cross(a, b)
        [
          a[1] * b[2] - a[2] * b[1],
          a[2] * b[0] - a[0] * b[2],
          a[0] * b[1] - a[1] * b[0]
        ]
      end

      def norm(a)
        Math.sqrt(dot(a, a))
      end

      def normalize(a)
        n = norm(a)
        [a[0] / n, a[1] / n, a[2] / n]
      end

      def add(a, b)
        [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
      end

      def sub(a, b)
        [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
      end

      def scale(a, s)
        [a[0] * s, a[1] * s, a[2] * s]
      end

      # Shortest angular distance between two unit vectors, in radians.
      def angle_between(a, b)
        Math.atan2(norm(cross(a, b)), dot(a, b))
      end

      # lat/lng in degrees -> unit vector. Longitude is not normalized first;
      # any equivalent representation (e.g. 190 vs -170) maps to the same
      # vector up to floating point rounding.
      def lat_lng_to_vector(lat, lng)
        phi = lat.to_f * Math::PI / 180.0
        lam = lng.to_f * Math::PI / 180.0
        c = Math.cos(phi)
        [c * Math.cos(lam), c * Math.sin(lam), Math.sin(phi)]
      end

      # unit vector -> [lat, lng] in degrees, lng in (-180, 180]
      def vector_to_lat_lng(v)
        z = v[2]
        z = 1.0 if z > 1.0
        z = -1.0 if z < -1.0
        [Math.asin(z) * 180.0 / Math::PI, Math.atan2(v[1], v[0]) * 180.0 / Math::PI]
      end
    end
  end
end
