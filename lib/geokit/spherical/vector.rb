module Geokit
  module Spherical
    # 3D vector helpers on the unit sphere. Vectors are plain [x, y, z]
    # arrays of floats. Working in 3D (instead of projecting lat/lng to a
    # plane) keeps every computation local and behaves identically across
    # the antimeridian and near the poles.
    module Vector
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
        raise ArgumentError, 'cannot normalize a zero vector' if n == 0.0
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

      def negate(a)
        [-a[0], -a[1], -a[2]]
      end

      # Angle between two vectors in radians, in [0, PI].
      def angle_between(a, b)
        Math.atan2(norm(cross(a, b)), dot(a, b))
      end

      # Unit vector for a lat/lng pair given in degrees.
      def from_lat_lng(lat_deg, lng_deg)
        lat = lat_deg.to_f * Math::PI / 180.0
        lng = lng_deg.to_f * Math::PI / 180.0
        cos_lat = Math.cos(lat)
        [cos_lat * Math.cos(lng), cos_lat * Math.sin(lng), Math.sin(lat)]
      end

      # [lat, lng] in degrees for a unit vector.
      def to_lat_lng(v)
        lat = Math.asin(v[2].clamp(-1.0, 1.0)) * 180.0 / Math::PI
        lng = Math.atan2(v[1], v[0]) * 180.0 / Math::PI
        [lat, lng]
      end
    end
  end
end
