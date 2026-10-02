module Geokit
  # Vector helpers for spherical cap calculations.
  module SphericalCapVectors
    private

    def dot(left, right)
      left[0] * right[0] + left[1] * right[1] + left[2] * right[2]
    end

    def cross(left, right)
      [
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
      ]
    end

    def add(left, right)
      [left[0] + right[0], left[1] + right[1], left[2] + right[2]]
    end

    def subtract(left, right)
      [left[0] - right[0], left[1] - right[1], left[2] - right[2]]
    end

    def scale(vector, factor)
      vector.map { |value| value * factor }
    end

    def norm_squared(vector)
      dot(vector, vector)
    end

    def norm(vector)
      Math.sqrt(norm_squared(vector))
    end

    def normalize(vector)
      length = norm(vector)
      raise ArgumentError, 'cannot normalize a zero vector' if length < 1e-12
      scale(vector, 1.0 / length)
    end

    def clamp(value, minimum, maximum)
      [minimum, [value, maximum].min].max
    end

    def acos_clamped(value)
      Math.acos(clamp(value, -1.0, 1.0))
    end

    def deg2rad(value)
      LatLng.deg2rad(value)
    end

    def rad2deg(value)
      LatLng.rad2deg(value)
    end

    def to_vector(location)
      latitude = deg2rad(location.lat)
      longitude = deg2rad(location.lng)
      cos_latitude = Math.cos(latitude)
      [
        cos_latitude * Math.cos(longitude),
        cos_latitude * Math.sin(longitude),
        Math.sin(latitude),
      ]
    end

    def to_lat_lng(vector)
      vector = normalize(vector)
      latitude = rad2deg(Math.asin(clamp(vector[2], -1.0, 1.0)))
      longitude = rad2deg(Math.atan2(vector[1], vector[0]))
      latitude = 0.0 if latitude.zero?
      longitude = 0.0 if longitude.zero?
      LatLng.new(latitude, longitude)
    end

    def lex_smaller?(left, right)
      3.times do |index|
        difference = left[index] - right[index]
        next if difference.abs <= 1e-12
        return difference.negative?
      end
      false
    end

    def same_vector?(left, right)
      norm(subtract(left, right)) < 1e-12
    end
  end
end
