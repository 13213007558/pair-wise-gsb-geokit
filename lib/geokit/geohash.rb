module Geokit
  # Encode and inspect standard base-32 Geohash prefixes.
  module Geohash
    BASE32 = '0123456789bcdefghjkmnpqrstuvwxyz'.freeze
    MAX_PRECISION = 12

    def self.encode(lat, lng, precision)
      validate_precision(precision)

      lat = lat.to_f
      lng = lng.to_f
      prefix = String.new(encoding: Encoding::UTF_8)
      lat_range = [-90.0, 90.0]
      lng_range = [-180.0, 180.0]
      value = 0

      (precision * 5).times do |global_bit|
        if global_bit.even?
          value = append_bit(value, lng, lng_range)
        else
          value = append_bit(value, lat, lat_range)
        end

        if global_bit % 5 == 4
          prefix << BASE32[value]
          value = 0
        end
      end

      prefix
    end

    def self.bounds(prefix = nil)
      prefix = prefix.to_s.downcase
      return default_bounds if prefix.empty?

      longitude_bits = 0
      latitude_bits = 0
      longitude_index = 0
      latitude_index = 0

      prefix.each_char.with_index do |character, character_index|
        value = BASE32.index(character)
        raise ArgumentError, "#{character.inspect} is not a Geohash character" if value.nil?

        4.downto(0) do |bit_offset|
          bit = (value >> bit_offset) & 1
          global_bit = (character_index * 5) + (4 - bit_offset)
          if global_bit.even?
            longitude_index = (longitude_index << 1) | bit
            longitude_bits += 1
          else
            latitude_index = (latitude_index << 1) | bit
            latitude_bits += 1
          end
        end
      end

      longitude_bounds(longitude_index, longitude_bits).merge(
        latitude_bounds(latitude_index, latitude_bits)
      )
    end

    def self.validate_precision(precision)
      precision = Integer(precision)
      return if precision >= 0 && precision <= MAX_PRECISION

      raise ArgumentError, "precision must be between 0 and #{MAX_PRECISION}"
    end

    def self.default_bounds
      { south: -90.0, west: -180.0, north: 90.0, east: 180.0 }
    end

    def self.longitude_bounds(index, bits)
      return { west: -180.0, east: 180.0 } if bits.zero?

      width = 360.0 / (1 << bits)
      west = -180.0 + index * width
      east = index == (1 << bits) - 1 ? 180.0 : west + width
      { west: west, east: east }
    end

    def self.latitude_bounds(index, bits)
      return { south: -90.0, north: 90.0 } if bits.zero?

      width = 180.0 / (1 << bits)
      south = -90.0 + index * width
      north = index == (1 << bits) - 1 ? 90.0 : south + width
      { south: south, north: north }
    end

    def self.append_bit(value, coordinate, range)
      middle = (range[0] + range[1]) / 2.0
      if coordinate >= middle
        range[0] = middle
        ((value << 1) | 1)
      else
        range[1] = middle
        (value << 1)
      end
    end

    private_class_method :append_bit, :default_bounds, :latitude_bounds,
                         :longitude_bounds
  end
end
