module Geokit
  # Validates that +value+ is a Numeric and returns it as a Float.
  # Boolean values are rejected on purpose.
  def self.assert_numeric(value, name = 'value')
    unless value.is_a?(Numeric)
      raise ArgumentError, "#{name} must be a Numeric, got #{value.class}"
    end
    Float(value)
  end

  # An Ellipsoid of revolution described by its semi-major axis +a+ (the
  # equatorial radius, measured in metres) and its inverse flattening +1/f+.
  # A sphere is represented by a nil inverse flattening.
  class Ellipsoid
    attr_reader :semi_major_axis, :inverse_flattening

    def initialize(semi_major_axis:, inverse_flattening: nil)
      @semi_major_axis = Geokit.assert_numeric(semi_major_axis,
                                              'semi_major_axis')
      unless @semi_major_axis.positive?
        raise ArgumentError, 'semi_major_axis must be greater than zero'
      end

      if inverse_flattening.nil?
        @inverse_flattening = nil
        @flattening = 0.0
      else
        @inverse_flattening = Geokit.assert_numeric(inverse_flattening,
                                                   'inverse_flattening')
        unless @inverse_flattening.positive?
          raise ArgumentError,
                'inverse_flattening must be greater than zero'
        end
        @flattening = 1.0 / @inverse_flattening
      end

      @semi_minor_axis = @semi_major_axis * (1.0 - @flattening)
      @first_eccentricity_squared = @flattening * (2.0 - @flattening)
      @second_eccentricity_squared =
        if @first_eccentricity_squared >= 1.0
          nil
        else
          @first_eccentricity_squared /
            (1.0 - @first_eccentricity_squared)
        end
    end

    # Flattening +f+; 0 for a sphere.
    def flattening
      @flattening
    end

    # Semi-minor axis +b+ (polar radius), in metres.
    def semi_minor_axis
      @semi_minor_axis
    end

    # First eccentricity squared, e^2 = f(2 - f).
    def first_eccentricity_squared
      @first_eccentricity_squared
    end

    # Second eccentricity squared, e'^2 = e^2 / (1 - e^2); nil for a sphere.
    def second_eccentricity_squared
      @second_eccentricity_squared
    end

    def sphere?
      @flattening.zero?
    end

    def ==(other)
      other.is_a?(Ellipsoid) &&
        semi_major_axis == other.semi_major_axis &&
        inverse_flattening == other.inverse_flattening
    end
    alias eql? ==

    def hash
      [semi_major_axis, inverse_flattening].hash
    end

    WGS84 = new(semi_major_axis: 6_378_137.0,
                inverse_flattening: 298.257_223_563)
    GRS80 = new(semi_major_axis: 6_378_137.0,
                inverse_flattening: 298.257_222_101)
    AIRY_1830 = new(semi_major_axis: 6_377_563.396,
                    inverse_flattening: 299.324_964_6)
    INTERNATIONAL_1924 = new(semi_major_axis: 6_378_388.0,
                             inverse_flattening: 297.0)
  end

end
