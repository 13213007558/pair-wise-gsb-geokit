module Geokit
  # A seven-parameter Helmert (Bursa-Wolf) similarity transformation
  # between two ECEF Cartesian coordinate systems:
  #
  #   X_target = t(epoch) + (1 + s(epoch)) * R(convention, epoch) * X_source
  #
  # +t+ is a translation in metres, +s+ the dimensionless scale change and
  # +R+ the small-angle rotation matrix about the X, Y and Z axes.
  #
  # Published parameter sets follow two opposite rotation conventions and
  # stating one is mandatory:
  # * +:position_vector+ (EPSG position-vector / Bursa-Wolf):
  #     R = I + [[0, -rz, ry], [rz, 0, -rx], [-ry, rx, 0]]
  # * +:coordinate_frame+ (coordinate-frame rotation): the rotation signs
  #   are reversed.
  #
  # Units are never guessed: translation, rotation and scale units must be
  # supplied, along with the rate units when parameter rates are used.
  # Epochs are decimal years (e.g. 2010.0); a time-dependent
  # transformation requires a +reference_epoch+ and an epoch at each call.
  class Helmert
    RAD = Math::PI / 180.0

    CONVENTIONS = %i[position_vector coordinate_frame].freeze

    # Multipliers converting the stated translation unit into metres.
    TRANSLATION_UNITS = {
      metres: 1.0,
      meters: 1.0,
      millimetres: 0.001,
      feet: 0.3048
    }.freeze

    # Multipliers converting the stated rotation unit into radians.
    ROTATION_UNITS = {
      radians: 1.0,
      degrees: RAD,
      arcseconds: RAD / 3600.0,
      milliarcseconds: RAD / 3_600_000.0
    }.freeze

    # Multipliers converting the stated scale unit into the dimensionless
    # scale change (ppm = parts per million).
    SCALE_UNITS = {
      ratio: 1.0,
      ppm: 1.0e-6,
      ppb: 1.0e-9
    }.freeze

    RATE_TRANSLATION_UNITS = {
      metres_per_year: 1.0,
      meters_per_year: 1.0,
      millimetres_per_year: 0.001,
      feet_per_year: 0.3048
    }.freeze

    RATE_ROTATION_UNITS = {
      radians_per_year: 1.0,
      degrees_per_year: RAD,
      arcseconds_per_year: RAD / 3600.0,
      milliarcseconds_per_year: RAD / 3_600_000.0
    }.freeze

    RATE_SCALE_UNITS = {
      ratio_per_year: 1.0,
      ppm_per_year: 1.0e-6,
      ppb_per_year: 1.0e-9
    }.freeze

    RATE_KEYS = %i[translation rotation scale].freeze

    attr_reader :name, :convention, :reference_epoch

    def initialize(convention:, translation:, rotation: [0.0, 0.0, 0.0],
                   scale: 0.0, translation_unit:, rotation_unit: :arcseconds,
                   scale_unit: :ppm, rates: nil, reference_epoch: nil,
                   rate_translation_unit: nil, rate_rotation_unit: nil,
                   rate_scale_unit: nil, name: nil)
      @name = name
      @convention = assert_convention(convention)

      @translation_unit = assert_unit(translation_unit, TRANSLATION_UNITS,
                                     'translation_unit')
      @rotation_unit = assert_unit(rotation_unit, ROTATION_UNITS,
                                  'rotation_unit')
      @scale_unit = assert_unit(scale_unit, SCALE_UNITS, 'scale_unit')

      @translation = assert_vector(translation, 'translation')
      @rotation = assert_vector(rotation, 'rotation')
      @scale = Geokit.assert_numeric(scale, 'scale')

      @rates = rates.nil? ? nil : prepare_rates(rates,
                                              rate_translation_unit,
                                              rate_rotation_unit,
                                              rate_scale_unit)
      @reference_epoch = prepare_reference_epoch(reference_epoch)
    end

    # True when parameter rates make the transformation epoch dependent.
    def time_dependent?
      !@rates.nil?
    end

    # Applies the forward transformation at +epoch+ (decimal years).
    def forward(point, epoch: nil)
      apply(matrix_at!(epoch), point)
    end

    # Applies the genuine inverse: the 4x4 forward matrix is formed at the
    # same epoch and inverted numerically. It is not built by flipping the
    # sign of every parameter (that is the exact inverse only for a zero
    # scale change).
    def backward(point, epoch: nil)
      apply(self.class.invert_matrix(matrix_at!(epoch)), point)
    end

    # Returns a Helmert wrapper whose forward is this numerical inverse.
    def inverse
      Inverse.new(self)
    end

    # Translation in metres at +epoch+.
    def translation_at(epoch)
      dt = delta_t(epoch)
      [0, 1, 2].map do |i|
        @translation[i] * @translation_unit +
          (@rates ? @rates[:translation][i] * dt : 0.0)
      end
    end

    # Rotations in radians at +epoch+.
    def rotation_at(epoch)
      dt = delta_t(epoch)
      [0, 1, 2].map do |i|
        @rotation[i] * @rotation_unit +
          (@rates ? @rates[:rotation][i] * dt : 0.0)
      end
    end

    # Dimensionless scale change at +epoch+ (e.g. 2.0e-5 for 20 ppm).
    def scale_at(epoch)
      @scale * @scale_unit +
        (@rates ? @rates[:scale] * delta_t(epoch) : 0.0)
    end

    # Builds the 4x4 homogeneous transformation matrix at +epoch+.
    def matrix_at(epoch)
      tx, ty, tz = translation_at(epoch)
      rx, ry, rz = rotation_at(epoch)
      rx, ry, rz = -rx, -ry, -rz if @convention == :coordinate_frame
      s = 1.0 + scale_at(epoch)

      [[s, s * -rz, s * ry, tx],
       [s * rz, s, s * -rx, ty],
       [s * -ry, s * rx, s, tz],
       [0.0, 0.0, 0.0, 1.0]]
    end

    # Inverts a 4x4 matrix using Gauss-Jordan elimination with partial
    # pivoting. Raises ArgumentError if the matrix is singular.
    def self.invert_matrix(matrix)
      n = matrix.size
      a = matrix.map.with_index do |row, i|
        row.map(&:to_f) + Array.new(n) { |j| i == j ? 1.0 : 0.0 }
      end

      n.times do |col|
        pivot = (col...n).max_by { |row| a[row][col].abs }
        if a[pivot][col].abs < 1.0e-15
          raise ArgumentError, 'Helmert matrix is singular; cannot invert'
        end
        a[col], a[pivot] = a[pivot], a[col] unless pivot == col

        pivot_value = a[col][col]
        a[col] = a[col].map { |value| value / pivot_value }
        n.times do |row|
          next if row == col
          factor = a[row][col]
          a[row] = a[row].map.with_index do |value, j|
            value - factor * a[col][j]
          end
        end
      end

      a.map { |row| row[n, n] }
    end

    # Transformation wrapper that applies the numerical inverse.
    class Inverse
      def initialize(helmert)
        @helmert = helmert
      end

      def forward(point, epoch: nil)
        @helmert.backward(point, epoch: epoch)
      end

      def backward(point, epoch: nil)
        @helmert.forward(point, epoch: epoch)
      end

      def time_dependent?
        @helmert.time_dependent?
      end

      def inverse
        @helmert
      end
    end

    private

    def matrix_at!(epoch)
      if @rates && epoch.nil?
        label = @name || 'Helmert'
        raise ArgumentError,
              "#{label} has parameter rates; an epoch (decimal years) " \
              'is required'
      end
      Geokit.assert_numeric(epoch, 'epoch') if @rates
      matrix_at(epoch)
    end

    def delta_t(epoch)
      return 0.0 unless @rates
      raise ArgumentError, 'an epoch (decimal years) is required' if epoch.nil?
      Geokit.assert_numeric(epoch, 'epoch') - @reference_epoch
    end

    def apply(matrix, point)
      x, y, z = point_coordinates(point)
      EcefPoint.new(
        matrix[0][0] * x + matrix[0][1] * y + matrix[0][2] * z + matrix[0][3],
        matrix[1][0] * x + matrix[1][1] * y + matrix[1][2] * z + matrix[1][3],
        matrix[2][0] * x + matrix[2][1] * y + matrix[2][2] * z + matrix[2][3]
      )
    end

    def point_coordinates(point)
      unless point.is_a?(EcefPoint)
        raise ArgumentError,
              "expected an EcefPoint, got #{point.class}; convert " \
              'geodetic coordinates with Datum#to_ecef first'
      end
      [point.x, point.y, point.z]
    end

    def assert_convention(value)
      unless value.is_a?(Symbol) && CONVENTIONS.include?(value)
        raise ArgumentError,
              'convention must be :position_vector or :coordinate_frame; ' \
              "got #{value.inspect}"
      end
      value
    end

    def assert_unit(value, units, label)
      unless value.is_a?(Symbol) && units.key?(value)
        raise ArgumentError,
              "#{label} must be one of #{units.keys.join(', ')}; " \
              "got #{value.inspect}"
      end
      units[value]
    end

    def assert_vector(value, label)
      unless value.respond_to?(:to_a) && value.to_a.size == 3
        raise ArgumentError, "#{label} must be an array of three components"
      end
      value.to_a.map { |component| Geokit.assert_numeric(component, label) }
    end

    def prepare_rates(rates, translation_rate_unit, rotation_rate_unit,
                    scale_rate_unit)
      unless rates.is_a?(Hash) && RATE_KEYS.all? { |key| rates.key?(key) }
        raise ArgumentError,
              "rates must include #{RATE_KEYS.join(', ')}"
      end
      @rate_translation_unit = assert_unit(translation_rate_unit,
                                          RATE_TRANSLATION_UNITS,
                                          'rate_translation_unit')
      @rate_rotation_unit = assert_unit(rotation_rate_unit,
                                       RATE_ROTATION_UNITS,
                                       'rate_rotation_unit')
      @rate_scale_unit = assert_unit(scale_rate_unit,
                                    RATE_SCALE_UNITS,
                                    'rate_scale_unit')

      prepared = {
        translation: assert_vector(rates[:translation], 'rates[:translation]'),
        rotation: assert_vector(rates[:rotation], 'rates[:rotation]'),
        scale: Geokit.assert_numeric(rates[:scale], 'rates[:scale]')
      }
      prepared[:translation].map! { |v| v * @rate_translation_unit }
      prepared[:rotation].map! { |v| v * @rate_rotation_unit }
      prepared[:scale] *= @rate_scale_unit
      prepared
    end

    def prepare_reference_epoch(value)
      if @rates && value.nil?
        raise ArgumentError,
              'reference_epoch (decimal years) is required when rates '
              'are supplied'
      end
      if @rates
        Geokit.assert_numeric(value, 'reference_epoch')
      elsif value
        raise ArgumentError,
              'reference_epoch is meaningless without a complete rates hash; '
              'omit it or supply translation, rotation and scale rates'
      end
      value
    end

    # Ordnance Survey Great Britain 1936 (Airy 1830) to ETRS89 / WGS 84,
    # position-vector convention. The EPSG:1314 registry lists the inverse
    # direction (-446.448, ...); these values are its sign-negated forward
    # parameters, which map OSGB36 ECEF to WGS 84 ECEF. A static set.
    OSGB36_TO_WGS84 = new(
      name: 'OSGB36 to WGS 84 (EPSG:1314)',
      convention: :position_vector,
      translation: [446.448, -125.157, 542.060],
      rotation: [0.1502, 0.2470, 0.8421],
      scale: 20.4894,
      translation_unit: :metres,
      rotation_unit: :arcseconds,
      scale_unit: :ppm
    )
  end
end
