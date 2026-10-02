module Geokit
  # A Datum couples a reference ellipsoid with an optional time-dependent
  # 7-parameter Helmert transformation to WGS 84. It provides geodetic
  # (latitude/longitude/height) to ECEF (earth-centered, earth-fixed)
  # conversion and datum-to-datum transformation of coordinates.
  #
  # Datum handling is always explicit: Geokit::LatLng carries no datum and
  # no online transformation service is consulted. All math is local.
  #
  # Example:
  #   osgb36 = Geokit::Datum::DATUMS[:osgb36]
  #   wgs84  = Geokit::Datum::DATUMS[:wgs84]
  #   lat, lng, h = osgb36.convert_to(wgs84, 52.65757, 1.71792, 53.0)
  class Datum
    # A reference ellipsoid defined by semi-major axis +a+ (meters) and
    # flattening +f+.
    class Ellipsoid
      attr_reader :name, :a, :f

      def initialize(name:, a:, f:)
        @name = name
        @a = Float(a)
        @f = Float(f)
      end

      # Semi-minor axis in meters.
      def b
        a * (1.0 - f)
      end

      # First eccentricity squared.
      def eccentricity_squared
        f * (2.0 - f)
      end
    end

    ELLIPSOIDS = {
      wgs84: Ellipsoid.new(name: 'WGS 84', a: 6_378_137.0, f: 1.0 / 298.257223563),
      grs80: Ellipsoid.new(name: 'GRS 80', a: 6_378_137.0, f: 1.0 / 298.257222101),
      airy1830: Ellipsoid.new(name: 'Airy 1830', a: 6_377_563.396, f: 1.0 / 299.3249646),
      international1924: Ellipsoid.new(name: 'International 1924', a: 6_378_388.0, f: 1.0 / 297.0)
    }.freeze

    # A 7-parameter Helmert (similarity) transformation between terrestrial
    # reference frames, optionally with per-year parameter rates for
    # reference-epoch propagation (a 14-parameter transformation).
    #
    # The rotation convention must be stated explicitly:
    #   :position_vector  -- EPSG "Position Vector transformation" (9606),
    #                        used by ITRF/ETRS transformations.
    #   :coordinate_frame -- EPSG "Coordinate Frame rotation" (9607), the
    #                        rotations are the negatives of position_vector.
    #
    # Units must be stated explicitly; nothing is guessed:
    #   translation_unit: :meters or :millimeters
    #   rotation_unit:    :arc_seconds, :degrees or :radians
    #   scale_unit:       :ppm or :unity
    #
    # When +rates+ (per-year parameter changes, in the same units) are
    # supplied, +reference_epoch+ (decimal year) is required and every
    # forward/inverse call must pass +epoch:+; a missing epoch raises
    # ArgumentError rather than guessing one.
    class HelmertTransform
      CONVENTIONS = %i[position_vector coordinate_frame].freeze
      TRANSLATION_UNITS = { meters: 1.0, millimeters: 1.0e-3 }.freeze
      ROTATION_UNITS = {
        arc_seconds: Math::PI / (180.0 * 3600.0),
        degrees: Math::PI / 180.0,
        radians: 1.0
      }.freeze
      SCALE_UNITS = { ppm: 1.0e-6, unity: 1.0 }.freeze
      RATE_KEYS = %i[dtx dty dtz drx dry drz dscale].freeze

      attr_reader :convention, :reference_epoch

      def initialize(tx:, ty:, tz:, rx:, ry:, rz:, scale:,
                     convention:, translation_unit:, rotation_unit:, scale_unit:,
                     rates: nil, reference_epoch: nil)
        @convention = validate_choice(:convention, convention, CONVENTIONS)
        translation_factor = validate_choice(:translation_unit, translation_unit, TRANSLATION_UNITS)
        rotation_factor = validate_choice(:rotation_unit, rotation_unit, ROTATION_UNITS)
        scale_factor = validate_choice(:scale_unit, scale_unit, SCALE_UNITS)

        @translation = [tx, ty, tz].map { |v| Float(v) * translation_factor }
        @rotation = [rx, ry, rz].map { |v| Float(v) * rotation_factor }
        @scale = Float(scale) * scale_factor

        if rates
          unknown = rates.keys.map(&:to_sym) - RATE_KEYS
          unless unknown.empty?
            raise ArgumentError, "unknown rate parameters: #{unknown.join(', ')}"
          end
          if reference_epoch.nil?
            raise ArgumentError, 'reference_epoch is required when parameter rates are given'
          end

          @rates = {
            translation: %i[dtx dty dtz].map { |k| Float(rates.fetch(k, 0.0)) * translation_factor },
            rotation: %i[drx dry drz].map { |k| Float(rates.fetch(k, 0.0)) * rotation_factor },
            scale: Float(rates.fetch(:dscale, 0.0)) * scale_factor
          }
          @reference_epoch = Float(reference_epoch)
        elsif !reference_epoch.nil?
          raise ArgumentError, 'reference_epoch was given but no parameter rates were given'
        end
      end

      def time_dependent?
        !@rates.nil?
      end

      # Applies the transformation to an ECEF [x, y, z] point (meters).
      def forward(point, epoch: nil)
        matrix, translation = forward_matrix(epoch)
        x, y, z = point.map { |v| Float(v) }
        [
          matrix[0][0] * x + matrix[0][1] * y + matrix[0][2] * z + translation[0],
          matrix[1][0] * x + matrix[1][1] * y + matrix[1][2] * z + translation[1],
          matrix[2][0] * x + matrix[2][1] * y + matrix[2][2] * z + translation[2]
        ]
      end

      # Applies the exact inverse of the forward transformation: the forward
      # linear map M = (1 + s) * R is inverted as a 3x3 matrix, so
      # inverse(forward(p)) == p to numerical precision. The inverse is NOT
      # computed by negating the parameters, which would only be a
      # first-order approximation.
      def inverse(point, epoch: nil)
        matrix, translation = forward_matrix(epoch)
        inv = invert3(matrix)
        x = Float(point[0]) - translation[0]
        y = Float(point[1]) - translation[1]
        z = Float(point[2]) - translation[2]
        [
          inv[0][0] * x + inv[0][1] * y + inv[0][2] * z,
          inv[1][0] * x + inv[1][1] * y + inv[1][2] * z,
          inv[2][0] * x + inv[2][1] * y + inv[2][2] * z
        ]
      end

      # The [matrix, translation] of the forward map at the given epoch.
      # Exposed so callers can verify invertibility directly.
      def forward_matrix(epoch = nil)
        params = parameters_at(epoch)
        sign = @convention == :position_vector ? 1.0 : -1.0
        rx, ry, rz = params[:rotation].map { |r| r * sign }
        k = 1.0 + params[:scale]
        matrix = [
          [k, -rz * k, ry * k],
          [rz * k, k, -rx * k],
          [-ry * k, rx * k, k]
        ]
        [matrix, params[:translation]]
      end

      private

      def validate_choice(name, value, allowed)
        key = value.to_sym
        return allowed[key] if allowed.is_a?(Hash) && allowed.key?(key)
        return key if allowed.is_a?(Array) && allowed.include?(key)

        list = allowed.is_a?(Array) ? allowed : allowed.keys
        raise ArgumentError,
              "#{name} must be one of #{list.join(', ')} (got #{value.inspect}); " \
              'units and conventions are never guessed'
      end

      def parameters_at(epoch)
        return { translation: @translation, rotation: @rotation, scale: @scale } unless @rates

        if epoch.nil?
          raise ArgumentError,
                'epoch is required for a time-dependent Helmert transform; refusing to guess one'
        end

        dt = Float(epoch) - @reference_epoch
        {
          translation: @translation.zip(@rates[:translation]).map { |p, rate| p + rate * dt },
          rotation: @rotation.zip(@rates[:rotation]).map { |p, rate| p + rate * dt },
          scale: @scale + @rates[:scale] * dt
        }
      end

      # Exact inverse of a 3x3 matrix via the adjugate over the determinant.
      def invert3(m)
        cofactors = [
          [m[1][1] * m[2][2] - m[1][2] * m[2][1], m[0][2] * m[2][1] - m[0][1] * m[2][2], m[0][1] * m[1][2] - m[0][2] * m[1][1]],
          [m[1][2] * m[2][0] - m[1][0] * m[2][2], m[0][0] * m[2][2] - m[0][2] * m[2][0], m[0][2] * m[1][0] - m[0][0] * m[1][2]],
          [m[1][0] * m[2][1] - m[1][1] * m[2][0], m[0][1] * m[2][0] - m[0][0] * m[2][1], m[0][0] * m[1][1] - m[0][1] * m[1][0]]
        ]
        det = m[0][0] * cofactors[0][0] + m[0][1] * cofactors[1][0] + m[0][2] * cofactors[2][0]
        raise ArgumentError, 'Helmert matrix is singular and cannot be inverted' if det.zero?

        cofactors.map { |row| row.map { |c| c / det } }
      end
    end

    attr_reader :name, :ellipsoid, :to_wgs84

    def initialize(name:, ellipsoid:, to_wgs84: nil)
      @name = name
      @ellipsoid = ellipsoid
      @to_wgs84 = to_wgs84
    end

    # Geodetic to ECEF. +lat+ and +lng+ are in degrees, +height+ is
    # ellipsoidal height in meters. Returns [x, y, z] in meters.
    def to_ecef(lat, lng, height)
      lat_rad = Float(lat) * Math::PI / 180.0
      lng_rad = Float(lng) * Math::PI / 180.0
      h = Float(height)
      a = ellipsoid.a
      e2 = ellipsoid.eccentricity_squared
      sin_lat = Math.sin(lat_rad)
      cos_lat = Math.cos(lat_rad)
      n = a / Math.sqrt(1.0 - e2 * sin_lat * sin_lat)
      [
        (n + h) * cos_lat * Math.cos(lng_rad),
        (n + h) * cos_lat * Math.sin(lng_rad),
        (n * (1.0 - e2) + h) * sin_lat
      ]
    end

    # ECEF to geodetic. Returns [lat, lng, height] with lat/lng in degrees
    # and height in meters. Bowring's closed form plus fixed-point
    # refinement; stable near the poles and for negative heights.
    def from_ecef(x, y, z)
      x = Float(x)
      y = Float(y)
      z = Float(z)
      a = ellipsoid.a
      b = ellipsoid.b
      e2 = ellipsoid.eccentricity_squared
      p = Math.sqrt(x * x + y * y)
      lng = Math.atan2(y, x)

      if p < 1.0e-12
        lat = z.negative? ? -Math::PI / 2.0 : Math::PI / 2.0
        return [lat * 180.0 / Math::PI, lng * 180.0 / Math::PI, z.abs - b]
      end

      ep2 = (a * a - b * b) / (b * b)
      theta = Math.atan2(z * a, p * b)
      sin_t = Math.sin(theta)
      cos_t = Math.cos(theta)
      lat = Math.atan2(z + ep2 * b * sin_t**3, p - e2 * a * cos_t**3)

      height = 0.0
      5.times do
        n = a / Math.sqrt(1.0 - e2 * Math.sin(lat)**2)
        height = p / Math.cos(lat) - n
        lat = Math.atan2(z, p * (1.0 - e2 * n / (n + height)))
      end
      n = a / Math.sqrt(1.0 - e2 * Math.sin(lat)**2)
      height = p / Math.cos(lat) - n

      [lat * 180.0 / Math::PI, lng * 180.0 / Math::PI, height]
    end

    # Transforms geodetic coordinates (degrees, meters) from this datum to
    # +target_datum+ via ECEF and each datum's Helmert transformation to
    # WGS 84. Returns [lat, lng, height] on the target datum. +epoch+ is
    # required whenever either datum's transformation is time-dependent.
    def convert_to(target_datum, lat, lng, height, epoch: nil)
      ecef = to_ecef(lat, lng, height)
      ecef = to_wgs84.forward(ecef, epoch: epoch) if to_wgs84
      ecef = target_datum.to_wgs84.inverse(ecef, epoch: epoch) if target_datum.to_wgs84
      target_datum.from_ecef(*ecef)
    end

    DATUMS = {
      wgs84: new(name: 'WGS 84', ellipsoid: ELLIPSOIDS[:wgs84]),
      osgb36: new(
        name: 'OSGB 36',
        ellipsoid: ELLIPSOIDS[:airy1830],
        # Ordnance Survey, "A guide to coordinate systems in Great Britain":
        # OSGB36 -> ETRS89, position-vector convention.
        to_wgs84: HelmertTransform.new(
          tx: 446.448, ty: -125.157, tz: 542.060,
          rx: 0.1502, ry: 0.2470, rz: 0.8421, scale: -20.4894,
          convention: :position_vector,
          translation_unit: :meters, rotation_unit: :arc_seconds, scale_unit: :ppm
        )
      )
    }.freeze
  end
end
