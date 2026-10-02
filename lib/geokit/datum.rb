module Geokit
  # A geodetic datum: a named reference Ellipsoid plus (optionally) a known
  # seven-parameter Helmert transformation linking it to a pivot datum.
  #
  # Datum objects are explicit: a plain LatLng never assumes one. They provide
  # geodetic (latitude/longitude/height) <-> ECEF (Earth-Centred,
  # Earth-Fixed XYZ) conversions on the datum ellipsoid. Transformations
  # between datums are performed by Helmert objects supplied by the caller.
  class Datum
    RAD = Math::PI / 180.0

    # Absolute tolerance for the ECEF -> geodetic iteration, in metres.
    HEIGHT_TOLERANCE = 1.0e-11
    # Cap on Bowring iterations (normally converges in 3-4 passes).
    MAX_ITERATIONS = 100

    attr_reader :name, :ellipsoid

    def initialize(name, ellipsoid)
      raise ArgumentError, 'ellipsoid is required' if ellipsoid.nil?
      unless ellipsoid.is_a?(Ellipsoid)
        raise ArgumentError, "expected an Ellipsoid, got #{ellipsoid.class}"
      end
      @name = name
      @ellipsoid = ellipsoid
    end

    # Converts a GeodeticPoint (or lat in degrees, lng in degrees and height
    # in metres) to an EcefPoint on this datum.
    def to_ecef(point_or_lat, lng = nil, height = nil)
      lat_deg, lng_deg, h =
        geodetic_arguments(point_or_lat, lng, height)

      lat = lat_deg * RAD
      lng = lng_deg * RAD
      sin_lat = Math.sin(lat)
      cos_lat = Math.cos(lat)
      prime_vertical = prime_vertical_radius(sin_lat)

      x = (prime_vertical + h) * cos_lat * Math.cos(lng)
      y = (prime_vertical + h) * cos_lat * Math.sin(lng)
      z = (prime_vertical * (1.0 - ellipsoid.first_eccentricity_squared) + h) *
          sin_lat
      EcefPoint.new(x, y, z)
    end

    # Converts an EcefPoint (metres) back to a GeodeticPoint on this datum
    # using Bowring's iterative formula. The iteration is robust at the poles
    # (the meridian plane is used when the planar radius vanishes) and for
    # negative ellipsoidal heights.
    def from_ecef(point_or_x, y = nil, z = nil)
      x, y, z = ecef_arguments(point_or_x, y, z)
      a = ellipsoid.semi_major_axis
      e2 = ellipsoid.first_eccentricity_squared

      planar_radius = Math.sqrt(x * x + y * y)

      if planar_radius < a * 1.0e-14
        lng_deg = 0.0
        z_sign = z.zero? ? 0.0 : (z > 0.0 ? 1.0 : -1.0)
        return GeodeticPoint.new(90.0 * z_sign, 0.0,
                                z.abs / Math.sqrt(1.0 - e2) - a,
                                datum: self)
      end

      lng_deg = Math.atan2(y, x) / RAD
      lat = Math.atan2(z, planar_radius * (1.0 - e2))
      height = 0.0

      MAX_ITERATIONS.times do
        sin_lat = Math.sin(lat)
        prime_vertical = prime_vertical_radius(sin_lat)
        new_height = planar_radius / Math.cos(lat) - prime_vertical
        lat = Math.atan2(z, planar_radius *
                             (1.0 - e2 * prime_vertical /
                                      (prime_vertical + new_height)))
        if (new_height - height).abs <= HEIGHT_TOLERANCE
          height = new_height
          break
        end
        height = new_height
      end

      GeodeticPoint.new(lat / RAD, lng_deg, height, datum: self)
    end

    # Transforms a geodetic point on this datum to a geodetic point on
    # +target_datum+ via +helmert+ (a Helmert expressed from self toward the
    # target). +epoch+ is required for time-dependent Helmert sets.
    def transform(point, helmert, target_datum, epoch: nil)
      source_ecef(point).transform(helmert, epoch: epoch)
                   .to_geodetic(target_datum)
    end

    def to_s
      "#<#{self.class.name} #{name}>"
    end

    def ==(other)
      other.is_a?(Datum) && name == other.name && ellipsoid == other.ellipsoid
    end
    alias eql? ==

    def hash
      [name, ellipsoid].hash
    end

    WGS84 = new('WGS 84', Ellipsoid::WGS84)
    # ETRS89 shares the GRS80 ellipsoid and is fixed to the stable part of
    # the Eurasian plate.
    ETRS89 = new('ETRS89', Ellipsoid::GRS80)
    # NAD83 uses the GRS80 ellipsoid as well.
    NAD83 = new('NAD83', Ellipsoid::GRS80)
    # Great Britain's historic horizontal datum (Airy 1830 ellipsoid).
    OSGB36 = new('OSGB36', Ellipsoid::AIRY_1830)
    # Many European mapping products (e.g. western Europe before ETRS89).
    ED50 = new('ED50', Ellipsoid::INTERNATIONAL_1924)

    private

    # Radius of curvature in the prime vertical, N = a / sqrt(1 - e^2 sin^2).
    def prime_vertical_radius(sin_lat)
      ellipsoid.semi_major_axis /
        Math.sqrt(1.0 - ellipsoid.first_eccentricity_squared *
                          sin_lat * sin_lat)
    end

    def source_ecef(point)
      return to_ecef(point) unless point.is_a?(GeodeticPoint) && point.datum
      unless point.datum == self
        raise ArgumentError,
              "point is on #{point.datum.name}, expected #{name}"
      end
      to_ecef(point)
    end

    def geodetic_arguments(point_or_lat, lng, height)
      if point_or_lat.is_a?(GeodeticPoint)
        unless lng.nil? && height.nil?
          raise ArgumentError,
                'pass a GeodeticPoint or lat/lng/height, not both'
        end
        [point_or_lat.lat, point_or_lat.lng, point_or_lat.height]
      else
        h = height.nil? ? 0.0 : Geokit.assert_numeric(height, 'height')
        [Geokit.assert_numeric(point_or_lat, 'lat'),
         Geokit.assert_numeric(lng, 'lng'), h]
      end
    end

    def ecef_arguments(point_or_x, y, z)
      if point_or_x.is_a?(EcefPoint)
        unless y.nil? && z.nil?
          raise ArgumentError,
                'pass an EcefPoint or x/y/z, not both'
        end
        [point_or_x.x, point_or_x.y, point_or_x.z]
      else
        [Geokit.assert_numeric(point_or_x, 'x'),
         Geokit.assert_numeric(y, 'y'),
         Geokit.assert_numeric(z, 'z')]
      end
    end
  end
end
