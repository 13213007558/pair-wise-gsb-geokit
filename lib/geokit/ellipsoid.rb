module Geokit
  # Geodesic calculations on a reference ellipsoid (oblate spheroid) using
  # Vincenty's formulae ("Direct and Inverse Solutions of Geodesics on the
  # Ellipsoid with Application of Nested Equations", Survey Review, 1975).
  #
  # This provides higher accuracy than the spherical/flat formulas in
  # Geokit::Mappable (which are left unchanged):
  #
  # * +inverse+ -- shortest distance between two points plus the azimuths
  #   (headings) at both endpoints.
  # * +direct+ -- endpoint and final azimuth given a start point, an initial
  #   azimuth and a distance.
  #
  # An ellipsoid is defined by its semi-major axis +a+ (in meters) and its
  # flattening +f+. Geokit::Ellipsoid::WGS84 is provided as a ready-made
  # preset and is the default used by the :ellipsoid formula in
  # Geokit::Mappable. Custom ellipsoids can be built for any reasonable
  # oblate spheroid (+a+ > 0, 0 <= f < 1; geodetic ellipsoids have
  # f ~ 1/300).
  #
  # Conventions:
  #
  # * Distances are computed in meters and converted with the same unit
  #   configuration as the rest of Geokit (:miles, :kms, :meters, :nms;
  #   Geokit.default_units when not specified). Both inverse and direct use
  #   the same unit handling.
  # * Azimuths are returned in degrees, clockwise from north, in [0, 360).
  # * Result longitudes are normalized to [-180, 180); longitude differences
  #   are wrapped across the antimeridian, so (10, 179.5) to (10, -179.5)
  #   crosses the Pacific, not the Atlantic.
  # * Coincident points (after longitude normalization) return a zero
  #   distance and zero azimuths; a zero-distance direct solution returns
  #   the start point.
  # * Antipodal (or nearly antipodal) points can have many equally short
  #   geodesics; Vincenty's iteration does not converge there. Instead of
  #   silently falling back to a spherical approximation, a
  #   Geokit::Ellipsoid::NotConvergedError is raised and the caller decides
  #   how to handle the ambiguity.
  #
  # No network access or geocoding service is used.
  class Ellipsoid
    # Raised when Vincenty's iteration fails to converge (typically for
    # (nearly) antipodal points where multiple equally short geodesics
    # exist). Never rescued internally: no spherical fallback is returned.
    class NotConvergedError < StandardError; end

    InverseResult = Struct.new(:distance, :initial_azimuth, :final_azimuth)
    DirectResult  = Struct.new(:endpoint, :final_azimuth)

    MAX_ITERATIONS = 200
    CONVERGENCE_THRESHOLD = 1e-12

    attr_reader :semi_major_axis, :flattening, :semi_minor_axis

    def initialize(semi_major_axis, flattening)
      a = Float(semi_major_axis)
      f = Float(flattening)
      unless a.positive? && a.finite?
        raise ArgumentError, "semi-major axis must be a positive finite number, got #{semi_major_axis.inspect}"
      end
      unless f >= 0 && f < 1
        raise ArgumentError, "flattening must be in [0, 1) for an oblate spheroid, got #{flattening.inspect}"
      end
      @semi_major_axis = a
      @flattening = f
      @semi_minor_axis = a * (1 - f)
    end

    # Inverse problem: shortest geodesic distance between +from+ and +to+
    # (Geokit::LatLng-compatible values) plus the azimuths at both ends.
    #
    # Returns an InverseResult with :distance (in +units+),
    # :initial_azimuth and :final_azimuth (degrees, [0, 360)).
    #
    # Raises NotConvergedError when no unique shortest geodesic can be
    # computed (e.g. antipodal points).
    def inverse(from, to, units: nil)
      units = self.class.normalize_units(units)
      from = Geokit::LatLng.normalize(from)
      to = Geokit::LatLng.normalize(to)

      phi1 = deg2rad(from.lat)
      phi2 = deg2rad(to.lat)
      l = deg2rad(self.class.normalize_lng(to.lng) - self.class.normalize_lng(from.lng))

      return InverseResult.new(0.0, 0.0, 0.0) if phi1 == phi2 && l.zero?

      f = flattening
      tan_u1 = (1 - f) * Math.sin(phi1) / Math.cos(phi1)
      tan_u2 = (1 - f) * Math.sin(phi2) / Math.cos(phi2)
      u1 = Math.atan(tan_u1)
      u2 = Math.atan(tan_u2)
      sin_u1 = Math.sin(u1)
      cos_u1 = Math.cos(u1)
      sin_u2 = Math.sin(u2)
      cos_u2 = Math.cos(u2)

      lambda = l
      sin_sigma = cos_sigma = sigma = sin_alpha = cos_sq_alpha = cos2_sigma_m = nil
      converged = false
      MAX_ITERATIONS.times do
        sin_lambda = Math.sin(lambda)
        cos_lambda = Math.cos(lambda)
        sin_sigma = Math.sqrt(
          (cos_u2 * sin_lambda)**2 +
          (cos_u1 * sin_u2 - sin_u1 * cos_u2 * cos_lambda)**2
        )
        # Coincident points (defensive; handled above after normalization).
        return InverseResult.new(0.0, 0.0, 0.0) if sin_sigma.zero?

        cos_sigma = sin_u1 * sin_u2 + cos_u1 * cos_u2 * cos_lambda
        sigma = Math.atan2(sin_sigma, cos_sigma)
        sin_alpha = cos_u1 * cos_u2 * sin_lambda / sin_sigma
        cos_sq_alpha = 1 - sin_alpha**2
        cos2_sigma_m = cos_sq_alpha.zero? ? 0.0 : cos_sigma - 2 * sin_u1 * sin_u2 / cos_sq_alpha
        c = f / 16 * cos_sq_alpha * (4 + f * (4 - 3 * cos_sq_alpha))
        previous = lambda
        lambda = l + (1 - c) * f * sin_alpha *
                 (sigma + c * sin_sigma *
                   (cos2_sigma_m + c * cos_sigma * (-1 + 2 * cos2_sigma_m**2)))
        if (lambda - previous).abs <= CONVERGENCE_THRESHOLD
          converged = true
          break
        end
      end

      unless converged
        raise NotConvergedError,
              "geodesic between #{from.ll} and #{to.ll} did not converge " +              "(points are (nearly) antipodal and have no unique shortest path)"
      end

      u_sq = cos_sq_alpha * (semi_major_axis**2 - semi_minor_axis**2) / semi_minor_axis**2
      big_a = 1 + u_sq / 16_384 * (4096 + u_sq * (-768 + u_sq * (320 - 175 * u_sq)))
      big_b = u_sq / 1024 * (256 + u_sq * (-128 + u_sq * (74 - 47 * u_sq)))
      delta_sigma = big_b * sin_sigma *
                    (cos2_sigma_m + big_b / 4 *
                      (cos_sigma * (-1 + 2 * cos2_sigma_m**2) -
                       big_b / 6 * cos2_sigma_m * (-3 + 4 * sin_sigma**2) * (-3 + 4 * cos2_sigma_m**2)))
      meters = semi_minor_axis * big_a * (sigma - delta_sigma)

      sin_lambda = Math.sin(lambda)
      cos_lambda = Math.cos(lambda)
      initial_azimuth = Math.atan2(cos_u2 * sin_lambda,
                                   cos_u1 * sin_u2 - sin_u1 * cos_u2 * cos_lambda)
      final_azimuth = Math.atan2(cos_u1 * sin_lambda,
                                 -sin_u1 * cos_u2 + cos_u1 * sin_u2 * cos_lambda)

      InverseResult.new(
        self.class.meters_to_units(meters, units),
        normalize_azimuth(rad2deg(initial_azimuth)),
        normalize_azimuth(rad2deg(final_azimuth))
      )
    end

    # Direct problem: given a start point, an initial azimuth (degrees,
    # clockwise from north) and a distance (in +units+), returns a
    # DirectResult with the :endpoint (a Geokit::LatLng, longitude
    # normalized to [-180, 180)) and the :final_azimuth (degrees, [0, 360)).
    def direct(from, azimuth, distance, units: nil)
      units = self.class.normalize_units(units)
      from = Geokit::LatLng.normalize(from)
      meters = self.class.units_to_meters(distance.to_f, units)
      alpha1 = deg2rad(azimuth.to_f)

      if meters.zero?
        return DirectResult.new(
          LatLng.new(from.lat, self.class.normalize_lng(from.lng)),
          normalize_azimuth(azimuth.to_f)
        )
      end

      f = flattening
      phi1 = deg2rad(from.lat)
      tan_u1 = (1 - f) * Math.sin(phi1) / Math.cos(phi1)
      u1 = Math.atan(tan_u1)
      sin_u1 = Math.sin(u1)
      cos_u1 = Math.cos(u1)
      sin_alpha1 = Math.sin(alpha1)
      cos_alpha1 = Math.cos(alpha1)

      sigma1 = Math.atan2(tan_u1, cos_alpha1)
      sin_alpha = cos_u1 * sin_alpha1
      cos_sq_alpha = 1 - sin_alpha**2
      u_sq = cos_sq_alpha * (semi_major_axis**2 - semi_minor_axis**2) / semi_minor_axis**2
      big_a = 1 + u_sq / 16_384 * (4096 + u_sq * (-768 + u_sq * (320 - 175 * u_sq)))
      big_b = u_sq / 1024 * (256 + u_sq * (-128 + u_sq * (74 - 47 * u_sq)))

      sigma = meters / (semi_minor_axis * big_a)
      converged = false
      sin_sigma = cos_sigma = cos2_sigma_m = delta_sigma = nil
      MAX_ITERATIONS.times do
        cos2_sigma_m = Math.cos(2 * sigma1 + sigma)
        sin_sigma = Math.sin(sigma)
        cos_sigma = Math.cos(sigma)
        delta_sigma = big_b * sin_sigma *
                      (cos2_sigma_m + big_b / 4 *
                        (cos_sigma * (-1 + 2 * cos2_sigma_m**2) -
                         big_b / 6 * cos2_sigma_m * (-3 + 4 * sin_sigma**2) * (-3 + 4 * cos2_sigma_m**2)))
        previous = sigma
        sigma = meters / (semi_minor_axis * big_a) + delta_sigma
        if (sigma - previous).abs <= CONVERGENCE_THRESHOLD
          converged = true
          break
        end
      end

      unless converged
        raise NotConvergedError,
              "direct geodesic from #{from.ll} did not converge"
      end

      tmp = sin_u1 * sin_sigma - cos_u1 * cos_sigma * cos_alpha1
      phi2 = Math.atan2(sin_u1 * cos_sigma + cos_u1 * sin_sigma * cos_alpha1,
                        (1 - f) * Math.sqrt(sin_alpha**2 + tmp**2))
      lambda = Math.atan2(sin_sigma * sin_alpha1,
                          cos_u1 * cos_sigma - sin_u1 * sin_sigma * cos_alpha1)
      c = f / 16 * cos_sq_alpha * (4 + f * (4 - 3 * cos_sq_alpha))
      l = lambda - (1 - c) * f * sin_alpha *
          (sigma + c * sin_sigma *
            (cos2_sigma_m + c * cos_sigma * (-1 + 2 * cos2_sigma_m**2)))
      lng2 = deg2rad(from.lng) + l
      alpha2 = Math.atan2(sin_alpha, -tmp)

      DirectResult.new(
        LatLng.new(rad2deg(phi2), self.class.normalize_lng(rad2deg(lng2))),
        normalize_azimuth(rad2deg(alpha2))
      )
    end

    # Validates and defaults the units option exactly like Geokit::Mappable.
    def self.normalize_units(units)
      units = Geokit.default_units if units.nil?
      unless Geokit::Mappable::ClassMethods::EARTH_RADIUS.key?(units)
        raise ArgumentError, "#{units} is an unsupported unit of length."
      end
      units
    end

    # Converts meters to the given unit using the same registry as
    # Geokit::Mappable, so ellipsoid and sphere/flat share one unit config.
    def self.meters_to_units(meters, units)
      registry = Geokit::Mappable::ClassMethods
      meters * registry::EARTH_RADIUS[units] / registry::EARTH_RADIUS_IN_METERS
    end

    def self.units_to_meters(distance, units)
      registry = Geokit::Mappable::ClassMethods
      distance * registry::EARTH_RADIUS_IN_METERS / registry::EARTH_RADIUS[units]
    end

    # Normalizes a longitude in degrees to [-180, 180).
    def self.normalize_lng(lng)
      ((lng.to_f + 180) % 360) - 180
    end

    private

    def deg2rad(degrees)
      degrees.to_f / 180.0 * Math::PI
    end

    def rad2deg(rad)
      rad.to_f * 180.0 / Math::PI
    end

    def normalize_azimuth(degrees)
      normalized = degrees % 360
      # Guard against tiny negative azimuths rounding up to exactly 360.0.
      normalized >= 360.0 ? 0.0 : normalized
    end

    # WGS84 ellipsoid preset (semi-major axis and flattening as defined by
    # the WGS84 datum).
    WGS84 = new(6_378_137.0, 1 / 298.257223563).freeze
  end
end
