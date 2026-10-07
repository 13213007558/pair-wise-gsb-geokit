require File.join(File.dirname(__FILE__), 'helper')

# Tests for Geokit::Ellipsoid (Vincenty geodesics on an oblate spheroid).
#
# Reference vectors:
# * The Bessel-ellipsoid example worked in Vincenty's 1975 paper
#   (Survey Review 23(176):88-93), the canonical published check.
# * WGS84 vectors computed with a 50-digit (mpmath) implementation of
#   Vincenty's inverse formula and cross-checked by an independent
#   high-precision numerical integration of the geodesic ODE (DOP853,
#   rtol 1e-13): agreement better than 5e-5 m in position and 4e-11 deg
#   in final azimuth for every vector. The two exactly polar vectors were
#   additionally verified against the analytic meridian arc
#   (agreement < 1e-6 m).
#
# Observed error of this implementation against the references:
# distance <= 3.6e-6 m, azimuth <= 2.4e-9 deg on all vectors.
# Test tolerances below (1e-4 m, 1e-6 deg) leave a wide margin while
# still proving sub-millimetre geodesic accuracy.
class EllipsoidTest < Test::Unit::TestCase #:nodoc: all
  WGS84 = Geokit::Ellipsoid::WGS84

  DISTANCE_TOLERANCE_METERS = 1e-4 # 0.1 mm
  AZIMUTH_TOLERANCE_DEG = 1e-6
  LATLNG_TOLERANCE_DEG = 1e-8 # ~1 mm

  # from(lat,lng), to(lat,lng), distance(m), initial azimuth, final azimuth
  WGS84_VECTORS = [
    # ordinary equatorial region
    [[0, 0], [10, 20], 2475724.5187498556, 62.879291486367884, 64.647060901043973],
    # London to New York
    [[51.5074, -0.1278], [40.7128, -74.0060], 5585233.5789467786, 288.36899811332881, 231.24122911945191],
    # Sydney to Tokyo
    [[-33.8688, 151.2093], [35.6762, 139.6503], 7792174.8271172564, 349.99751778795417, 349.77466835910226],
    # short hop (San Francisco to Oakland)
    [[37.7749, -122.4194], [37.8044, -122.2711], 13466.783133967620, 75.882786031662420, 75.973659013495645],
    # near antipodal (0.5 deg off exact antipodes)
    [[0, 0], [0.5, 179.5], 19936288.5789786176, 25.671872853894430, 154.32708548433972],
    # near antipodal (0.1 deg off exact antipodes)
    [[10, 20], [-9.9, -159.0], 19925527.3944032240, 283.41085492925382, 256.51625368853857],
    # near antipodal (southern hemisphere)
    [[-40, 50], [39.5, -129.5], 19939061.0930086128, 205.17932113423258, 335.01454805027059],
    # near pole, same latitude
    [[89.0, 0], [89.0, 120], 193456.99325045545, 30.003779080548825, 149.99622091945118],
    # exactly over the north pole
    [[85.0, 0], [85.0, 180], 1116911.1772935498, 0.0, 180.0],
    # exactly over the pole, oblique longitudes
    [[88.5, -30], [86.2, 150], 591971.41924970577, 0.0, 180.0],
    # across the antimeridian, short
    [[10, 179.5], [10, -179.5], 109639.32210565462, 89.913173759580020, 90.086826240419980],
    # across the antimeridian, long
    [[-20, 170], [30, -160], 6398977.3365716170, 30.955295858579059, 33.908928793048440]
  ].freeze

  # Exact analytical WGS84 values.
  EQUATORIAL_QUADRANT_METERS = 6378137.0 * Math::PI / 2 # a * pi/2, exact
  MERIDIAN_QUADRANT_METERS = 10001965.7293127228 # published WGS84 quarter meridian

  # Vincenty's 1975 paper example (Bessel 1841 ellipsoid).
  BESSEL = Geokit::Ellipsoid.new(6377397.155, 1 / 299.1528128)
  PAPER_FROM = Geokit::LatLng.new(55.75, 0)            # 55 45 00 N, 0 E
  PAPER_TO = Geokit::LatLng.new(-33.4333333333, 108.2166666667) # 33 26 S, 108 13 E
  PAPER_DISTANCE = 14110526.170                        # meters
  PAPER_INITIAL_AZIMUTH = 96.6024443333                # 96 36 08.7996
  PAPER_FINAL_AZIMUTH = 137.8727818056                 # 137 52 22.0145

  def test_vincenty_paper_inverse
    res = BESSEL.inverse(PAPER_FROM, PAPER_TO, units: :meters)
    assert_in_delta PAPER_DISTANCE, res.distance, 0.001
    assert_in_delta PAPER_INITIAL_AZIMUTH, res.initial_azimuth, 1e-6
    assert_in_delta PAPER_FINAL_AZIMUTH, res.final_azimuth, 1e-6
  end

  def test_vincenty_paper_direct
    res = BESSEL.direct(PAPER_FROM, PAPER_INITIAL_AZIMUTH, PAPER_DISTANCE, units: :meters)
    assert_in_delta PAPER_TO.lat, res.endpoint.lat, 1e-8
    assert_in_delta PAPER_TO.lng, res.endpoint.lng, 1e-8
    assert_in_delta PAPER_FINAL_AZIMUTH, res.final_azimuth, 1e-6
  end

  def test_wgs84_inverse_vectors
    WGS84_VECTORS.each do |from, to, distance, initial_azimuth, final_azimuth|
      res = WGS84.inverse(from, to, units: :meters)
      assert_in_delta distance, res.distance, DISTANCE_TOLERANCE_METERS,
                      "distance off for #{from.inspect} -> #{to.inspect}"
      assert_in_delta initial_azimuth, res.initial_azimuth, AZIMUTH_TOLERANCE_DEG,
                      "initial azimuth off for #{from.inspect} -> #{to.inspect}"
      assert_in_delta final_azimuth, res.final_azimuth, AZIMUTH_TOLERANCE_DEG,
                      "final azimuth off for #{from.inspect} -> #{to.inspect}"
    end
  end

  def test_wgs84_direct_vectors
    WGS84_VECTORS.each do |from, to, distance, initial_azimuth, final_azimuth|
      res = WGS84.direct(from, initial_azimuth, distance, units: :meters)
      assert_in_delta to[0], res.endpoint.lat, LATLNG_TOLERANCE_DEG,
                      "endpoint lat off for #{from.inspect}"
      assert_lng_in_delta to[1], res.endpoint.lng,
                          "endpoint lng off for #{from.inspect}"
      assert_in_delta final_azimuth, res.final_azimuth, AZIMUTH_TOLERANCE_DEG,
                      "final azimuth off for #{from.inspect}"
    end
  end

  def test_equatorial_quadrant_is_exact
    res = WGS84.inverse([0, 0], [0, 90], units: :meters)
    assert_in_delta EQUATORIAL_QUADRANT_METERS, res.distance, DISTANCE_TOLERANCE_METERS
    assert_in_delta 90.0, res.initial_azimuth, AZIMUTH_TOLERANCE_DEG
    assert_in_delta 90.0, res.final_azimuth, AZIMUTH_TOLERANCE_DEG
  end

  def test_meridian_quadrant
    res = WGS84.inverse([0, 0], [90, 0], units: :meters)
    assert_in_delta MERIDIAN_QUADRANT_METERS, res.distance, DISTANCE_TOLERANCE_METERS
    assert_in_delta 0.0, res.initial_azimuth, AZIMUTH_TOLERANCE_DEG
  end

  def test_inverse_direct_round_trip
    WGS84_VECTORS.each do |from, to, distance, initial_azimuth, _final_azimuth|
      forward = WGS84.direct(from, initial_azimuth, distance, units: :meters)
      back = WGS84.inverse(from, forward.endpoint, units: :meters)
      assert_in_delta distance, back.distance, DISTANCE_TOLERANCE_METERS
      assert_in_delta initial_azimuth, back.initial_azimuth, AZIMUTH_TOLERANCE_DEG
    end
  end

  def test_coincident_points_return_zero
    res = WGS84.inverse([45, 90], [45, 90], units: :meters)
    assert_equal 0.0, res.distance
    assert_equal 0.0, res.initial_azimuth
    assert_equal 0.0, res.final_azimuth
  end

  def test_coincident_across_antimeridian
    # (0, 180) and (0, -180) are the same physical point.
    res = WGS84.inverse([0, 180], [0, -180], units: :meters)
    assert_equal 0.0, res.distance
  end

  def test_direct_with_zero_distance_returns_start
    res = WGS84.direct([45, 90], 123, 0, units: :meters)
    assert_in_delta 45, res.endpoint.lat, 1e-12
    assert_in_delta 90, res.endpoint.lng, 1e-12
    assert_equal 123.0, res.final_azimuth
  end

  def test_antipodal_points_raise_not_converged
    # Exact antipodes have infinitely many equally short geodesics; the
    # iteration cannot converge and must raise instead of silently
    # returning a spherical approximation.
    assert_raise(Geokit::Ellipsoid::NotConvergedError) do
      WGS84.inverse([0, 0], [0, 180], units: :meters)
    end
    assert_raise(Geokit::Ellipsoid::NotConvergedError) do
      WGS84.inverse([10, 20], [-10, -160], units: :meters)
    end
  end

  def test_antimeridian_longitude_wrapping
    # Shortest path crosses the antimeridian, not the long way around.
    res = WGS84.inverse([10, 179.5], [10, -179.5], units: :meters)
    assert res.distance < 120_000, 'should take the short path across 180'
    direct = WGS84.direct([10, 179.5], res.initial_azimuth, res.distance, units: :meters)
    assert_in_delta(-179.5, direct.endpoint.lng, LATLNG_TOLERANCE_DEG)
    assert_operator direct.endpoint.lng, :>=, -180
    assert_operator direct.endpoint.lng, :<, 180
  end

  def test_result_longitudes_are_normalized
    res = WGS84.direct([0, 179], 90, 500_000, units: :meters)
    assert_operator res.endpoint.lng, :>=, -180
    assert_operator res.endpoint.lng, :<, 180
    res = WGS84.direct([0, -179], 270, 500_000, units: :meters)
    assert_operator res.endpoint.lng, :>=, -180
    assert_operator res.endpoint.lng, :<, 180
  end

  def test_units_are_shared_between_inverse_and_direct
    meters = WGS84.inverse([0, 0], [10, 20], units: :meters).distance
    kms = WGS84.inverse([0, 0], [10, 20], units: :kms).distance
    miles = WGS84.inverse([0, 0], [10, 20], units: :miles).distance
    assert_in_delta meters / 1000.0, kms, 1e-9
    assert_in_delta meters / 1609.0, miles, 1e-9

    # direct consumes the same units inverse produces
    kms_back = WGS84.direct([0, 0], 62.879291486367884, kms, units: :kms)
    assert_in_delta 10, kms_back.endpoint.lat, LATLNG_TOLERANCE_DEG
    assert_in_delta 20, kms_back.endpoint.lng, LATLNG_TOLERANCE_DEG
  end

  def test_default_units_follow_geokit_default_units
    old = Geokit.default_units
    Geokit.default_units = :kms
    begin
      res = WGS84.inverse([0, 0], [10, 20])
      assert_in_delta 2475.7245187498556, res.distance, 1e-6
    ensure
      Geokit.default_units = old
    end
  end

  def test_unsupported_units_raise
    assert_raise(ArgumentError) { WGS84.inverse([0, 0], [1, 1], units: :furlongs) }
    assert_raise(ArgumentError) { WGS84.direct([0, 0], 45, 100, units: :furlongs) }
  end

  def test_ellipsoid_parameter_validation
    assert_equal 6378137.0, WGS84.semi_major_axis
    assert_in_delta 1 / 298.257223563, WGS84.flattening, 1e-18
    assert_raise(ArgumentError) { Geokit::Ellipsoid.new(0, 0.003) }
    assert_raise(ArgumentError) { Geokit::Ellipsoid.new(-6378137, 0.003) }
    assert_raise(ArgumentError) { Geokit::Ellipsoid.new(6378137, -0.001) }
    assert_raise(ArgumentError) { Geokit::Ellipsoid.new(6378137, 1.0) }
    assert_raise(ArgumentError) { Geokit::Ellipsoid.new(6378137, 1.5) }
    # a sphere (f = 0) is a valid degenerate ellipsoid
    sphere = Geokit::Ellipsoid.new(6378137.0, 0.0)
    res = sphere.inverse([0, 0], [0, 90], units: :meters)
    assert_in_delta 6378137.0 * Math::PI / 2, res.distance, DISTANCE_TOLERANCE_METERS
  end

  def test_mappable_distance_between_with_ellipsoid_formula
    from = Geokit::LatLng.new(51.5074, -0.1278)
    to = Geokit::LatLng.new(40.7128, -74.0060)
    assert_in_delta 5585.2335789467786,
                    Geokit::LatLng.distance_between(from, to, units: :kms, formula: :ellipsoid),
                    1e-7
    assert_in_delta 5585.2335789467786,
                    from.distance_to(to, units: :kms, formula: :ellipsoid),
                    1e-7
  end

  def test_mappable_endpoint_with_ellipsoid_formula
    from = Geokit::LatLng.new(51.5074, -0.1278)
    endpoint = from.endpoint(288.36899811332881, 5585.2335789467786,
                             units: :kms, formula: :ellipsoid)
    assert_in_delta 40.7128, endpoint.lat, LATLNG_TOLERANCE_DEG
    assert_in_delta(-74.0060, endpoint.lng, LATLNG_TOLERANCE_DEG)
  end

  def test_default_formula_remains_sphere
    from = Geokit::LatLng.new(51.5074, -0.1278)
    to = Geokit::LatLng.new(40.7128, -74.0060)
    default = Geokit::LatLng.distance_between(from, to, units: :kms)
    sphere = Geokit::LatLng.distance_between(from, to, units: :kms, formula: :sphere)
    ellipsoid = Geokit::LatLng.distance_between(from, to, units: :kms, formula: :ellipsoid)
    assert_equal sphere, default
    # sphere and ellipsoid differ by ~0.3% here; they must not be conflated
    assert (sphere - ellipsoid).abs > 1.0
  end

  # Compares longitudes modulo 360: -180 and 180 are the same meridian.
  def assert_lng_in_delta(expected, actual, message = nil)
    wrapped = ((expected - actual + 180) % 360) - 180
    assert_in_delta 0, wrapped, LATLNG_TOLERANCE_DEG, message
  end

  def test_default_ellipsoid_is_wgs84_and_configurable
    assert_same Geokit::Ellipsoid::WGS84, Geokit.default_ellipsoid
    old = Geokit.default_ellipsoid
    begin
      Geokit.default_ellipsoid = BESSEL
      res = Geokit::LatLng.distance_between(PAPER_FROM, PAPER_TO,
                                            units: :meters, formula: :ellipsoid)
      assert_in_delta PAPER_DISTANCE, res, 0.001
    ensure
      Geokit.default_ellipsoid = old
    end
  end
end
