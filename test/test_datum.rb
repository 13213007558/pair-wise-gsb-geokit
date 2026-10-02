require File.join(File.dirname(__FILE__), 'helper')

class DatumEllipsoidTest < Test::Unit::TestCase
  def test_ellipsoid_geometry
    wgs = Geokit::Ellipsoid::WGS84
    assert_in_delta 1.0 / 298.257_223_563, wgs.flattening, 1.0e-15
    assert_in_delta 6_356_752.314_245, wgs.semi_minor_axis, 1.0e-6
    assert_in_delta 6.694_379_990_14e-3,
                    wgs.first_eccentricity_squared, 1.0e-13
  end

  def test_ellipsoid_rejects_bad_parameters
    assert_raise(ArgumentError) do
      Geokit::Ellipsoid.new(semi_major_axis: -1.0,
                           inverse_flattening: 298.0)
    end
    assert_raise(ArgumentError) do
      Geokit::Ellipsoid.new(semi_major_axis: '6378137',
                           inverse_flattening: 298.0)
    end
  end

  def test_sphere_has_zero_flattening
    sphere = Geokit::Ellipsoid.new(semi_major_axis: 1.0)
    assert sphere.sphere?
    assert_equal 0.0, sphere.flattening
    assert_equal 0.0, sphere.first_eccentricity_squared
  end
end

class GeodeticConversionTest < Test::Unit::TestCase
  D = Geokit::Datum::WGS84

  # Fixed control values computed independently from the closed-form
  # geodetic -> ECEF formulas on the WGS84 ellipsoid.
  FIXED_ECEF = {
    [45.0, 45.0, 0.0] => [3_194_419.145061, 3_194_419.145061,
                           4_487_348.408866],
    [45.0, 45.0, -1500.0] => [3_193_669.145061, 3_193_669.145061,
                               4_486_287.748694]
  }.freeze

  def test_geodetic_to_ecef_fixed_control_values
    FIXED_ECEF.each do |(lat, lng, height), xyz|
      point = D.to_ecef(Geokit::GeodeticPoint.new(lat, lng, height))
      assert_array_in_delta(xyz, point.to_a, 1.0e-6)
    end
  end

  def test_geodetic_to_ecef_accepts_lat_lng_height
    point = D.to_ecef(45.0, 45.0, 0.0)
    assert_in_delta 3_194_419.145061, point.x, 1.0e-6
    assert_in_delta 3_194_419.145061, point.y, 1.0e-6
    assert_in_delta 4_487_348.408866, point.z, 1.0e-6
  end

  def test_negative_height_round_trip
    point = Geokit::GeodeticPoint.new(45.0, 45.0, -1500.0)
    back = D.from_ecef(D.to_ecef(point))
    assert_in_delta point.lat, back.lat, 1.0e-10
    assert_in_delta point.lng, back.lng, 1.0e-10
    assert_in_delta point.height, back.height, 1.0e-7
  end

  def test_near_north_pole_round_trip
    point = Geokit::GeodeticPoint.new(89.9999, 123.0, -500.0)
    back = D.from_ecef(D.to_ecef(point))
    assert_in_delta point.lat, back.lat, 1.0e-10
    assert_in_delta point.lng, back.lng, 1.0e-10
    assert_in_delta point.height, back.height, 1.0e-7
  end

  def test_near_south_pole_round_trip
    point = Geokit::GeodeticPoint.new(-89.9999, -45.0, 10.0)
    back = D.from_ecef(D.to_ecef(point))
    assert_in_delta point.lat, back.lat, 1.0e-10
    assert_in_delta point.lng, back.lng, 1.0e-10
    assert_in_delta point.height, back.height, 1.0e-7
  end

  def test_exact_poles
    a = D.ellipsoid.semi_major_axis
    b = D.ellipsoid.semi_minor_axis
    polar_factor = 1.0 - D.ellipsoid.flattening

    north_height = 25.0
    north = D.from_ecef(Geokit::EcefPoint.new(0.0, 0.0,
                                              b + polar_factor * north_height))
    assert_in_delta 90.0, north.lat, 1.0e-12
    assert_in_delta north_height, north.height, 1.0e-7

    south_height = -40.0
    south = D.from_ecef(Geokit::EcefPoint.new(0.0, 0.0,
                                              -b + polar_factor * south_height.abs))
    assert_in_delta(-90.0, south.lat, 1.0e-12)
    assert_in_delta(south_height, south.height, 1.0e-7)
  end

  def test_sphere_round_trip
    sphere = Geokit::Datum.new('Sphere',
                              Geokit::Ellipsoid.new(semi_major_axis: 6_371_000.0))
    point = Geokit::GeodeticPoint.new(30.0, 60.0, -1000.0, datum: sphere)
    ecef = sphere.to_ecef(point)
    assert_array_in_delta([2_758_290.911053, 4_777_500.0, 3_185_000.0],
                          ecef.to_a, 1.0e-6)
    back = sphere.from_ecef(ecef)
    assert_in_delta 30.0, back.lat, 1.0e-10
    assert_in_delta 60.0, back.lng, 1.0e-10
    assert_in_delta(-1000.0, back.height, 1.0e-7)
  end

  def test_lat_lng_has_no_datum_and_never_converts_implicitly
    point = Geokit::LatLng.new(51.5, -1.0)
    assert !point.respond_to?(:datum)
    geodetic = Geokit::GeodeticPoint.from_lat_lng(point, D, height: 50.0)
    assert_same D, geodetic.datum
    assert !geodetic.to_lat_lng.respond_to?(:datum)
    assert_equal '51.5,-1.0', geodetic.to_lat_lng.ll
  end

  def test_geodetic_point_without_datum_refuses_ecef_conversion
    point = Geokit::GeodeticPoint.new(51.5, -1.0, 50.0)
    assert_nil point.datum
    assert_raise(ArgumentError) { point.to_ecef }
  end

  def test_transform_requires_matching_source_datum
    point = Geokit::GeodeticPoint.new(51.5, -1.0, 50.0,
                                     datum: Geokit::Datum::ED50)
    error = assert_raise(ArgumentError) do
      Geokit::Datum::OSGB36.transform(point,
                                     Geokit::Helmert::OSGB36_TO_WGS84,
                                     Geokit::Datum::WGS84)
    end
    assert_match(/ED50/, error.message)
    assert_match(/OSGB36/, error.message)
  end

  def test_geodetic_point_validates_latitude
    assert_raise(ArgumentError) do
      Geokit::GeodeticPoint.new(90.1, 0.0)
    end
  end

  def test_osgb36_to_wgs84_round_trip
    osgb = Geokit::Datum::OSGB36
    helmert = Geokit::Helmert::OSGB36_TO_WGS84
    point = Geokit::GeodeticPoint.new(51.5, -1.0, 50.0)

    wgs = osgb.transform(point, helmert, Geokit::Datum::WGS84)
    # Independently computed EPSG:1314 forward control values
    # (position-vector convention).
    wgs_ecef = Geokit::Datum::WGS84.to_ecef(wgs)
    assert_array_in_delta([3_978_220.949567, -69_544.738451,
                           4_968_676.621115], wgs_ecef.to_a, 1.0e-3)
    assert_in_delta 51.500489947, wgs.lat, 1.0e-7
    assert_in_delta(-1.001506510, wgs.lng, 1.0e-7)
    assert_in_delta 358.070064439, wgs.height, 1.0e-6

    back = Geokit::Datum::WGS84.transform(wgs, helmert.inverse, osgb)
    assert_in_delta point.lat, back.lat, 1.0e-9
    assert_in_delta point.lng, back.lng, 1.0e-9
    assert_in_delta point.height, back.height, 1.0e-6
  end
end
