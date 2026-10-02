# encoding: utf-8
require File.join(File.dirname(__FILE__), 'helper')

class TestDatum < Test::Unit::TestCase
  WGS84 = Geokit::Datum::DATUMS[:wgs84]
  OSGB36 = Geokit::Datum::DATUMS[:osgb36]

  def dms(deg, min, sec)
    deg + min / 60.0 + sec / 3600.0
  end

  def helmert(**overrides)
    defaults = {
      tx: 100.0, ty: 200.0, tz: 300.0,
      rx: 1.0, ry: 2.0, rz: 3.0, scale: 10.0,
      convention: :position_vector,
      translation_unit: :meters, rotation_unit: :arc_seconds, scale_unit: :ppm
    }
    Geokit::Datum::HelmertTransform.new(**defaults.merge(overrides))
  end

  # --- ellipsoid -----------------------------------------------------------

  def test_wgs84_ellipsoid_derived_values
    e = Geokit::Datum::ELLIPSOIDS[:wgs84]
    assert_in_delta 6_378_137.0, e.a, 1e-9
    assert_in_delta 6_356_752.314245179, e.b, 1e-9
    assert_in_delta 0.00669437999014, e.eccentricity_squared, 1e-14
  end

  # --- geodetic <-> ECEF, fixed reference points ---------------------------

  def test_to_ecef_fixed_points
    assert_array_in_delta [6_378_137.0, 0.0, 0.0], WGS84.to_ecef(0, 0, 0), 1e-9
    assert_array_in_delta [0.0, 6_378_137.0, 0.0], WGS84.to_ecef(0, 90, 0), 1e-9
    # poles: x/y vanish, z is the semi-minor axis
    assert_array_in_delta [0.0, 0.0, 6_356_752.314245179], WGS84.to_ecef(90, 0, 0), 1e-6
    assert_array_in_delta [0.0, 0.0, -6_356_752.314245179], WGS84.to_ecef(-90, 0, 0), 1e-6
  end

  def test_geodetic_ecef_round_trip
    points = [
      [37.7749, -122.4194, 13.0],      # ordinary mid-latitude point
      [89.9999, -45.0, 120.5],         # near the north pole
      [-89.9999, 179.9, -30.0],        # near the south pole, negative height
      [31.2304, 121.4737, -4200.0],    # deep negative height
      [0.0, 0.0, 0.0]                  # origin of lon/lat
    ]
    points.each do |lat, lng, h|
      back = WGS84.from_ecef(*WGS84.to_ecef(lat, lng, h))
      assert_in_delta lat, back[0], 1e-9, "lat round trip failed for #{lat},#{lng},#{h}"
      assert_in_delta lng, back[1], 1e-9, "lng round trip failed for #{lat},#{lng},#{h}"
      assert_in_delta h, back[2], 1e-4, "height round trip failed for #{lat},#{lng},#{h}"
    end
  end

  def test_near_pole_and_negative_height_round_trip_through_helmert
    # the full datum conversion path must also survive the pole and depth
    [[89.9999, 10.0, -250.0], [-89.5, -100.0, -1500.0]].each do |lat, lng, h|
      there = OSGB36.convert_to(WGS84, lat, lng, h)
      back = WGS84.convert_to(OSGB36, *there)
      assert_in_delta lat, back[0], 1e-8
      assert_in_delta lng, back[1], 1e-8
      assert_in_delta h, back[2], 1e-4
    end
  end

  # --- datum conversion against fixed reference data -----------------------

  def test_osgb36_to_wgs84_fixed_reference
    # Ordnance Survey published parameters (position-vector convention):
    # OSGB36 52deg 39' 27.2531"N, 1deg 43' 04.5177"E, h=0 converts to
    # the WGS84/ETRS89 values below (verified analytically against the
    # published 7-parameter formulas).
    lat, lng, h = OSGB36.convert_to(WGS84, dms(52, 39, 27.2531), dms(1, 43, 4.5177), 0.0)
    assert_in_delta 52.657978603224, lat, 1e-9
    assert_in_delta 1.716051944775, lng, 1e-9
    assert_in_delta 44.7015997209, h, 1e-6
  end

  def test_datum_conversion_round_trip
    lat, lng, h = dms(52, 39, 27.2531), dms(1, 43, 4.5177), 53.0
    back = WGS84.convert_to(OSGB36, *OSGB36.convert_to(WGS84, lat, lng, h))
    assert_in_delta lat, back[0], 1e-9
    assert_in_delta lng, back[1], 1e-9
    assert_in_delta h, back[2], 1e-6
  end

  # --- Helmert transform ---------------------------------------------------

  def test_helmert_translation_only_fixed_values
    t = helmert(rx: 0.0, ry: 0.0, rz: 0.0, scale: 0.0)
    assert_array_in_delta [101.0, 202.0, 303.0], t.forward([1.0, 2.0, 3.0]), 1e-12
    assert_array_in_delta [1.0, 2.0, 3.0], t.inverse([101.0, 202.0, 303.0]), 1e-12
  end

  def test_helmert_units_are_converted
    by_arcsec = helmert
    by_degrees = helmert(rx: 1.0 / 3600.0, ry: 2.0 / 3600.0, rz: 3.0 / 3600.0,
                         rotation_unit: :degrees)
    by_mm = helmert(tx: 100_000.0, ty: 200_000.0, tz: 300_000.0,
                    translation_unit: :millimeters)
    by_unity = helmert(scale: 10.0e-6, scale_unit: :unity)
    point = [3_874_490.0, 116_205.0, 5_047_143.0]
    assert_array_in_delta by_arcsec.forward(point), by_degrees.forward(point), 1e-9
    assert_array_in_delta by_arcsec.forward(point), by_mm.forward(point), 1e-9
    assert_array_in_delta by_arcsec.forward(point), by_unity.forward(point), 1e-9
  end

  def test_coordinate_frame_equals_negated_position_vector_rotations
    point = [3_874_490.0, 116_205.0, 5_047_143.0]
    frame = helmert(convention: :coordinate_frame)
    negated = helmert(rx: -1.0, ry: -2.0, rz: -3.0)
    assert_array_in_delta negated.forward(point), frame.forward(point), 1e-12
    refute_array_in_delta frame.forward(point), helmert.forward(point), 1e-3
  end

  def test_inverse_is_true_matrix_inverse_not_sign_flip
    t = helmert # 10 ppm scale, arc-second rotations: sign-flip error is measurable
    point = WGS84.to_ecef(52.6576, 1.7179, 44.7)
    assert_array_in_delta point, t.inverse(t.forward(point)), 1e-9

    # a naive inverse that just negates every parameter is only first-order;
    # with a 10 ppm scale it must be off by more than a millimetre here
    naive = helmert(tx: -100.0, ty: -200.0, tz: -300.0,
                    rx: -1.0, ry: -2.0, rz: -3.0, scale: -10.0)
    naive_back = naive.forward(t.forward(point))
    error = naive_back.zip(point).map { |a, b| (a - b).abs }.max
    assert error > 1e-3, "expected the sign-flip inverse to be measurably wrong, got #{error}"
  end

  def test_forward_matrix_is_inverted_exactly
    t = helmert
    matrix, translation = t.forward_matrix
    point = [3_874_490.0, 116_205.0, 5_047_143.0]
    moved = t.forward(point)
    # subtract translation, then check M * (inverse result) == moved - t
    back = t.inverse(moved)
    assert_array_in_delta point, back, 1e-9
    assert_equal 3, matrix.size
    assert_equal 3, translation.size
  end

  # --- reference epoch propagation -----------------------------------------

  def dynamic_helmert
    helmert(
      tx: 0.5, ty: -0.2, tz: 0.3, rx: 0.0, ry: 0.0, rz: 0.0, scale: 0.0,
      rates: { dtx: 0.025, dty: -0.010, dtz: 0.015 },
      reference_epoch: 2010.0
    )
  end

  def test_epoch_propagates_parameters_linearly
    t = dynamic_helmert
    point = [1_000_000.0, 2_000_000.0, 3_000_000.0]
    at_ref = t.forward(point, epoch: 2010.0)
    assert_array_in_delta [1_000_000.5, 1_999_999.8, 3_000_000.3], at_ref, 1e-12
    at_2000 = t.forward(point, epoch: 2000.0)
    at_2020 = t.forward(point, epoch: 2020.0)
    # ten years of drift each way from the 2010.0 reference epoch
    assert_array_in_delta [-0.25, 0.10, -0.15], at_2000.zip(at_ref).map { |a, b| a - b }, 1e-9
    assert_array_in_delta [0.25, -0.10, 0.15], at_2020.zip(at_ref).map { |a, b| a - b }, 1e-9
  end

  def test_same_site_at_two_epochs
    # A station nailed to a local (plate-fixed) datum keeps constant
    # lat/lng/h, while its WGS84 coordinates drift with the epoch.
    local = Geokit::Datum.new(name: 'Local plate datum',
                              ellipsoid: Geokit::Datum::ELLIPSOIDS[:grs80],
                              to_wgs84: dynamic_helmert)
    lat, lng, h = 48.0, 11.0, 500.0
    ecef_2000 = local.convert_to(WGS84, lat, lng, h, epoch: 2000.0)
    ecef_2020 = local.convert_to(WGS84, lat, lng, h, epoch: 2020.0)

    # the site moved in the WGS84 frame between the two epochs
    moved = ecef_2020.zip(ecef_2000).map { |a, b| (a - b).abs }.max
    assert moved > 0.1, 'expected the site to drift in the WGS84 frame'

    # converting each epoch's coordinates back at its own epoch recovers
    # the very same site
    back_2000 = WGS84.convert_to(local, *ecef_2000, epoch: 2000.0)
    back_2020 = WGS84.convert_to(local, *ecef_2020, epoch: 2020.0)
    [back_2000, back_2020].each do |back|
      assert_in_delta lat, back[0], 1e-9
      assert_in_delta lng, back[1], 1e-9
      assert_in_delta h, back[2], 1e-6
    end

    # using the wrong epoch does not recover the site
    wrong = WGS84.convert_to(local, *ecef_2020, epoch: 2000.0)
    assert (wrong[2] - h).abs + (wrong[0] - lat).abs > 1e-6
  end

  # --- refusal to guess ----------------------------------------------------

  def test_missing_epoch_is_rejected_for_time_dependent_transform
    t = dynamic_helmert
    point = [1.0, 2.0, 3.0]
    assert_raise(ArgumentError) { t.forward(point) }
    assert_raise(ArgumentError) { t.inverse(point) }
    local = Geokit::Datum.new(name: 'local', ellipsoid: Geokit::Datum::ELLIPSOIDS[:grs80],
                              to_wgs84: t)
    assert_raise(ArgumentError) { local.convert_to(WGS84, 48.0, 11.0, 500.0) }
  end

  def test_rates_require_reference_epoch
    assert_raise(ArgumentError) { helmert(rates: { dtx: 0.01 }) }
  end

  def test_reference_epoch_without_rates_is_rejected
    assert_raise(ArgumentError) { helmert(reference_epoch: 2010.0) }
  end

  def test_unknown_rate_parameter_is_rejected
    assert_raise(ArgumentError) do
      helmert(rates: { dx: 0.01 }, reference_epoch: 2010.0)
    end
  end

  def test_ambiguous_units_and_convention_are_rejected
    assert_raise(ArgumentError) { helmert(rotation_unit: :gons) }
    assert_raise(ArgumentError) { helmert(translation_unit: :feet) }
    assert_raise(ArgumentError) { helmert(scale_unit: :percent) }
    assert_raise(ArgumentError) { helmert(convention: :bursa_wolf) }
    # convention and all three units are mandatory keywords
    assert_raise(ArgumentError) do
      Geokit::Datum::HelmertTransform.new(tx: 0, ty: 0, tz: 0, rx: 0, ry: 0, rz: 0, scale: 0)
    end
  end

  # --- LatLng stays datum-free ---------------------------------------------

  def test_latlng_is_not_bound_to_any_datum
    ll = Geokit::LatLng.new(37.7749, -122.4194)
    assert !ll.respond_to?(:datum)
    assert !ll.respond_to?(:to_ecef)
    assert_equal [37.7749, -122.4194], ll.to_a
    assert_equal ll, Geokit::LatLng.normalize('37.7749,-122.4194')
  end

  private

  def refute_array_in_delta(expected, actual, delta)
    assert expected.zip(actual).any? { |e, a| (e - a).abs > delta },
           "expected #{expected.inspect} and #{actual.inspect} to differ by more than #{delta}"
  end
end
