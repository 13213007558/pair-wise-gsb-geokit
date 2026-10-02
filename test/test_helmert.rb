require File.join(File.dirname(__FILE__), 'helper')

class HelmertTest < Test::Unit::TestCase
  def base_helmert(**options)
    Geokit::Helmert.new(
      **{
        convention: :position_vector,
        translation: [100.0, -200.0, 300.0],
        rotation: [1.0, -2.0, 3.0],
        scale: 2.0,
        translation_unit: :metres,
        rotation_unit: :arcseconds,
        scale_unit: :ppm
      }.merge(options)
    )
  end

  def test_forward_matrix_position_vector
    helmert = base_helmert
    m = helmert.matrix_at(nil)
    assert_in_delta 1.000002, m[0][0], 1.0e-15
    assert_in_delta(-3.0 / 3600 * Math::PI / 180 * 1.000002,
                    m[0][1], 1.0e-15)
    assert_in_delta(-1.0 / 3600 * Math::PI / 180 * 1.000002,
                    m[1][2], 1.0e-15)
    assert_in_delta 300.0, m[2][3], 1.0e-12
  end

  def test_coordinate_frame_reverses_rotation_signs
    position = base_helmert.matrix_at(nil)
    frame = base_helmert(convention: :coordinate_frame).matrix_at(nil)
    [0, 1, 2].product([0, 1, 2]) do |row, col|
      next if row == col
      assert_in_delta(-position[row][col], frame[row][col], 1.0e-18)
    end
    assert_in_delta position[0][3], frame[0][3], 1.0e-18
  end

  def test_convention_is_required
    error = assert_raise(ArgumentError) do
      base_helmert(convention: nil)
    end
    assert_match(/position_vector/, error.message)
  end

  def test_unknown_convention_rejected
    assert_raise(ArgumentError) do
      base_helmert(convention: :bursa_wolf)
    end
  end

  def test_translation_unit_required_and_applied
    error = assert_raise(ArgumentError) do
      base_helmert(translation_unit: nil)
    end
    assert_match(/translation_unit/, error.message)

    helmert = base_helmert(translation_unit: :millimetres,
                          translation: [1000.0, 0.0, 0.0])
    assert_in_delta 1.0, helmert.matrix_at(nil)[0][3], 1.0e-15
  end

  def test_rotation_unit_required_and_applied
    assert_raise(ArgumentError) { base_helmert(rotation_unit: nil) }

    radians = 1.0e-6
    helmert = base_helmert(rotation_unit: :radians,
                          rotation: [radians, 0.0, 0.0])
    assert_in_delta(-radians * 1.000002, helmert.matrix_at(nil)[1][2], 1.0e-12)

    degree = base_helmert(rotation_unit: :degrees, rotation: [1.0, 0.0, 0.0],
                          scale: 0.0)
    assert_in_delta(-Math::PI / 180, degree.matrix_at(nil)[1][2], 1.0e-15)
  end

  def test_scale_unit_required_and_applied
    assert_raise(ArgumentError) { base_helmert(scale_unit: nil) }

    ratio = base_helmert(scale_unit: :ratio, scale: 0.0001)
    assert_in_delta 1.0001, ratio.matrix_at(nil)[0][0], 1.0e-15

    ppb = base_helmert(scale_unit: :ppb, scale: 1000.0)
    assert_in_delta 1.000001, ppb.matrix_at(nil)[0][0], 1.0e-15
  end

  def test_feet_and_milliarcsecond_units
    # 3.280839895... feet is exactly one metre; 1 mas is 1/1000 arcsec.
    helmert = Geokit::Helmert.new(
      convention: :position_vector,
      translation: [3.2808398950131235, 0.0, 0.0],
      translation_unit: :feet,
      rotation: [1000.0, 0.0, 0.0],
      rotation_unit: :milliarcseconds,
      scale: 1000.0,
      scale_unit: :ppb
    )
    matrix = helmert.matrix_at(nil)
    assert_in_delta 1.0, matrix[0][3], 1.0e-7
    assert_in_delta(-1.0 / 3600 * Math::PI / 180,
                    matrix[1][2], 1.0e-9)
    assert_in_delta 1.000001, matrix[0][0], 1.0e-12
  end

  def test_forward_application
    helmert = base_helmert(translation: [10.0, 20.0, 30.0],
                           rotation: [0.0, 0.0, 0.0],
                           scale: 0.0)
    point = Geokit::EcefPoint.new(1000.0, 2000.0, 3000.0)
    transformed = helmert.forward(point)
    assert_array_in_delta([1010.0, 2020.0, 3030.0], transformed.to_a,
                          1.0e-9)
  end

  def test_backward_is_numerical_inverse_with_nonzero_scale
    helmert = base_helmert
    point = Geokit::EcefPoint.new(3_977_000.0, -70_000.0,
                                 4_968_000.0)
    forward = helmert.forward(point)
    back = helmert.backward(forward)
    assert_array_in_delta(point.to_a, back.to_a, 1.0e-6)

    identity_error = product_error(helmert)
    assert identity_error < 1.0e-12
  end

  def test_sign_flipping_is_not_the_inverse_when_scale_is_nonzero
    # x' = (1 + 500 ppm) x. Negating every parameter gives x'' =
    # (1 - 500 ppm)(1 + 500 ppm) x, which leaves a -0.25 ppm x residual;
    # the genuine inverse uses the inverted matrix (1 / (1 + 500 ppm)).
    helmert = Geokit::Helmert.new(
      convention: :position_vector,
      translation: [0.0, 0.0, 0.0],
      rotation: [0.0, 0.0, 0.0],
      scale: 500.0,
      translation_unit: :metres,
      rotation_unit: :radians,
      scale_unit: :ppm
    )
    point = Geokit::EcefPoint.new(1_000_000.0, 0.0, 0.0)
    sign_flipped = Geokit::Helmert.new(
      convention: :position_vector,
      translation: [0.0, 0.0, 0.0],
      rotation: [0.0, 0.0, 0.0],
      scale: -500.0,
      translation_unit: :metres,
      rotation_unit: :radians,
      scale_unit: :ppm
    )
    wrong_back = sign_flipped.forward(helmert.forward(point))
    residual = wrong_back.to_a.zip(point.to_a).map { |a, b| (a - b).abs }.max
    assert_in_delta 0.25, residual, 1.0e-9

    correct_back = helmert.backward(helmert.forward(point))
    assert_in_delta point.x, correct_back.x, 1.0e-6
  end

  def test_inverse_wrapper_round_trips
    helmert = base_helmert
    point = Geokit::EcefPoint.new(1_000_000.0, -2_000_000.0,
                                 5_000_000.0)
    wrapped = helmert.inverse
    round_trip = wrapped.forward(helmert.forward(point))
    assert_array_in_delta(point.to_a, round_trip.to_a, 1.0e-6)
    assert_same helmert, wrapped.inverse
  end

  def test_non_ecef_point_rejected
    assert_raise(ArgumentError) do
      base_helmert.forward(Geokit::GeodeticPoint.new(0, 0, 0))
    end
  end

  private

  def product_error(helmert)
    inverse = Geokit::Helmert.invert_matrix(helmert.matrix_at(nil))
    forward = helmert.matrix_at(nil)
    [0, 1, 2, 3].map do |row|
      [0, 1, 2, 3].map do |col|
        value = [0, 1, 2, 3].map { |k| inverse[row][k] * forward[k][col] }.sum
        value - (row == col ? 1.0 : 0.0)
      end
    end.flatten.map(&:abs).max
  end
end

class HelmertRatesTest < Test::Unit::TestCase
  def rates_helmert(**options)
    Geokit::Helmert.new(
      **{
        convention: :position_vector,
        translation: [-1.0, 2.0, 3.0],
        rotation: [0.02, -0.03, 0.05],
        scale: -0.5,
        translation_unit: :metres,
        rotation_unit: :milliarcseconds,
        scale_unit: :ppb,
        rates: {
          translation: [21.0, -8.0, 5.0],
          rotation: [0.1, -0.2, 0.15],
          scale: 0.02
        },
        rate_translation_unit: :millimetres_per_year,
        rate_rotation_unit: :milliarcseconds_per_year,
        rate_scale_unit: :ppb_per_year,
        reference_epoch: 2010.0
      }.merge(options)
    )
  end

  def wgs84_origin
    Geokit::Datum::WGS84.to_ecef(Geokit::GeodeticPoint.new(0.0, 0.0, 0.0))
  end

  # Independent control values: 6378137 m on the equator transformed with
  # the linear rates at 2005.0 and 2015.0.
  EPOCH_2005 = [6_378_135.891173117, 2.018354543,
                 2.945005582].freeze
  EPOCH_2015 = [6_378_136.102448746, 1.984737665,
                 3.056849743].freeze

  def test_values_at_two_epochs
    helmert = rates_helmert
    assert helmert.time_dependent?
    assert_array_in_delta(EPOCH_2005, helmert.forward(wgs84_origin,
                                                       epoch: 2005.0).to_a,
                          1.0e-6)
    assert_array_in_delta(EPOCH_2015, helmert.forward(wgs84_origin,
                                                       epoch: 2015.0).to_a,
                          1.0e-6)
  end

  def test_same_physical_point_evolves_between_epochs
    helmert = rates_helmert
    at_2005 = helmert.forward(wgs84_origin, epoch: 2005.0)
    at_2015 = helmert.forward(wgs84_origin, epoch: 2015.0)
    delta = at_2015.to_a.zip(at_2005.to_a).map { |later, earlier| later - earlier }
    assert_array_in_delta([0.211275629, -0.033616879, 0.111844162],
                          delta, 1.0e-6)
  end

  def test_round_trip_at_each_epoch
    helmert = rates_helmert
    [2005.0, 2015.0].each do |epoch|
      forward = helmert.forward(wgs84_origin, epoch: epoch)
      back = helmert.backward(forward, epoch: epoch)
      assert_array_in_delta(wgs84_origin.to_a, back.to_a, 1.0e-6)
    end
  end

  def test_same_physical_place_converted_at_two_epochs
    # One geodetic place (a benchmark on the source datum). Its source
    # ECEF coordinates are fixed; the rate-carrying Helmert yields a
    # different target ECEF position at each observation epoch.
    wgs = Geokit::Datum::WGS84
    place = Geokit::GeodeticPoint.new(0.0, 0.0, 0.0, datum: wgs)
    source_ecef = wgs.to_ecef(place)
    helmert = rates_helmert

    target_2005 = helmert.forward(source_ecef, epoch: 2005.0)
    target_2015 = helmert.forward(source_ecef, epoch: 2015.0)
    drift = target_2015.to_a.zip(target_2005.to_a).map { |later, earlier|
      later - earlier
    }
    assert_array_in_delta([0.211275629, -0.033616879, 0.111844162],
                          drift, 1.0e-6)

    # Both target positions refer back to the same source place, but only
    # when each is inverted at its own epoch.
    [2005.0, 2015.0].each do |epoch|
      target = helmert.forward(source_ecef, epoch: epoch)
      back_source = helmert.backward(target, epoch: epoch)
      assert_array_in_delta(source_ecef.to_a, back_source.to_a, 1.0e-6)
      back_place = wgs.from_ecef(back_source)
      assert_in_delta 0.0, back_place.lat, 1.0e-9
      assert_in_delta 0.0, back_place.lng, 1.0e-9
      assert_in_delta 0.0, back_place.height, 1.0e-6
    end

    # A target position at 2015 must not be inverted using 2005
    # parameters; the 10 years of rates are then left uncompensated.
    wrong_epoch = helmert.backward(target_2015, epoch: 2005.0)
    mismatch = wrong_epoch.to_a.zip(source_ecef.to_a).map { |a, b|
      (a - b).abs
    }.max
    assert mismatch > 0.01
  end

  def test_rate_units_are_explicit
    error = assert_raise(ArgumentError) do
      rates_helmert(rate_translation_unit: nil)
    end
    assert_match(/rate_translation_unit/, error.message)
    assert_raise(ArgumentError) { rates_helmert(rate_rotation_unit: nil) }
    assert_raise(ArgumentError) { rates_helmert(rate_scale_unit: nil) }
  end

  def test_reference_epoch_required_with_rates
    error = assert_raise(ArgumentError) do
      rates_helmert(reference_epoch: nil)
    end
    assert_match(/reference_epoch/, error.message)
  end

  def test_epoch_required_when_applying_rates
    helmert = rates_helmert
    assert_raise(ArgumentError) { helmert.forward(wgs84_origin) }
    assert_raise(ArgumentError) { helmert.backward(wgs84_origin) }
  end

  def test_incomplete_rates_rejected
    assert_raise(ArgumentError) do
      rates_helmert(rates: { translation: [0, 0, 0] })
    end
  end

  def test_static_transform_needs_no_epoch
    helmert = Geokit::Helmert.new(
      convention: :coordinate_frame,
      translation: [1.0, 2.0, 3.0],
      rotation: [0.1, 0.2, 0.3],
      scale: 1.0,
      translation_unit: :metres,
      rotation_unit: :degrees,
      scale_unit: :ppm
    )
    assert !helmert.time_dependent?
    transformed = helmert.forward(wgs84_origin)
    assert transformed.x > 0.0
  end

  def test_reference_epoch_without_rates_is_ambiguous
    assert_raise(ArgumentError) do
      Geokit::Helmert.new(
        convention: :position_vector,
        translation: [0.0, 0.0, 0.0],
        translation_unit: :metres,
        reference_epoch: 2010.0
      )
    end
  end

end
