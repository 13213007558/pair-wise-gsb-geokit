require File.join(File.dirname(__FILE__), 'helper')

class SphericalCapTest < Test::Unit::TestCase #:nodoc: all
  def point(lat, lng, id = nil)
    location = Geokit::LatLng.new(lat, lng)
    location.instance_variable_set(:@identity, id) if id
    location.define_singleton_method(:identity) { @identity }
    location
  end

  def cap(points, options = {})
    Geokit::LatLng.minimum_spherical_cap(points, options)
  end

  def assert_contains(cap, points)
    points.each do |location|
      distance = cap.center.distance_to(
        location, units: cap.units, formula: :sphere
      )
      assert distance <= cap.radius + 1e-6
    end
  end

  def assert_boundary(cap, points, support_ids)
    support = points.values_at(*support_ids)
    assert_equal support_ids, cap.support_ids
    support.each do |location|
      distance = cap.center.distance_to(
        location, units: cap.units, formula: :sphere
      )
      assert_in_delta cap.radius, distance, 1e-6
    end
  end

  def test_single_point_is_zero_radius_and_unique
    only = point(31.2, 121.4, :only)
    result = cap([only])

    assert_equal only, result.center
    assert_in_delta 0.0, result.angular_radius, 1e-12
    assert_equal [:only], result.support_ids
    assert_equal [only], result.support_points
    assert result.unique?
    assert_contains result, [only]
  end

  def test_repeated_points_report_every_boundary_identity
    first = point(10.0, 20.0, :first)
    second = point(10.0, 20.0, :second)
    inside = point(10.5, 20.0, :inside)
    result = cap([first, second, inside], ids: %i[first second inside])

    assert_in_delta 0.25, result.radius_degrees, 1e-9
    assert_equal %i[first second inside], result.support_ids
    assert_contains result, [first, second, inside]
  end

  def test_uses_point_identity_when_explicit_ids_are_omitted
    first = point(10.0, 20.0, :first)
    second = point(10.0, 20.0, :second)
    inside = point(10.5, 20.0, :inside)
    result = cap([first, second, inside])

    assert_equal %i[first second inside], result.support_ids
  end

  def test_two_points_use_their_spherical_midpoint
    points = [point(0.0, 0.0), point(0.0, 90.0)]
    result = cap(points, units: :kms)

    assert_equal 0.0, result.center.lat
    assert_in_delta 45.0, result.center.lng, 1e-10
    assert_in_delta Math::PI / 4.0, result.angular_radius, 1e-12
    assert_in_delta(
      result.angular_radius * Geokit::LatLng.units_sphere_multiplier(:kms),
      result.radius,
      1e-9
    )
    assert_equal [0, 1], result.support_ids
    assert result.unique?
    assert_contains result, points
    assert_boundary result, points, [0, 1]
  end

  def test_points_on_same_great_circle_inside_semicircle
    points = [point(0.0, 0.0), point(0.0, 30.0), point(0.0, 60.0)]
    result = cap(points)

    assert_in_delta 30.0, result.radius_degrees, 1e-9
    assert_equal [0, 2], result.support_ids
    assert result.unique?
    assert_contains result, points
    assert_boundary result, points, [0, 2]
  end

  def test_antipodal_pair_has_nonunique_hemisphere
    points = [point(0.0, 0.0, :a), point(0.0, 180.0, :b)]
    result = cap(points, ids: %i[a b], units: :kms)

    assert_in_delta Math::PI / 2.0, result.angular_radius, 1e-12
    assert_in_delta(
      Geokit::LatLng.units_sphere_multiplier(:kms) * Math::PI / 2.0,
      result.radius,
      1e-9
    )
    assert_equal %i[a b], result.support_ids
    assert !result.unique?
    assert_contains result, points
  end

  def test_three_points_spanning_full_great_circle_have_two_centers
    points = [point(0.0, 0.0), point(0.0, 120.0), point(0.0, 240.0)]
    result = cap(points)

    assert_in_delta 90.0, result.radius_degrees, 1e-9
    assert_equal [0, 1, 2], result.support_ids
    assert !result.unique?
    assert_contains result, points
    assert_boundary result, points, [0, 1, 2]
  end

  def test_great_circle_with_blocking_pole_is_unique
    points = [point(0.0, 0.0), point(0.0, 120.0),
              point(0.0, 240.0), point(90.0, 0.0)]
    result = cap(points)

    assert_in_delta 90.0, result.center.lat, 1e-9
    assert_in_delta 0.0, result.center.lng, 1e-9
    assert_in_delta 90.0, result.radius_degrees, 1e-9
    assert_equal [0, 1, 2], result.support_ids
    assert result.unique?
    assert_contains result, points
  end

  def test_set_requiring_more_than_hemisphere
    southern_latitude = -19.47122063449069
    points = [point(90.0, 0.0),
              point(southern_latitude, 0.0),
              point(southern_latitude, 120.0),
              point(southern_latitude, 240.0)]
    result = cap(points)

    assert result.angular_radius > Math::PI / 2.0
    assert_in_delta 109.47122063449069, result.radius_degrees, 1e-9
    assert_equal 3, result.support_ids.length
    assert !result.unique?
    assert_contains result, points
    assert_boundary result, points, result.support_ids
  end

  def test_shuffling_inputs_preserves_optimal_radius
    points = [point(0.0, 0.0), point(0.0, 120.0), point(0.0, 240.0),
              point(10.0, 20.0), point(-20.0, 70.0), point(30.0, -90.0)]

    radii = Array.new(20) do
      cap(points.shuffle).angular_radius
    end

    radii.each do |radius|
      assert_in_delta cap(points).angular_radius, radius, 1e-11
    end
  end

  def test_rejects_empty_input
    error = assert_raise(ArgumentError) { cap([]) }
    assert_equal 'expected at least one point', error.message
  end

  def test_rejects_mismatched_ids
    error = assert_raise(ArgumentError) do
      cap([point(0, 0), point(1, 1)], ids: [:only])
    end
    assert_equal 'ids must contain one entry per point', error.message
  end
end
