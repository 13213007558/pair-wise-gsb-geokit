require File.join(File.dirname(__FILE__), 'helper')

class SphericalCapTest < Test::Unit::TestCase
  def test_single_point
    point = [10, 20]
    cap = Geokit::LatLng.minimal_spherical_cap(point, units: :kms)

    assert_equal 0.0, cap.angular_radius
    assert_equal [10.0, 20.0], cap.center.to_a
    assert_equal [0], cap.support_indices
    assert_equal [0], cap.boundary_indices
    assert cap.unique_center?
  end

  def test_duplicate_points
    points = [[10, 0], [10, 0], [-11, 0]]
    cap = Geokit::LatLng.minimal_spherical_cap(points, units: :kms)

    assert_in_degrees 10.5, cap.angular_radius_degrees
    assert_equal [0, 2], cap.support_indices
    assert_equal [0, 1, 2], cap.boundary_indices
    assert_cover cap, points
  end

  def test_two_points
    points = [[0, -45], [0, 45]]
    cap = Geokit::LatLng.minimal_spherical_cap(points, units: :kms)

    assert_in_degrees 45, cap.angular_radius_degrees
    assert_equal [0, 0], cap.center.to_a
    assert_equal [0, 1], cap.support_indices
    assert_cover cap, points
  end

  def test_antipodal_points_are_a_nonunique_hemisphere
    points = [[0, 0], [0, 180]]
    cap = Geokit::LatLng.minimal_spherical_cap(points, units: :kms)

    assert_in_degrees 90, cap.angular_radius_degrees
    assert_equal [0, 1], cap.support_indices
    assert_equal [0, 1], cap.boundary_indices
    assert cap.nonunique_center?
    assert_equal :hemisphere_center_continuum, cap.nonunique_reason
    assert_cover cap, points
  end

  def test_three_points_on_same_great_circle
    points = [[0, 0], [0, 120], [0, 240]]
    cap = Geokit::LatLng.minimal_spherical_cap(points, units: :kms)

    assert_in_degrees 90, cap.angular_radius_degrees
    assert_equal 2, cap.support_indices.size
    assert_equal [0, 1, 2], cap.boundary_indices
    assert cap.nonunique_center?
    assert_equal :opposite_hemisphere_centers, cap.nonunique_reason
    assert_cover cap, points
  end

  def test_unique_hemisphere_with_two_boundary_directions
    points = [[0, 0], [0, 180], [0, 90], [45, 0]]
    cap = Geokit::LatLng.minimal_spherical_cap(points, units: :kms)

    assert_in_degrees 90, cap.angular_radius_degrees
    assert_equal [0, 1, 2], cap.boundary_indices
    assert cap.unique_center?
    assert_cover cap, points
  end

  def test_three_noncollinear_supports
    points = [[0, 0], [0, 90], [90, 0]]
    cap = Geokit::LatLng.minimal_spherical_cap(points, units: :kms)

    assert_in_degrees 54.735610317245346, cap.angular_radius_degrees
    assert_equal [0, 1, 2], cap.support_indices
    assert_equal [0, 1, 2], cap.boundary_indices
    assert cap.unique_center?
    assert_cover cap, points
  end

  def test_cap_larger_than_hemisphere
    points = rotated_tetrahedron_points
    cap = Geokit::LatLng.minimal_spherical_cap(points, units: :kms)

    assert cap.angular_radius > Math::PI / 2
    assert_in_degrees 109.47122063449069, cap.angular_radius_degrees
    assert_equal 3, cap.support_indices.size
    assert_equal 3, cap.boundary_indices.size
    assert cap.nonunique_center?
    assert_equal :multiple_spherical_caps, cap.nonunique_reason
    assert_cover cap, points
  end

  def test_antipodal_set_requiring_cap_larger_than_hemisphere
    points = [[0, 0], [0, 180], [90, 0], [-90, 0], [0, 90], [0, -90]]
    cap = Geokit::LatLng.minimal_spherical_cap(points, units: :kms)

    assert cap.angular_radius > Math::PI / 2
    assert_in_degrees 125.26438968275465, cap.angular_radius_degrees
    assert_equal 3, cap.support_indices.size
    assert_equal 3, cap.boundary_indices.size
    assert cap.nonunique_center?
    assert_cover cap, points
  end

  def test_radius_unit_conversions
    points = [[0, -45], [0, 45]]
    cap = Geokit::LatLng.minimal_spherical_cap(points, units: :kms)

    earth_radius_kms = Geokit::LatLng.units_sphere_multiplier(:kms)
    assert_in_delta cap.angular_radius * earth_radius_kms, cap.radius, 1e-9
    assert_in_delta(
      cap.angular_radius * Geokit::LatLng.units_sphere_multiplier(:miles),
      cap.radius_in_units(:miles),
      1e-9
    )
    assert_in_delta(
      Geokit::LatLng.distance_between(points[0], cap.center, units: :meters),
      cap.radius_in_units(:meters),
      0.001
    )
  end

  def test_order_does_not_change_optimal_radius
    points = rotated_tetrahedron_points
    shuffled = points.shuffle(random: Random.new(17))
    original = Geokit::LatLng.minimal_spherical_cap(points, units: :kms)
    reordered = Geokit::LatLng.minimal_spherical_cap(shuffled, units: :kms)

    assert_in_delta original.angular_radius, reordered.angular_radius, 1e-10
    assert_equal original.center, reordered.center
    assert_cover reordered, shuffled

    expected_support = original.support_indices.map { |index| points[index] }.sort
    actual_support = reordered.support_indices.map { |index| shuffled[index] }.sort
    assert_equal expected_support, actual_support
  end

  def test_small_sets_are_checked_by_independent_candidate_enumeration
    [
      [[10, 20]],
      [[10, 0], [10, 0], [-11, 0]],
      [[0, -45], [0, 45]],
      [[0, 0], [0, 180]],
      [[0, 0], [0, 120], [0, 240]],
      [[0, 0], [0, 180], [0, 90], [45, 0]],
      [[0, 0], [0, 90], [90, 0]],
      rotated_tetrahedron_points,
      [[0, 0], [0, 180], [90, 0], [-90, 0], [0, 90], [0, -90]]
    ].each do |points|
      cap = Geokit::LatLng.minimal_spherical_cap(points, units: :kms)
      reference_radius = brute_force_minimum_radius(points)

      assert_in_delta reference_radius, cap.angular_radius, 1e-9
      assert_cover cap, points
      cap.support_indices.each do |index|
        distance = cap.angular_distance(points[index])
        assert_in_delta cap.angular_radius, distance, 1e-8
      end
    end
  end

  private

  def assert_in_degrees(expected, actual)
    assert_in_delta expected, actual, 1e-9
  end

  def assert_cover(cap, points)
    points.each do |point|
      assert cap.contains?(point)
      assert cap.angular_distance(point) <= cap.angular_radius + 1e-9
    end
    cap.boundary_indices.each do |index|
      assert_in_delta cap.angular_radius,
                      cap.angular_distance(points[index]), 1e-8
    end
  end

  def rotated_tetrahedron_points
    vectors = [[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]]
    vectors.map do |x, y, z|
      [(x - z) / Math.sqrt(2), y, (x + z) / Math.sqrt(2)]
    end.map { |vector| point_from_test_vector(vector) }
  end

  def brute_force_minimum_radius(points)
    vectors = points.map { |point| test_vector(point) }.uniq { |v| v.map { |x| x.round(10) } }
    candidates = []

    vectors.each { |v| candidates << unit(v) }
    vectors.combination(2).each do |a, b|
      midpoint = unit(a.zip(b).map { |x, y| x + y })
      candidates.concat([midpoint, midpoint.map { |x| -x }]) if midpoint
      if dot3(a, b) <= -1 + 1e-10
        axis = [[1, 0, 0], [0, 1, 0], [0, 0, 1]].max_by do |candidate|
          (dot3(a, candidate).abs - 1).abs
        end
        perpendicular = unit(cross3(a, axis))
        candidates.concat([perpendicular, perpendicular.map { |x| -x }])
      end
      normal = unit(cross3(a, b))
      candidates.concat([normal, normal.map { |x| -x }]) if normal
    end
    vectors.combination(3).each do |a, b, c|
      normal = unit(cross3(b.zip(a).map { |x, y| x - y },
                          c.zip(a).map { |x, y| x - y }))
      candidates.concat([normal, normal.map { |x| -x }]) if normal
    end

    candidates.map do |center|
      vectors.map { |point| Math.acos([[dot3(center, point), -1].max, 1].min) }.max
    end.min
  end

  def test_vector(point)
    latitude = point[0] * Math::PI / 180
    longitude = point[1] * Math::PI / 180
    [
      Math.cos(latitude) * Math.cos(longitude),
      Math.cos(latitude) * Math.sin(longitude),
      Math.sin(latitude)
    ]
  end

  def point_from_test_vector(vector)
    vector = unit(vector)
    [
      Math.asin(vector[2]) * 180 / Math::PI,
      Math.atan2(vector[1], vector[0]) * 180 / Math::PI
    ]
  end

  def unit(vector)
    magnitude = Math.sqrt(dot3(vector, vector))
    return if magnitude < 1e-12
    vector.map { |value| value / magnitude }
  end

  def dot3(first, second)
    first.zip(second).map { |x, y| x * y }.inject(:+)
  end

  def cross3(first, second)
    [
      first[1] * second[2] - first[2] * second[1],
      first[2] * second[0] - first[0] * second[2],
      first[0] * second[1] - first[1] * second[0]
    ]
  end
end
