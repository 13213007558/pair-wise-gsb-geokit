require File.expand_path(File.dirname(__FILE__) + '/helper')

class SpatialIndexTest < Test::Unit::TestCase
  def test_bulk_load_and_duplicate_coordinates
    index = Geokit::SpatialIndex.new(
      'c' => [10, 20],
      'a' => [10, 20],
      'b' => [11, 20]
    )

    assert_equal %w[a c b], index.nearest([10, 20], 3).map(&:id)
    assert_equal [10.0, 20.0], index.location_for('a').to_a
  end

  def test_insert_update_and_delete_do_not_leave_old_point
    index = Geokit::SpatialIndex.new([['one', [0, 0]]])
    index.add('two', 10, 10)
    index.add('one', 20, 20)

    assert_equal 2, index.size
    assert_equal [20.0, 20.0], index.location_for('one').to_a
    assert_equal %w[one],
                 index.within_radius([20, 20], 1, units: :kms).map(&:id)

    assert index.delete('two')
    assert_equal 1, index.size
    assert !index.include?('two')
    assert_equal false, index.delete('two')
  end

  def test_delete_last_point_at_shared_coordinate
    index = Geokit::SpatialIndex.new([['a', [0, 0]], ['b', [0, 0]]])
    index.delete('a')

    assert_equal ['b'], index.nearest([0, 0], 2).map(&:id)
    assert_equal 1, index.size
  end

  def test_equal_distances_are_stable_by_id
    index = Geokit::SpatialIndex.new(
      'z' => [0, 1],
      'y' => [0, -1],
      'x' => [1, 0]
    )

    assert_equal %w[x y z], index.nearest([0, 0], 3).map(&:id)
    assert_equal %w[x y z], index.within_radius([0, 0], 100).map(&:id)
  end

  def test_radius_uses_existing_sphere_distance
    index = Geokit::SpatialIndex.new([['a', [1, 0]]], units: :kms)
    expected = Geokit::LatLng.distance_between(
      [0, 0], [1, 0], units: :kms, formula: :sphere
    )

    assert_equal [], index.within_radius([0, 0], expected - 0.001).map(&:id)
    assert_equal ['a'], index.within_radius([0, 0], expected).map(&:id)
    assert_in_delta expected, index.nearest([0, 0]).first.distance, 1e-9
  end

  def test_antimeridian_poles_and_huge_radius
    points = {
      north: [90, 45],
      south: [-90, -120],
      east: [0, 179.999999],
      west: [0, -179.999999],
      middle: [45, 0]
    }
    index = Geokit::SpatialIndex.new(points, units: :kms)

    assert_equal %i[east west],
                 index.within_radius([0, 180], 0.2).map(&:id)
    assert_equal :north, index.nearest([89.999999, 100], 1).first.id
    assert_equal :south, index.nearest([-89.999999, 100], 1).first.id

    earth_circumference = 2 * Math::PI *
                          Geokit::LatLng.units_sphere_multiplier(:kms)
    assert_equal points.keys.sort,
                 index.within_radius(
                   [0, 0], earth_circumference, units: :kms
                 ).map(&:id).sort
  end

  def test_zero_and_antipodal_radius_do_not_lose_numeric_matches
    points = {
      same: [0, 0],
      nearby: [1e-12, 1e-12],
      antipode: [0, 180],
      other: [10, 20]
    }
    index = Geokit::SpatialIndex.new(points, units: :meters)

    assert_equal %i[nearby same],
                 index.within_radius([0, 0], 0).map(&:id)
    antipodal_distance = Geokit::LatLng.distance_between(
      [0, 0], [0, 180], units: :meters, formula: :sphere
    )
    assert_include index.within_radius(
      [0, 0], antipodal_distance
    ).map(&:id), :antipode
  end

  def test_nearest_with_fewer_than_k_points
    index = Geokit::SpatialIndex.new({ a: [0, 0], b: [1, 1] })

    assert_equal %i[a b], index.nearest([0, 0], 5).map(&:id)
  end

  def test_updating_shared_coordinate_removes_old_reference
    index = Geokit::SpatialIndex.new(
      [['a', [0, 0]], ['b', [0, 0]]]
    )
    index.add('a', 45, 45)

    assert_equal [[0.0, 0.0]],
                 index.within_radius([0, 0], 1, units: :kms)
                   .map(&:location).map(&:to_a)
    assert_equal ['b'], index.within_radius([0, 0], 1, units: :kms)
                                .map(&:id)
  end

  def test_random_operation_sequence_matches_brute_force
    random = Random.new(742_391)
    index = Geokit::SpatialIndex.new([], units: :meters)
    reference = {}

    500.times do
      id = random.rand(25)
      location = random_point(random)
      case random.rand(4)
      when 0, 1
        index.add(id, location)
        reference[id] = location
      when 2
        index.delete(id)
        reference.delete(id)
      else
        point = random_point(random)
        if random.rand(2).zero?
          k = random.rand(1..6)
          assert_brute_force_knn(index, reference, point, k)
        else
          radius = random.rand(0.0..22_000_000)
          assert_brute_force_radius(index, reference, point, radius)
        end
      end
    end
  end

  def test_large_query_reports_distance_calculations
    random = Random.new(426_426)
    entries = Array.new(10_000) do |id|
      [id, [random.rand(-90.0..90.0), random.rand(-180.0...180.0)]]
    end
    index = Geokit::SpatialIndex.new(entries, units: :kms)
    origin = [35.0, -95.0]
    results = index.nearest(origin, 10)
    expected = entries.map do |id, location|
      Geokit::SpatialIndex::Result.new(
        id,
        Geokit::LatLng.new(*location),
        Geokit::LatLng.distance_between(
          origin, location, units: :kms, formula: :sphere
        )
      )
    end.sort_by { |result| [result.distance, result.id] }.first(10)

    assert_equal 10, results.size
    assert_equal expected.map(&:id), results.map(&:id)
    expected.zip(results) do |expected_result, actual_result|
      assert_in_delta expected_result.distance,
                      actual_result.distance, 1e-9
    end
    assert_operator index.distance_calculation_count, :>, 0
    assert_operator index.distance_calculation_count, :<, entries.size
    recorded_count = index.distance_calculation_count
    index.nearest(origin, 10)
    assert_equal recorded_count, index.distance_calculation_count
  end

  private

  def random_point(random)
    [
      random.rand(-90.0..90.0),
      random.rand(-180.0...180.0)
    ]
  end

  def brute_force_results(reference, point, units = :meters)
    reference.map do |id, location|
      distance = Geokit::LatLng.distance_between(
        point, location, units: units, formula: :sphere
      )
      [id, distance]
    end
  end

  def assert_brute_force_knn(index, reference, point, k)
    expected = brute_force_results(reference, point)
                       .sort_by { |id, distance| [distance, id] }
                       .first(k)
    actual = index.nearest(point, k).map { |result| [result.id, result.distance] }

    assert_equal expected.map(&:first), actual.map(&:first)
    expected.zip(actual) do |expected_pair, actual_pair|
      assert_in_delta expected_pair.last, actual_pair.last, 1e-9
    end
  end

  def assert_brute_force_radius(index, reference, point, radius)
    expected = brute_force_results(reference, point)
                       .select { |_, distance| distance <= radius }
                       .sort_by { |id, distance| [distance, id] }
    actual = index.within_radius(point, radius).map do |result|
      [result.id, result.distance]
    end

    assert_equal expected.map(&:first), actual.map(&:first)
    expected.zip(actual) do |expected_pair, actual_pair|
      assert_in_delta expected_pair.last, actual_pair.last, 1e-9
    end
  end
end
