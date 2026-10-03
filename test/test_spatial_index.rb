require File.join(File.dirname(__FILE__), 'helper')

class SpatialIndexTest < Test::Unit::TestCase
  def setup
    @index = Geokit::SpatialIndex.new(units: :kms)
  end

  def test_bulk_insert_and_lookup
    @index = Geokit::SpatialIndex.build(
      1 => [10.0, 20.0],
      2 => [-10.0, -20.0]
    )

    assert_equal 2, @index.size
    assert_equal Geokit::LatLng.new(10.0, 20.0), @index[1]
    assert @index.include?(2)
  end

  def test_update_removes_old_position_and_keeps_one_id
    @index.insert(1, [0.0, 1.0])
    @index.insert(1, [45.0, 45.0])

    assert_equal 1, @index.size
    assert_equal [[1, 0.0]], @index.within_radius(
      [45.0, 45.0], 1, units: :kms, with_distances: true
    )
    assert_empty @index.within_radius([0.0, 0.0], 1000, units: :kms)
  end

  def test_delete_removes_id
    @index.insert(1, [0.0, 0.0])
    @index.delete(1)

    assert_equal 0, @index.size
    assert !@index.include?(1)
    assert_empty @index.within_radius([0.0, 0.0], 10_000, units: :kms)
  end

  def test_duplicate_coordinates_and_tie_order
    1.upto(5) { |id| @index.insert(id, [0.0, 0.0]) }

    assert_equal [1, 2, 3, 4, 5], @index.within_radius([0.0, 0.0], 1, units: :kms)
    assert_equal [1, 2, 3], @index.nearest([1.0, 0.0], 3, units: :kms)
  end

  def test_results_are_sorted_by_distance_and_id
    points = {
      2 => [0.0, 0.2],
      1 => [0.0, 0.1],
      3 => [0.0, 0.1]
    }
    @index = Geokit::SpatialIndex.build(points, units: :kms)

    ids = @index.within_radius([0.0, 0.0], 50, units: :kms)
    assert_equal [1, 3, 2], ids

    nearest = @index.nearest([0.0, 0.0], 2, units: :kms)
    assert_equal [1, 3], nearest
  end

  def test_with_distances_uses_existing_sphere_model
    @index.insert(1, [0.0, 1.0])

    expected = Geokit::LatLng.distance_between(
      [0.0, 0.0], [0.0, 1.0], units: :kms, formula: :sphere
    )
    actual = @index.within_radius(
      [0.0, 0.0], 200, units: :kms, with_distances: true
    ).first.last

    assert_in_delta expected, actual, 1e-9
  end

  def test_antimeridian_radius_query
    @index.insert(1, [0.0, 179.95])
    @index.insert(2, [0.0, -179.95])
    @index.insert(3, [0.0, 0.0])

    assert_equal [1, 2], @index.within_radius(
      [0.0, -180.0], 20, units: :kms
    ).sort
  end

  def test_pole_queries
    points = {
      1 => [89.95, 0.0],
      2 => [89.95, 90.0],
      3 => [89.95, -90.0],
      4 => [-89.95, 0.0]
    }
    @index = Geokit::SpatialIndex.build(points, units: :kms)

    assert_equal [1, 2, 3], @index.within_radius(
      [90.0, 0.0], 10, units: :kms
    ).sort
    assert_equal [1, 2, 3], @index.nearest([90.0, 0.0], 3, units: :kms)
  end

  def test_very_large_radius_returns_every_result
    1.upto(20) { |id| @index.insert(id, [rand * 180 - 90, rand * 360 - 180]) }

    assert_equal (1..20).to_a, @index.within_radius(
      [0.0, 0.0], 40_000, units: :kms
    ).sort
  end

  def test_knn_returns_exact_ties
    @index.insert(1, [1.0, 0.0])
    @index.insert(2, [-1.0, 0.0])
    @index.insert(3, [0.0, 1.0])

    assert_equal [1, 2], @index.nearest([0.0, 0.0], 2, units: :kms)
  end

  def test_random_operation_sequence_matches_brute_force
    srand(20261003)
    points = {}

    300.times do |operation|
      id = rand(1..80)
      point = [rand * 180 - 90, wrap_longitude(rand * 720 - 360)]

      case operation % 4
      when 0, 1
        @index.insert(id, point)
        points[id] = point
      when 2
        @index.delete(id)
        points.delete(id)
      else
        @index.insert(id, point)
        points[id] = point
      end

      verify_query(points, [rand * 180 - 90, wrap_longitude(rand * 720 - 360)], rand * 1000)
      verify_knn(points, [rand * 180 - 90, wrap_longitude(rand * 720 - 360)], rand(1..8))
    end
  end

  def test_large_collection_prunes_without_approximating
    srand(4242)
    points = {}
    50_000.times do |id|
      points[id + 1] = [rand * 180 - 90, rand * 360 - 180]
    end
    index = Geokit::SpatialIndex.build(points, units: :kms)
    center = [20.0, 30.0]
    radius = 100.0

    expected = brute_force_radius(points, center, radius)
    actual = index.within_radius(center, radius, units: :kms)

    assert_equal expected, actual
    assert index.distance_computations < points.size
    assert_equal 11, index.distance_computations
    puts "large collection: #{points.size} points, " \
         "#{index.distance_computations} exact distance computations"
  end

  def test_cell_lower_bound_never_underestimates_sphere_distance
    srand(11)
    250.times do
      lat_min, lat_max = [rand(-90.0..90.0), rand(-90.0..90.0)].minmax
      lng_min, lng_max = [rand(-180.0...180.0), rand(-180.0...180.0)].minmax
      center = [rand(-90.0..90.0), rand(-180.0...180.0)]
      index = Geokit::SpatialIndex.new(units: :kms, lat_bins: 1, lng_bins: 1)
      center_vector = index.send(:vector_for, Geokit::LatLng.normalize(center))
      minimum_distance = nil

      12.times do |lat_step|
        12.times do |lng_step|
          lat = lat_min + ((lat_max - lat_min) * lat_step / 11.0)
          lng = lng_min + ((lng_max - lng_min) * lng_step / 11.0)
          vector = index.send(:vector_for, Geokit::LatLng.new(lat, lng))
          distance = index.send(:angular_distance_between, center_vector, vector)
          minimum_distance = distance if minimum_distance.nil? || distance < minimum_distance
        end
      end

      index.bulk_insert(
        1 => [lat_min, lng_min],
        2 => [lat_min, lng_max],
        3 => [lat_max, lng_min],
        4 => [lat_max, lng_max]
      )
      cell = index.instance_variable_get(:@cells).first
      lower_bound = index.send(:lower_angular_distance, center_vector, cell)

      assert lower_bound <= minimum_distance + 1e-8
    end
  end

  private

  def wrap_longitude(lng)
    (((lng + 180.0) % 360.0) + 360.0) % 360.0 - 180.0
  end

  def verify_query(points, center, radius)
    assert_equal brute_force_radius(points, center, radius),
                 @index.within_radius(center, radius, units: :kms)
  end

  def verify_knn(points, center, count)
    distances = points.map do |id, point|
      [distance(center, point), id]
    end
    expected = stable_pairs(distances).first(count).map(&:last)

    assert_equal expected, @index.nearest(center, count, units: :kms)
  end

  def brute_force_radius(points, center, radius)
    points.map { |id, point| [distance(center, point), id] }
          .select { |distance, _| distance <= radius + 1e-9 }
          .then { |pairs| stable_pairs(pairs) }
          .map(&:last)
  end

  def distance(from, to)
    Geokit::LatLng.distance_between(from, to, units: :kms, formula: :sphere)
  end

  def stable_pairs(pairs)
    pairs.sort do |left, right|
      difference = left.first - right.first
      next left.last <=> right.last if difference.abs <= 1e-12 * 6376.77271

      difference <=> 0
    end
  end
end
