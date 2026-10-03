require File.join(File.dirname(__FILE__), 'helper')

class ContinuousFrechetTest < Test::Unit::TestCase
  PROJECTION = ->(point) { [point.lng.to_f, point.lat.to_f] }

  def points(*coordinates)
    coordinates.map { |x, y| Geokit::LatLng.new(y, x) }
  end

  def assert_monotonic_witness(result)
    assert result.witness.monotonic?
    assert_equal result.witness.pairs.first.first, 0.0
    assert_equal result.witness.pairs.first.last, 0.0
    assert_in_delta result.witness.arc_length_a,
                    result.witness.pairs.last.first, 1e-9
    assert_in_delta result.witness.arc_length_b,
                    result.witness.pairs.last.last, 1e-9
  end

  def assert_witness_distance(result, path_a, path_b, expected,
                             extra_points = [])
    projected_a = path_a.map(&PROJECTION)
    projected_b = path_b.map(&PROJECTION)
    samples = result.witness.pairs + extra_points

    samples.each do |arc_a, arc_b|
      point_a = point_at(projected_a, arc_a, result.witness.arc_length_a)
      point_b = point_at(projected_b, arc_b, result.witness.arc_length_b)
      distance = Math.hypot(point_a[0] - point_b[0],
                            point_a[1] - point_b[1])
      assert distance <= result.upper_bound + 1e-9,
             "witness distance #{distance} exceeds #{expected}"
    end
  end

  def point_at(path, arc_length, total_length)
    return path.first if total_length.zero?

    remaining = arc_length
    path.each_cons(2) do |first, second|
      segment_length = Math.hypot(second[0] - first[0],
                                  second[1] - first[1])
      if remaining <= segment_length
        ratio = segment_length.zero? ? 0.0 : remaining / segment_length
        return [
          first[0] + (second[0] - first[0]) * ratio,
          first[1] + (second[1] - first[1]) * ratio,
        ]
      end
      remaining -= segment_length
    end
    path.last
  end

  def midpoint_witness_samples(result)
    result.witness.pairs.each_cons(2).map do |first, second|
      [
        (first[0] + second[0]) / 2.0,
        (first[1] + second[1]) / 2.0,
      ]
    end
  end

  def test_parallel_segments_are_continuously_evaluated
    path_a = points([0, 0], [10, 0])
    path_b = points([0, 1], [10, 1])

    result = Geokit::ContinuousFrechet.distance_between(
      path_a, path_b, PROJECTION, absolute_tolerance: 1e-7
    )

    assert_in_delta 1.0, result.distance, 1e-7
    assert result.lower_bound <= 1.0
    assert result.upper_bound >= 1.0
    assert result.upper_bound - result.lower_bound < 1e-7
    assert_monotonic_witness(result)
    assert_witness_distance(result, path_a, path_b, 1.0,
                            midpoint_witness_samples(result))
  end

  def test_arbitrary_subdivision_does_not_change_distance
    path_a = points([0, 0], [10, 0])
    path_b = points([0, 1], [10, 1])
    subdivided_a = points([0, 0], [2, 0], [7, 0], [10, 0])
    subdivided_b = points([0, 1], [4, 1], [5, 1], [9, 1], [10, 1])

    result = Geokit::ContinuousFrechet.distance_between(
      path_a, path_b, PROJECTION, absolute_tolerance: 1e-7
    )
    subdivided_result = Geokit::ContinuousFrechet.distance_between(
      subdivided_a, subdivided_b, PROJECTION, absolute_tolerance: 1e-7
    )

    assert_in_delta result.distance, subdivided_result.distance, 1e-6
    assert_in_delta 1.0, subdivided_result.distance, 1e-6
    assert_monotonic_witness(subdivided_result)
    assert_witness_distance(subdivided_result, subdivided_a,
                            subdivided_b, 1.0,
                            midpoint_witness_samples(subdivided_result))
  end

  def test_reordered_equal_point_set_is_not_hausdorff_distance
    forward = points([0, 0], [10, 0], [10, 10])
    reordered = points([10, 10], [0, 0], [10, 0])

    result = Geokit::ContinuousFrechet.distance_between(
      forward, reordered, PROJECTION, absolute_tolerance: 1e-7
    )

    assert_in_delta 10 * Math.sqrt(2), result.lower_bound, 1e-6
    assert result.lower_bound > 1.0
    assert_monotonic_witness(result)
  end

  def test_interior_match_beats_discrete_vertex_dynamic_programming
    path_a = points([0, 0], [10, 0])
    path_b = points([0, 0], [5, 0.1], [10, 0])

    result = Geokit::ContinuousFrechet.distance_between(
      path_a, path_b, PROJECTION, absolute_tolerance: 1e-7
    )

    assert_in_delta 0.1, result.distance, 1e-7
    assert_monotonic_witness(result)
    assert_witness_distance(result, path_a, path_b, 0.1,
                            midpoint_witness_samples(result))
  end

  def test_repeated_points_and_zero_length_segments
    path_a = points([0, 0], [0, 0], [10, 0], [10, 0])
    path_b = points([0, 1], [5, 1], [5, 1], [10, 1])

    result = Geokit::ContinuousFrechet.distance_between(
      path_a, path_b, PROJECTION, absolute_tolerance: 1e-7
    )

    assert_in_delta 1.0, result.distance, 1e-7
    assert_monotonic_witness(result)
  end

  def test_subdivision_is_invariant_on_both_trajectories
    base_a = points([0, 0], [10, 0])
    base_b = points([0, 1], [10, 1])
    subdivided_a = points([0, 0], [1, 0], [4, 0], [8, 0], [10, 0])
    subdivided_b = points([0, 1], [3, 1], [6, 1], [10, 1])

    expected = Geokit::ContinuousFrechet.distance_between(
      base_a, base_b, PROJECTION, absolute_tolerance: 1e-7
    )
    result = Geokit::ContinuousFrechet.distance_between(
      subdivided_a, subdivided_b, PROJECTION,
      absolute_tolerance: 1e-7
    )

    assert_in_delta expected.distance, result.distance, 1e-6
    assert_monotonic_witness(result)
  end

  def test_single_point_trajectory
    segment = points([0, 0], [10, 0])
    point = points([5, 5])

    result = Geokit::ContinuousFrechet.distance_between(
      point, segment, PROJECTION
    )

    assert_in_delta 5.0 * Math.sqrt(2), result.distance, 1e-9
    assert result.exact
    assert_monotonic_witness(result)
  end

  def test_collapsed_repeated_trajectory
    collapsed = points([2, 0], [2, 0], [2, 0])
    segment = points([0, 0], [5, 0], [10, 0])

    result = Geokit::ContinuousFrechet.distance_between(
      collapsed, segment, PROJECTION
    )

    assert_equal 8.0, result.distance
    assert result.exact
    assert(result.witness.pairs.all? { |arc_a, _| arc_a.zero? })
    assert_monotonic_witness(result)
  end

  def test_empty_trajectories
    both_empty = Geokit::ContinuousFrechet.distance_between(
      [], [], PROJECTION
    )
    one_empty = Geokit::ContinuousFrechet.distance_between(
      [Geokit::LatLng.new(0, 0)], [], PROJECTION
    )

    assert_equal 0.0, both_empty.distance
    assert both_empty.exact
    assert_equal Float::INFINITY, one_empty.distance
    assert one_empty.exact
  end

  def test_requires_projection
    assert_raise ArgumentError do
      Geokit::ContinuousFrechet.distance_between(
        [Geokit::LatLng.new(0, 0)], [Geokit::LatLng.new(1, 1)]
      )
    end
  end

  def test_requires_projection_to_return_two_metric_coordinates
    assert_raise ArgumentError do
      Geokit::ContinuousFrechet.distance_between(
        [Geokit::LatLng.new(0, 0)], [Geokit::LatLng.new(1, 1)],
        ->(_point) { [0, 0, 1] }
      )
    end
  end

  def test_block_projection_and_returned_bounds
    result = Geokit::ContinuousFrechet.distance_between(
      points([0, 0], [10, 0]),
      points([0, 2], [10, 2]),
      nil,
      absolute_tolerance: 1e-7
    ) { |point| [point.lng, point.lat] }

    assert_in_delta 2.0, result.distance, 1e-7
    assert result.lower_bound <= result.distance
    assert result.upper_bound >= result.distance
    assert result.absolute_error >= 0.0
    assert result.iterations > 0
  end
end
