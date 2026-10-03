# encoding: utf-8
require File.join(File.dirname(__FILE__), 'helper')

class FrechetTest < Test::Unit::TestCase #:nodoc: all
  # Caller-chosen local "projection": treats lat/lng degrees as planar
  # kilometers. The library never picks this for us.
  PLANAR = ->(ll) { [ll.lng * 1000.0, ll.lat * 1000.0] }

  def ll(x, y)
    Geokit::LatLng.new(y, x)
  end

  def subdivide(path, pieces)
    points = []
    path.each_cons(2) do |p, q|
      pieces.times do |k|
        points << Geokit::LatLng.new(
          p.lat + (q.lat - p.lat) * k.to_f / pieces,
          p.lng + (q.lng - p.lng) * k.to_f / pieces
        )
      end
    end
    points << path.last
    points
  end

  def point_at(path, arc_length)
    remaining = arc_length
    path.each_cons(2) do |p, q|
      segment = Math.hypot((q.lng - p.lng) * 1000.0, (q.lat - p.lat) * 1000.0)
      next if segment.zero?

      if remaining <= segment
        ratio = remaining / segment
        return [
          (p.lng + (q.lng - p.lng) * ratio) * 1000.0,
          (p.lat + (q.lat - p.lat) * ratio) * 1000.0
        ]
      end
      remaining -= segment
    end
    [path.last.lng * 1000.0, path.last.lat * 1000.0]
  end

  def path_length(path)
    length = 0.0
    path.each_cons(2) do |p, q|
      length += Math.hypot((q.lng - p.lng) * 1000.0, (q.lat - p.lat) * 1000.0)
    end
    length
  end

  def vertex_hausdorff(path_a, path_b)
    [path_a, path_b].each_with_index.map do |src, i|
      dst = i.zero? ? path_b : path_a
      src.map do |p|
        dst.map do |q|
          Math.hypot((q.lng - p.lng) * 1000.0, (q.lat - p.lat) * 1000.0)
        end.min
      end.max
    end.max
  end

  def assert_witness_valid(result, path_a, path_b)
    witness = result.witness
    assert_equal [0.0, 0.0], witness.first, 'witness must start at [0, 0]'
    assert_in_delta path_length(path_a), witness.last[0], 1e-6
    assert_in_delta path_length(path_b), witness.last[1], 1e-6
    witness.each_cons(2) do |(s0, t0), (s1, t1)|
      assert s1 >= s0 - 1e-9, "witness moves backwards on first curve: #{s0} -> #{s1}"
      assert t1 >= t0 - 1e-9, "witness moves backwards on second curve: #{t0} -> #{t1}"
    end
    witness.each do |s, t|
      pa = point_at(path_a, s)
      pb = point_at(path_b, t)
      leash = Math.hypot(pa[0] - pb[0], pa[1] - pb[1])
      assert leash <= result.upper_bound + 1e-6,
             "witness leash #{leash} exceeds upper bound #{result.upper_bound}"
    end
  end

  def test_continuous_not_discrete_vertex_dynamic_programming
    # A vertex-to-vertex DP must pair Q's middle vertex (2, 1) with an
    # endpoint of P, costing sqrt(5). The continuous distance is 1:
    # the middle of P's segment pairs with (2, 1).
    path_a = [ll(0, 0), ll(4, 0)]
    path_b = [ll(0, 1), ll(2, 1), ll(4, 1)]
    result = Geokit::Frechet.distance(path_a, path_b, projection: PLANAR)
    assert_in_delta 1000.0, result.distance, 1e-3
    assert result.upper_bound < Math.sqrt(5) * 1000.0 - 1.0,
           'must beat the discrete vertex-to-vertex answer'
    assert_witness_valid(result, path_a, path_b)
  end

  def test_subdividing_a_segment_does_not_change_distance
    path_a = [ll(0, 0), ll(3, 4), ll(6, 0)]
    path_b = [ll(0, 1), ll(3, 5), ll(6, 1)]
    coarse = Geokit::Frechet.distance(path_a, path_b, projection: PLANAR)
    fine = Geokit::Frechet.distance(subdivide(path_a, 50), subdivide(path_b, 7),
                                    projection: PLANAR)
    assert_in_delta coarse.distance, fine.distance,
                    coarse.precision + fine.precision
    assert_witness_valid(fine, subdivide(path_a, 50), subdivide(path_b, 7))
  end

  def test_two_vehicles_same_route_different_sampling_density
    # The two vehicles drive the same geometric route; one logs GPS
    # points far more densely than the other.
    sparse = [ll(0, 0), ll(2, 3), ll(5, 3), ll(8, 0)]
    dense = subdivide(sparse, 40)
    result = Geokit::Frechet.distance(sparse, dense, projection: PLANAR)
    assert_in_delta 0.0, result.distance, 1e-3
    assert_witness_valid(result, sparse, dense)
  end

  def test_same_point_set_different_direction_is_not_hausdorff
    # Identical point sets (even identical traced curves), opposite
    # direction: Hausdorff distance is 0, Frechet is not, because
    # neither walk may go backwards.
    forward = [ll(0, 0), ll(4, 0), ll(10, 0)]
    backward = forward.reverse
    assert_in_delta 0.0, vertex_hausdorff(forward, backward), 1e-9
    result = Geokit::Frechet.distance(forward, backward, projection: PLANAR)
    assert_in_delta 10_000.0, result.distance, 1e-3
    assert_witness_valid(result, forward, backward)
  end

  def test_same_vertices_different_order_is_not_hausdorff
    path_a = [ll(0, 0), ll(2, 2), ll(4, 0)]
    path_b = [ll(0, 0), ll(4, 0), ll(2, 2)]
    assert_in_delta 0.0, vertex_hausdorff(path_a, path_b), 1e-9
    result = Geokit::Frechet.distance(path_a, path_b, projection: PLANAR)
    assert result.distance > 100.0, 'order of traversal must matter'
    assert_witness_valid(result, path_a, path_b)
  end

  def test_empty_trajectory_raises
    assert_raise(ArgumentError) do
      Geokit::Frechet.distance([], [ll(0, 0)], projection: PLANAR)
    end
    assert_raise(ArgumentError) do
      Geokit::Frechet.distance([ll(0, 0)], [], projection: PLANAR)
    end
  end

  def test_projection_is_chosen_by_caller_not_guessed
    assert_raise(ArgumentError) do
      Geokit::Frechet.distance([ll(0, 0)], [ll(0, 0)])
    end
    assert_raise(ArgumentError) do
      Geokit::Frechet.distance([ll(0, 0)], [ll(0, 0)], projection: :not_callable)
    end
    assert_raise(ArgumentError) do
      Geokit::Frechet.distance([ll(0, 0)], [ll(0, 0)], projection: ->(_p) { :junk })
    end
  end

  def test_zero_length_segments_and_repeated_points
    with_duplicates = [ll(0, 0), ll(0, 0), ll(1, 0), ll(1, 0), ll(1, 0)]
    plain = [ll(0, 0), ll(1, 0)]
    result = Geokit::Frechet.distance(with_duplicates, plain, projection: PLANAR)
    assert_in_delta 0.0, result.distance, 1e-3
    assert_witness_valid(result, with_duplicates, plain)
  end

  def test_curve_collapsed_to_a_point_by_zero_length_segments
    collapsed = [ll(1, 1), ll(1, 1), ll(1, 1)]
    path = [ll(0, 0), ll(3, 4)]
    from_collapsed = Geokit::Frechet.distance(collapsed, path, projection: PLANAR)
    from_point = Geokit::Frechet.distance([ll(1, 1)], path, projection: PLANAR)
    assert_in_delta from_point.distance, from_collapsed.distance, 1e-3
  end

  def test_single_point_trajectories
    both = Geokit::Frechet.distance([ll(0, 0)], [ll(3, 4)], projection: PLANAR)
    assert_in_delta 5000.0, both.distance, 1e-9
    assert_equal 0.0, both.precision
    assert_equal [[0.0, 0.0]], both.witness

    point_vs_line = Geokit::Frechet.distance([ll(0, 0)], [ll(0, 0), ll(3, 4)],
                                             projection: PLANAR)
    assert_in_delta 5000.0, point_vs_line.distance, 1e-9
    assert_equal [0.0, 0.0], point_vs_line.witness.first
    assert_in_delta 5000.0, point_vs_line.witness.last[1], 1e-9
    assert_equal 0.0, point_vs_line.witness.last[0]

    line_vs_point = Geokit::Frechet.distance([ll(0, 0), ll(3, 4)], [ll(0, 0)],
                                             projection: PLANAR)
    assert_in_delta 5000.0, line_vs_point.distance, 1e-9
  end

  def test_identical_points_give_zero
    result = Geokit::Frechet.distance([ll(1, 1), ll(1, 1)], [ll(1, 1), ll(1, 1)],
                                      projection: PLANAR)
    assert_in_delta 0.0, result.distance, 1e-9
  end

  def test_reports_certified_bounds_and_precision
    path_a = [ll(0, 0), ll(3, 4), ll(6, 0)]
    path_b = [ll(0, 1), ll(3, 5), ll(6, 1)]
    result = Geokit::Frechet.distance(path_a, path_b, projection: PLANAR,
                                      tolerance: 1e-3)
    assert result.lower_bound <= result.distance
    assert result.distance <= result.upper_bound
    assert_in_delta result.upper_bound - result.lower_bound, result.precision, 1e-12
    assert result.precision <= 1e-3, "precision #{result.precision} exceeds tolerance"
    # The bracket is honest: it contains the true continuous distance,
    # which the witness achieves up to the upper bound.
    assert_witness_valid(result, path_a, path_b)
  end

  def test_equirectangular_projection_helper
    origin = Geokit::LatLng.new(37.0, -122.0)
    projection = Geokit::Frechet.equirectangular_projection(origin)
    assert_equal [0.0, 0.0], projection.call(origin)
    _x, y = projection.call(Geokit::LatLng.new(38.0, -122.0))
    assert_in_delta 111_295.68, y, 1.0
    x, = projection.call(Geokit::LatLng.new(37.0, -121.0))
    assert_in_delta 111_295.68 * Math.cos(37.0 * Math::PI / 180.0), x, 1.0
  end

  def test_works_with_caller_supplied_projection_over_lat_lng
    origin = Geokit::LatLng.new(37.0, -122.0)
    projection = Geokit::Frechet.equirectangular_projection(origin)
    path_a = [Geokit::LatLng.new(37.0, -122.0), Geokit::LatLng.new(37.001, -122.0)]
    path_b = [Geokit::LatLng.new(37.0, -122.0), Geokit::LatLng.new(37.001, -122.001)]
    result = Geokit::Frechet.distance(path_a, path_b, projection: projection)
    assert_in_delta 88.9, result.distance, 1.0
    assert_witness_valid_lat_lng(result, path_a, path_b, projection)
  end

  def assert_witness_valid_lat_lng(result, path_a, path_b, projection)
    projected_a = path_a.map { |p| projection.call(p) }
    projected_b = path_b.map { |p| projection.call(p) }
    result.witness.each do |s, t|
      pa = planar_point_at(projected_a, s)
      pb = planar_point_at(projected_b, t)
      leash = Math.hypot(pa[0] - pb[0], pa[1] - pb[1])
      assert leash <= result.upper_bound + 1e-6
    end
  end

  def planar_point_at(points, arc_length)
    remaining = arc_length
    points.each_cons(2) do |p, q|
      segment = Math.hypot(q[0] - p[0], q[1] - p[1])
      next if segment.zero?

      if remaining <= segment
        ratio = remaining / segment
        return [p[0] + (q[0] - p[0]) * ratio, p[1] + (q[1] - p[1]) * ratio]
      end
      remaining -= segment
    end
    points.last
  end
end
