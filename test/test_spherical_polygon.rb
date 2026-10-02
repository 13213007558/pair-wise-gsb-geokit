# encoding: utf-8
require File.join(File.dirname(__FILE__), 'helper')

class SphericalPolygonTest < Test::Unit::TestCase
  def setup
    @octant = [
      Geokit::LatLng.new(0, 0),
      Geokit::LatLng.new(0, 90),
      Geokit::LatLng.new(90, 0)
    ]
    @subdivided_octant = [
      Geokit::LatLng.new(0, 0),
      Geokit::LatLng.new(0, 45),
      Geokit::LatLng.new(0, 90),
      Geokit::LatLng.new(90, 0)
    ]
    @triangle_one = [
      Geokit::LatLng.new(0, 0),
      Geokit::LatLng.new(0, 45),
      Geokit::LatLng.new(90, 0)
    ]
    @triangle_two = [
      Geokit::LatLng.new(0, 45),
      Geokit::LatLng.new(0, 90),
      Geokit::LatLng.new(90, 0)
    ]
    @outer_north_triangle = [
      Geokit::LatLng.new(30, 0),
      Geokit::LatLng.new(30, 120),
      Geokit::LatLng.new(30, 240)
    ]
    @north_triangle_hole = [
      Geokit::LatLng.new(30, 0),
      Geokit::LatLng.new(30, 240),
      Geokit::LatLng.new(30, 120)
    ]
    @opposite_octant = [
      Geokit::LatLng.new(0, 180),
      Geokit::LatLng.new(-90, 0),
      Geokit::LatLng.new(0, -90)
    ]
    @octant_hole = [
      Geokit::LatLng.new(0, 20),
      Geokit::LatLng.new(30, 45),
      Geokit::LatLng.new(0, 70)
    ]
    @octant_hole_region = [
      Geokit::LatLng.new(0, 20),
      Geokit::LatLng.new(0, 70),
      Geokit::LatLng.new(30, 45)
    ]
  end

  def test_spherical_triangle_area_and_first_moment
    polygon = Geokit::SphericalPolygon.new(@octant, radius: 1)

    assert_in_delta Math::PI / 2, polygon.solid_angle, 1e-14
    assert_in_delta Math::PI / 2, polygon.area, 1e-14
    assert_vector_in_delta(
      [Math::PI / 4, Math::PI / 4, Math::PI / 4],
      polygon.first_moment
    )
  end

  def test_center_is_first_moment_direction
    center = Geokit::SphericalPolygon.new(@octant, radius: 1).center

    assert_in_delta Math.asin(1.0 / Math.sqrt(3)), center.lat * Math::PI / 180, 1e-14
    assert_in_delta 45, center.lng, 1e-14
  end

  def test_symmetric_region_center
    center = Geokit::SphericalPolygon.new(@outer_north_triangle, radius: 1).center

    assert_in_delta 90, center.lat, 1e-12
    assert_in_delta 0, center.lng, 1e-12
  end

  def test_reversed_ring_selects_larger_left_side_region
    polygon = Geokit::SphericalPolygon.new(@octant.reverse, radius: 1)

    assert_in_delta 7 * Math::PI / 2, polygon.solid_angle, 1e-14
    assert_vector_in_delta(
      [-Math::PI / 4, -Math::PI / 4, -Math::PI / 4],
      polygon.first_moment
    )
  end

  def test_hole_subtracts_signed_area_and_first_moment
    polygon = Geokit::SphericalPolygon.new(
      @outer_north_triangle,
      [@north_triangle_hole],
      radius: 1
    )

    assert_in_delta 0, polygon.solid_angle, 1e-14
    assert_vector_in_delta([0, 0, 0], polygon.first_moment)
  end

  def test_hole_inside_a_ring_subtracts_contained_region
    outer = Geokit::SphericalPolygon.new(@octant, radius: 1)
    contained = Geokit::SphericalPolygon.new(@octant_hole_region, radius: 1)
    polygon = Geokit::SphericalPolygon.new(@octant, [@octant_hole], radius: 1)

    assert_in_delta(
      outer.solid_angle - contained.solid_angle,
      polygon.solid_angle,
      1e-14
    )
    assert_vector_in_delta(
      outer.first_moment.zip(contained.first_moment).map { |left, right| left - right },
      polygon.first_moment
    )
  end

  def test_multiple_outer_rings_add
    polygon = Geokit::SphericalPolygon.new(
      [@octant, @opposite_octant],
      radius: 1
    )

    assert_in_delta Math::PI, polygon.solid_angle, 1e-14
    assert_vector_in_delta([0, 0, 0], polygon.first_moment)
    assert_nil polygon.center
  end

  def test_area_and_unnormalized_moment_add_after_split
    original = Geokit::SphericalPolygon.new(@octant, radius: 1)
    first = Geokit::SphericalPolygon.new(@triangle_one, radius: 1)
    second = Geokit::SphericalPolygon.new(@triangle_two, radius: 1)
    combined_moment = first.first_moment.zip(second.first_moment).map do |left, right|
      left + right
    end

    assert_in_delta original.solid_angle,
                    first.solid_angle + second.solid_angle,
                    1e-14
    assert_vector_in_delta(original.first_moment, combined_moment)
  end

  def test_subdividing_great_circle_edges_does_not_change_result
    polygon = Geokit::SphericalPolygon.new(@subdivided_octant, radius: 1)
    original = Geokit::SphericalPolygon.new(@octant, radius: 1)

    assert_in_delta original.solid_angle, polygon.solid_angle, 1e-14
    assert_vector_in_delta(original.first_moment, polygon.first_moment)
  end

  def test_radius_scales_area_and_moment
    polygon = Geokit::SphericalPolygon.new(@octant, radius: 2)

    assert_in_delta 2 * Math::PI, polygon.area, 1e-14
    assert_vector_in_delta(
      [2 * Math::PI, 2 * Math::PI, 2 * Math::PI],
      polygon.first_moment
    )
  end

  def test_polygon_wrapper
    polygon = Geokit::Polygon.new(@octant)

    assert_in_delta Math::PI / 2, polygon.spherical_area(radius: 1), 1e-14
    assert_vector_in_delta(
      [Math::PI / 4, Math::PI / 4, Math::PI / 4],
      polygon.spherical_first_moment(radius: 1)
    )
    assert_in_delta 45, polygon.spherical_center.lng, 1e-14
  end

  def test_antipodal_adjacent_points_are_rejected
    assert_raise ArgumentError do
      Geokit::SphericalPolygon.new(
        [
          Geokit::LatLng.new(0, 0),
          Geokit::LatLng.new(0, 180),
          Geokit::LatLng.new(90, 0)
        ],
        radius: 1
      )
    end
  end

  def test_coordinate_pair_ring_and_single_hole_are_accepted
    polygon = Geokit::SphericalPolygon.new(
      [[0, 0], [0, 90], [90, 0]],
      [[30, 0], [30, 240], [30, 120]],
      radius: 1
    )

    assert_in_delta(
      Math::PI / 2 - 2.000839033510993,
      polygon.solid_angle,
      1e-14
    )
  end

  def test_multiple_rings_of_latlng_are_accepted
    polygon = Geokit::SphericalPolygon.new([@octant, @opposite_octant])

    assert_equal 2, polygon.outer_rings.size

    polygon_with_holes = Geokit::SphericalPolygon.new(
      @octant,
      [@north_triangle_hole, @octant_hole]
    )
    assert_equal 2, polygon_with_holes.holes.size
  end

  def assert_vector_in_delta(expected, actual, delta = 1e-14)
    expected.zip(actual) do |expected_component, actual_component|
      assert_in_delta expected_component, actual_component, delta
    end
  end
end
