# encoding: utf-8
require File.join(File.dirname(__FILE__), 'helper')

class SphericalPolygonTest < Test::Unit::TestCase #:nodoc: all
  def setup
    @a = Geokit::LatLng.new(0, 0)
    @b = Geokit::LatLng.new(0, 90)
    @c = Geokit::LatLng.new(90, 0)
    @d = Geokit::LatLng.new(45, 90)

    @triangle = Geokit::SphericalPolygon.new([@a, @b, @c])
  end

  def test_analytic_spherical_triangle_area_and_symmetric_center
    assert_in_delta Math::PI / 2.0, @triangle.solid_angle, 1e-12
    assert_in_delta Math::PI / 2.0, @triangle.area, 1e-12

    expected_moment = Math::PI / 4.0
    assert_in_delta expected_moment, @triangle.first_moment.x, 1e-12
    assert_in_delta expected_moment, @triangle.first_moment.y, 1e-12
    assert_in_delta expected_moment, @triangle.first_moment.z, 1e-12

    assert_in_delta Math.asin(1.0 / Math.sqrt(3.0)) * 180.0 / Math::PI,
                    @triangle.center.lat, 1e-12
    assert_in_delta 45.0, @triangle.center.lng, 1e-12
  end

  def test_small_geodesic_triangle_has_positive_area
    small = Geokit::SphericalPolygon.new([
      Geokit::LatLng.new(10.0, 20.0),
      Geokit::LatLng.new(10.0, 20.001),
      Geokit::LatLng.new(10.001, 20.0005)
    ])

    assert_operator small.solid_angle, :>, 0.0
    assert_in_delta 10.0003, small.center.lat, 1e-4
    assert_in_delta 20.0005, small.center.lng, 1e-4
  end

  def test_center_aliases_spherical_centroid
    assert_equal @triangle.center, @triangle.spherical_centroid
  end

  def test_area_and_moment_add_when_region_is_split
    first = Geokit::SphericalPolygon.new([@a, @b, @d])
    second = Geokit::SphericalPolygon.new([@a, @d, @c])
    combined = Geokit::SphericalPolygon.new([[@a, @b, @d], [@a, @d, @c]])

    assert_in_delta @triangle.solid_angle,
                    first.solid_angle + second.solid_angle, 1e-12
    assert_in_delta @triangle.solid_angle, combined.solid_angle, 1e-12

    summed_moment = first.first_moment + second.first_moment
    @triangle.first_moment.to_a.each_with_index do |expected, index|
      assert_in_delta expected, summed_moment.to_a[index], 1e-12
      assert_in_delta expected, combined.first_moment.to_a[index], 1e-12
    end
  end

  def test_subdividing_a_great_circle_edge_does_not_change_geometry
    subdivided = Geokit::SphericalPolygon.new([
      @a,
      Geokit::LatLng.new(0, 22.5),
      Geokit::LatLng.new(0, 45.0),
      Geokit::LatLng.new(0, 67.5),
      @b,
      @c
    ])

    assert_in_delta @triangle.solid_angle, subdivided.solid_angle, 1e-12
    @triangle.first_moment.to_a.each_with_index do |expected, index|
      assert_in_delta expected, subdivided.first_moment.to_a[index], 1e-12
    end
  end

  def test_multiple_outer_rings_are_added
    opposite = Geokit::SphericalPolygon.new([
      [@a, @b, @c],
      [
        Geokit::LatLng.new(0, 180),
        Geokit::LatLng.new(-90, 0),
        Geokit::LatLng.new(0, -90)
      ]
    ])

    assert_in_delta Math::PI, opposite.solid_angle, 1e-12
    opposite.first_moment.to_a.each do |component|
      assert_in_delta 0.0, component, 1e-12
    end
    assert_nil opposite.center
  end

  def test_holes_subtract_area_and_moment
    outer = [0, 120, 240].map { |longitude| Geokit::LatLng.new(30, longitude) }
    hole = [240, 120, 0].map { |longitude| Geokit::LatLng.new(60, longitude) }
    outer_polygon = Geokit::SphericalPolygon.new(outer)
    hole_polygon = Geokit::SphericalPolygon.new(hole.reverse)
    with_hole = Geokit::SphericalPolygon.new(outer, holes: hole)

    assert_in_delta outer_polygon.solid_angle - hole_polygon.solid_angle,
                    with_hole.solid_angle, 1e-12
    expected_moment = outer_polygon.first_moment - hole_polygon.first_moment
    expected_moment.to_a.each_with_index do |expected, index|
      assert_in_delta expected, with_hole.first_moment.to_a[index], 1e-12
    end
  end

  def test_reversed_ring_selects_the_side_larger_than_a_hemisphere
    complement = Geokit::SphericalPolygon.new([@a, @c, @b])

    assert_in_delta 4.0 * Math::PI - Math::PI / 2.0,
                    complement.solid_angle, 1e-12
    @triangle.first_moment.to_a.each_with_index do |component, index|
      assert_in_delta(-component, complement.first_moment.to_a[index], 1e-12)
    end
    assert_in_delta(-35.264389682754654, complement.center.lat, 1e-12)
    assert_in_delta(-135.0, complement.center.lng, 1e-12)
  end

  def test_hemisphere_has_a_unique_center
    hemisphere = Geokit::SphericalPolygon.new([
      Geokit::LatLng.new(0, 0),
      Geokit::LatLng.new(0, 90),
      Geokit::LatLng.new(0, 180),
      Geokit::LatLng.new(0, 270)
    ])

    assert_in_delta 2.0 * Math::PI, hemisphere.solid_angle, 1e-12
    assert_in_delta 0.0, hemisphere.first_moment.x, 1e-12
    assert_in_delta 0.0, hemisphere.first_moment.y, 1e-12
    assert_in_delta Math::PI, hemisphere.first_moment.z, 1e-12
    assert_in_delta 90.0, hemisphere.center.lat, 1e-12
  end

  def test_whole_sphere_has_zero_first_moment_and_no_unique_center
    whole = Geokit::SphericalPolygon.whole_sphere

    assert_in_delta 4.0 * Math::PI, whole.solid_angle, 1e-12
    whole.first_moment.to_a.each { |component| assert_in_delta 0.0, component, 1e-12 }
    assert_nil whole.center
  end

  def test_closed_rings_and_repeated_edge_points_are_supported
    closed = Geokit::SphericalPolygon.new([@a, @b, @c, @a])
    repeated = Geokit::SphericalPolygon.new([@a, @a, @b, @b, @c])

    assert_in_delta @triangle.solid_angle, closed.solid_angle, 1e-12
    assert_in_delta @triangle.solid_angle, repeated.solid_angle, 1e-12
  end

  def test_antipodal_adjacent_vertices_are_rejected
    error = assert_raise ArgumentError do
      Geokit::SphericalPolygon.new([
        Geokit::LatLng.new(0, 0),
        Geokit::LatLng.new(0, 180),
        Geokit::LatLng.new(90, 0)
      ])
    end
    assert_match(/unique great-circle arc/, error.message)
  end

  def test_radius_scales_area_and_moment
    radius = 10.0

    assert_in_delta @triangle.solid_angle * radius * radius,
                    @triangle.area(radius: radius), 1e-12
    expected_moment = Math::PI / 4.0 * radius**3
    @triangle.first_moment(radius: radius).to_a.each do |component|
      assert_in_delta expected_moment, component, 1e-9
    end
  end

  def test_units_select_a_named_sphere_radius
    radius = Geokit::LatLng.units_sphere_multiplier(:meters)

    assert_in_delta @triangle.area(radius: radius),
                    @triangle.area(units: :meters), 1e-3
  end
end
