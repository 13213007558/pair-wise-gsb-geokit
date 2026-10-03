require File.join(File.dirname(__FILE__), 'helper')

class SphericalPolygonTest < Test::Unit::TestCase
  DEG = Math::PI / 180.0

  def setup
    @ll = lambda { |lat, lng| Geokit::LatLng.new(lat, lng) }
  end

  # A geodesic rectangle: edges are great-circle arcs between the corners.
  def rect(lat1, lat2, lng1, lng2)
    Geokit::SphericalPolygon.new([[lat1, lng1], [lat1, lng2], [lat2, lng2], [lat2, lng1]])
  end

  # Spherical wedge between two meridians, from the equator to the north
  # pole. All edges are meridians or the equator, so the area is exactly
  # (lng2 - lng1) in radians.
  def wedge(lng1, lng2)
    Geokit::SphericalPolygon.new([[0, lng1], [0, lng2], [90, 0]])
  end

  # Canonical representation of a polygon's rings, independent of vertex
  # order, start vertex and part order. Used to compare results for exact
  # geometric equality (up to rounding).
  def canonical(poly)
    poly.parts.map do |part|
      rings = [part.shell] + part.holes
      canon_rings = rings.map { |ring| canonical_ring(ring) }
      [canon_rings[0], canon_rings[1..-1].sort]
    end.sort_by(&:inspect)
  end

  def canonical_ring(ring)
    pts = ring.map { |p| [p.lat.round(7), p.lng.round(7)] }
    start = pts.each_index.min_by { |k| pts[k] }
    pts.rotate(start)
  end

  # --- point membership -------------------------------------------------

  def test_contains_across_antimeridian
    poly = rect(-10, 10, 170, -170)
    assert poly.contains?(@ll.call(0, 175))
    assert poly.contains?(@ll.call(0, -175))
    assert poly.contains?(@ll.call(0, 180))
    assert poly.contains?(@ll.call(-9, 179))
    assert !poly.contains?(@ll.call(0, 0))
    assert !poly.contains?(@ll.call(0, 90))
    assert !poly.contains?(@ll.call(0, 160))
    assert !poly.contains?(@ll.call(20, 180))
  end

  def test_contains_near_north_pole
    poly = rect(75, 85, 170, -170)
    assert poly.contains?(@ll.call(80, 180))
    assert poly.contains?(@ll.call(80, 175))
    assert poly.contains?(@ll.call(80, -175))
    assert !poly.contains?(@ll.call(80, 160))
    assert !poly.contains?(@ll.call(70, 180))
    assert !poly.contains?(@ll.call(86, 180))
    assert !poly.contains?(@ll.call(80, 0))
  end

  def test_contains_triangle_around_pole
    tri = Geokit::SphericalPolygon.new([[80, 0], [80, 120], [80, 240]])
    assert tri.contains?(@ll.call(89, 10))
    assert tri.contains?(@ll.call(89.9, -100))
    assert tri.contains?(@ll.call(86, 60))
    assert !tri.contains?(@ll.call(83, 60))
    assert !tri.contains?(@ll.call(70, 0))
  end

  def test_boundary_points_are_contained
    poly = rect(0, 10, 0, 10)
    assert poly.on_boundary?(@ll.call(0, 5))
    assert poly.on_boundary?(@ll.call(0, 0))
    assert poly.contains?(@ll.call(0, 5))
    assert poly.contains?(@ll.call(0, 0))
  end

  def test_equivalent_longitude_representations
    a = rect(-10, 10, 170, 190)
    b = rect(-10, 10, 170, -170)
    assert_in_delta a.area, b.area, 1e-12
    assert b.contains?(@ll.call(0, 180))
  end

  def test_wide_ring_spanning_most_longitudes
    # A band 240 degrees wide: membership and area must not depend on any
    # projection and must agree with the exact fan solid angle.
    band = Geokit::SphericalPolygon.new(
      [[-10, -120], [-10, 0], [-10, 120], [10, 120], [10, 0], [10, -120]]
    )
    assert_in_delta 2.3372859459660535, band.area, 1e-9
    assert band.contains?(@ll.call(0, -90))
    assert band.contains?(@ll.call(0, 90))
    assert band.contains?(@ll.call(0, 0))
    assert !band.contains?(@ll.call(0, 180))
    assert !band.contains?(@ll.call(20, 180))
    assert !band.contains?(@ll.call(-50, 100))

    cutter = rect(-20, 20, 90, 260)
    inter = band & cutter
    assert_in_delta band.area + cutter.area, (band | cutter).area + inter.area, 1e-12
    assert inter.contains?(@ll.call(0, 100))
    assert !inter.contains?(@ll.call(0, 0))
  end

  # --- exact areas with meridian wedges ---------------------------------

  def test_exact_wedge_areas
    a = wedge(170, -170) # 170..190, 20 degrees
    b = wedge(175, -175) # 175..185, 10 degrees, inside a
    assert_in_delta 20 * DEG, a.area, 1e-12
    assert_in_delta 10 * DEG, b.area, 1e-12
    assert_in_delta 10 * DEG, (a & b).area, 1e-12
    assert_in_delta 20 * DEG, (a | b).area, 1e-12
    assert_in_delta 10 * DEG, (a - b).area, 1e-12 # two 5-degree wedges
    assert_equal 2, (a - b).parts.size
  end

  def test_exact_shared_meridian_union
    a = wedge(170, 190)
    b = wedge(190, 200)
    union = a | b
    assert_in_delta 30 * DEG, union.area, 1e-12
    assert_equal 1, union.parts.size
    assert (a & b).empty?
  end

  # --- boolean algebra laws ----------------------------------------------

  def test_intersection_is_commutative
    a = Geokit::SphericalPolygon.new([[0, 170], [0, -170], [10, -160], [10, 175], [5, 172]])
    b = Geokit::SphericalPolygon.new([[2, 175], [2, -165], [12, -165], [12, 172]])
    assert_in_delta (a & b).area, (b & a).area, 1e-15
    assert_equal canonical(a & b), canonical(b & a)
  end

  def test_union_is_commutative
    a = rect(0, 10, 170, -170)
    b = rect(5, 15, 175, -165)
    assert_in_delta (a | b).area, (b | a).area, 1e-15
    assert_equal canonical(a | b), canonical(b | a)
  end

  def test_area_conservation
    a = rect(-10, 10, 170, -170)
    b = rect(-5, 15, 175, -165)
    assert_in_delta a.area + b.area, (a | b).area + (a & b).area, 1e-12
    assert_in_delta a.area, (a - b).area + (a & b).area, 1e-12
    assert_in_delta b.area, (b - a).area + (a & b).area, 1e-12
  end

  def test_area_conservation_near_pole
    a = rect(75, 85, 170, -170)
    b = rect(80, 88, 160, -160)
    assert_in_delta a.area + b.area, (a | b).area + (a & b).area, 1e-12
    assert_in_delta a.area, (a - b).area + (a & b).area, 1e-12
  end

  # --- degenerate contact between the two operands -----------------------

  def test_shared_partial_arc
    a = rect(0, 3, 0, 3)
    b = rect(1, 4, 3, 6) # shares the arc lat 1..3 on lng=3
    union = a | b
    assert_in_delta a.area + b.area, union.area, 1e-12
    assert_equal 1, union.parts.size
    assert (a & b).empty?
  end

  def test_shared_vertex_only
    a = rect(0, 2, 0, 2)
    b = rect(2, 4, 2, 4) # touches a at the single corner (2, 2)
    union = a | b
    assert_in_delta a.area + b.area, union.area, 1e-12
    assert (a & b).empty?
    assert union.contains?(@ll.call(1, 1))
    assert union.contains?(@ll.call(3, 3))
  end

  def test_tangent_collinear_edges
    a = wedge(170, 190)
    b = Geokit::SphericalPolygon.new([[-10, 190], [-10, 200], [0, 190]])
    union = a | b
    assert_in_delta a.area + b.area, union.area, 1e-12
    assert (a & b).empty?
  end

  def test_shared_edge_difference_keeps_boundary
    a = rect(0, 1, 0, 1)
    b = rect(0, 1, 1, 2) # adjacent, shares edge lng=1
    diff = a - b
    assert_in_delta a.area, diff.area, 1e-12
  end

  # --- holes and multiple shells ------------------------------------------

  def test_difference_carves_a_hole
    big = rect(0, 10, 0, 10)
    small = rect(3, 5, 3, 5)
    diff = big - small
    assert_equal 1, diff.parts.size
    assert_equal 1, diff.parts[0].holes.size
    assert_in_delta big.area - small.area, diff.area, 1e-12
    assert diff.contains?(@ll.call(1, 1))
    assert !diff.contains?(@ll.call(4, 4))
  end

  def test_intersection_with_multiple_shells
    donut = Geokit::SphericalPolygon.new([
      { shell: [[0, 0], [0, 10], [10, 10], [10, 0]],
        holes: [[[4, 4], [4, 6], [6, 6], [6, 4]]] }
    ])
    band = rect(4.5, 5.5, -1, 11)
    inter = donut & band
    assert_equal 2, inter.parts.size
    assert_in_delta donut.area + band.area, (donut | band).area + inter.area, 1e-12
    assert inter.contains?(@ll.call(5, 2))
    assert !inter.contains?(@ll.call(5, 5))
    assert inter.contains?(@ll.call(5, 8))
  end

  def test_union_can_create_a_hole
    u_shape = Geokit::SphericalPolygon.new(
      [[0, 0], [0, 10], [10, 10], [10, 0], [8, 0], [8, 8], [2, 8], [2, 0]]
    )
    plug = rect(1, 9, -2, 1)
    union = u_shape | plug
    assert_equal 1, union.parts.size
    assert_equal 1, union.parts[0].holes.size
    assert_in_delta u_shape.area + plug.area, union.area + (u_shape & plug).area, 1e-12
    assert !union.contains?(@ll.call(5, 5))
    assert union.contains?(@ll.call(5, 0.5))
    assert union.contains?(@ll.call(5, -1))
  end

  def test_boolean_results_can_be_chained
    u_shape = Geokit::SphericalPolygon.new(
      [[0, 0], [0, 10], [10, 10], [10, 0], [8, 0], [8, 8], [2, 8], [2, 0]]
    )
    union = u_shape | rect(1, 9, -2, 1)
    band = rect(4, 6, -3, 11)
    inter = union & band
    assert_equal 2, inter.parts.size
    assert_in_delta union.area + band.area, (union | band).area + inter.area, 1e-12
    assert inter.contains?(@ll.call(5, 0.5))
    assert !inter.contains?(@ll.call(5, 5))
  end

  def test_polygon_with_hole_contains
    donut = Geokit::SphericalPolygon.new([
      { shell: [[0, 0], [0, 10], [10, 10], [10, 0]],
        holes: [[[4, 4], [4, 6], [6, 6], [6, 4]]] }
    ])
    assert donut.contains?(@ll.call(2, 2))
    assert !donut.contains?(@ll.call(5, 5))
    assert !donut.contains?(@ll.call(20, 20))
    shell_only = Geokit::SphericalPolygon.new([[[0, 0], [0, 10], [10, 10], [10, 0]]])
    hole_only = Geokit::SphericalPolygon.new([[[4, 4], [4, 6], [6, 6], [6, 4]]])
    assert_in_delta shell_only.area - hole_only.area, donut.area, 1e-12
  end

  # --- input normalization -------------------------------------------------

  def test_orientation_is_normalized
    ccw = rect(0, 10, 0, 10)
    cw = Geokit::SphericalPolygon.new([[0, 0], [10, 0], [10, 10], [0, 10]])
    assert_in_delta ccw.area, cw.area, 1e-15
    assert cw.contains?(@ll.call(5, 5))
  end

  def test_start_vertex_rotation_does_not_change_geometry
    ring = [[0, 0], [0, 10], [10, 10], [10, 0]]
    other = rect(5, 15, 5, 15)
    base = Geokit::SphericalPolygon.new(ring)
    rotated = Geokit::SphericalPolygon.new(ring.rotate(2))
    assert_in_delta base.area, rotated.area, 1e-12
    assert_equal canonical(base & other), canonical(rotated & other)
    assert_equal canonical(base | other), canonical(rotated | other)
    assert_equal canonical(base - other), canonical(rotated - other)
  end

  def test_closing_vertex_is_dropped
    open = rect(0, 10, 0, 10)
    closed = Geokit::SphericalPolygon.new([[0, 0], [0, 10], [10, 10], [10, 0], [0, 0]])
    assert_in_delta open.area, closed.area, 1e-15
    assert_equal 4, closed.parts[0].shell.size
  end

  def test_empty_polygon
    empty = Geokit::SphericalPolygon.new([])
    poly = rect(0, 10, 0, 10)
    assert empty.empty?
    assert_equal 0.0, empty.area
    assert !empty.contains?(@ll.call(0, 0))
    assert (empty & poly).empty?
    assert (poly - poly).empty?
    assert_in_delta poly.area, (empty | poly).area, 1e-15
    assert_in_delta poly.area, (poly - empty).area, 1e-15
  end

  # --- invalid input is reported ------------------------------------------

  def test_self_intersecting_ring_is_rejected
    assert_raise(Geokit::SphericalPolygon::InvalidRingError) do
      Geokit::SphericalPolygon.new([[0, 0], [1, 1], [1, 0], [0, 1]])
    end
  end

  def test_antipodal_edge_is_rejected
    assert_raise(Geokit::SphericalPolygon::InvalidRingError) do
      Geokit::SphericalPolygon.new([[0, 0], [0, 180], [1, 1]])
    end
  end

  def test_half_circumference_edge_is_rejected
    assert_raise(Geokit::SphericalPolygon::InvalidRingError) do
      Geokit::SphericalPolygon.new([[10, 10], [-10, 190], [0, 45]])
    end
  end

  def test_too_few_vertices_is_rejected
    assert_raise(Geokit::SphericalPolygon::InvalidRingError) do
      Geokit::SphericalPolygon.new([[0, 0], [1, 1]])
    end
  end

  def test_hemisphere_ring_is_rejected
    # All four vertices lie on one great circle: the ring bounds exactly a
    # hemisphere and the interior is ambiguous.
    assert_raise(Geokit::SphericalPolygon::InvalidRingError) do
      Geokit::SphericalPolygon.new([[-20, 90], [-20, 270], [20, 270], [20, 90]])
    end
  end

  def test_hole_outside_shell_is_rejected
    assert_raise(Geokit::SphericalPolygon::InvalidRingError) do
      Geokit::SphericalPolygon.new([
        { shell: [[0, 0], [0, 10], [10, 10], [10, 0]],
          holes: [[[20, 20], [20, 21], [21, 21], [21, 20]]] }
      ])
    end
  end

  def test_crossing_rings_are_rejected
    assert_raise(Geokit::SphericalPolygon::InvalidRingError) do
      Geokit::SphericalPolygon.new([
        { shell: [[0, 0], [0, 10], [10, 10], [10, 0]],
          holes: [[[5, 5], [5, 15], [15, 15], [15, 5]]] }
      ])
    end
  end

  def test_overlapping_shells_are_rejected
    assert_raise(Geokit::SphericalPolygon::InvalidRingError) do
      Geokit::SphericalPolygon.new([
        [[0, 0], [0, 10], [10, 10], [10, 0]],
        [[2, 2], [2, 3], [3, 3], [3, 2]]
      ])
    end
  end
end
