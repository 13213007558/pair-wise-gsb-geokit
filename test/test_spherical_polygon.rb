# encoding: utf-8
require File.join(File.dirname(__FILE__), 'helper')

class SphericalPolygonTest < Test::Unit::TestCase #:nodoc: all
  Spherical = Geokit::Spherical

  # Builds a polygon from [lat, lng] arrays: first ring is the outer ring,
  # the rest are holes.
  def poly(*rings)
    Spherical::Polygon.new(rings[0], rings[1..] || [])
  end

  # Square-ish region, counterclockwise (interior on the left).
  def box(lat0, lng0, lat1, lng1)
    poly([[lat0, lng0], [lat0, lng1], [lat1, lng1], [lat1, lng0]])
  end

  # Canonical, start-rotation-invariant signature of a MultiPolygon.
  def signature(multi)
    multi.polygons.map do |polygon|
      rings = ([polygon.outer] + polygon.holes).map do |ring|
        pts = ring.map { |ll| [ll.lat.round(9), ll.lng.round(9)] }
        start = pts.each_with_index.min_by { |pt, _i| pt }.last
        pts.rotate(start)
      end
      [rings[0], (rings[1..] || []).sort]
    end.sort_by { |outer, _holes| outer }
  end

  def antimeridian_box
    # Lat band -10..10, lng band 170 -> 180 -> -170.
    poly([[-10, 170], [-10, -170], [10, -170], [10, 170]])
  end

  # --- point membership -------------------------------------------------

  def test_contains_across_antimeridian
    region = antimeridian_box
    assert region.contains?([0, 180])
    assert region.contains?([0, 179.9])
    assert region.contains?([0, -179.9])
    assert region.contains?([5, 175])
    assert region.contains?([-5, -175])
    assert !region.contains?([0, 150])
    assert !region.contains?([0, -150])
    assert !region.contains?([20, 180])
    assert !region.contains?([-20, 180])
  end

  def test_contains_near_north_pole
    cap = poly([[80, 0], [80, 120], [80, -120]])
    assert cap.contains?([90, 0])   # the pole itself
    assert cap.contains?([89, 45])
    assert cap.contains?([85, 200])
    assert !cap.contains?([70, 0])
    assert !cap.contains?([0, 0])
    assert !cap.contains?([-80, 0])
  end

  def test_contains_is_start_rotation_invariant
    ring = [[-10, 170], [-10, -170], [10, -170], [10, 170]]
    rotated = ring.rotate(2)
    assert_equal poly(ring).contains?([0, 180]), poly(rotated).contains?([0, 180])
    assert_equal poly(ring).contains?([0, 150]), poly(rotated).contains?([0, 150])
  end

  def test_closed_ring_input_is_accepted
    ring = [[-10, 170], [-10, -170], [10, -170], [10, 170], [-10, 170]]
    assert poly(ring).contains?([0, 180])
  end

  # --- boolean operations ------------------------------------------------

  def test_intersection_is_commutative
    a = antimeridian_box
    b = poly([[-5, 160], [-5, -160], [5, -160], [5, 160]])
    ab = a.intersection(b)
    ba = b.intersection(a)
    assert_in_delta ab.area, ba.area, 1e-12
    assert_equal signature(ab), signature(ba)
    assert ab.contains?([0, 175])
    assert !ab.contains?([8, 175])
  end

  def test_area_conservation
    a = antimeridian_box
    b = poly([[-5, 160], [-5, -160], [5, -160], [5, 160]])
    total = a.area + b.area
    assert_in_delta total, a.union(b).area + a.intersection(b).area, 1e-9
    assert_in_delta a.area, a.difference(b).area + a.intersection(b).area, 1e-9
    assert_in_delta b.area, b.difference(a).area + b.intersection(a).area, 1e-9
  end

  def test_area_conservation_near_pole
    a = poly([[80, 0], [80, 120], [80, -120]])
    b = poly([[75, 60], [75, 180], [75, -60]])
    total = a.area + b.area
    assert_in_delta total, a.union(b).area + a.intersection(b).area, 1e-9
    assert_in_delta a.area, a.difference(b).area + a.intersection(b).area, 1e-9
  end

  def test_boolean_result_is_start_rotation_invariant
    ring_a = [[-10, 170], [-10, -170], [10, -170], [10, 170]]
    b = poly([[-5, 160], [-5, -160], [5, -160], [5, 160]])
    base = poly(ring_a)
    rotated = poly(ring_a.rotate(2))
    assert_equal signature(base.intersection(b)), signature(rotated.intersection(b))
    assert_equal signature(base.union(b)), signature(rotated.union(b))
    assert_equal signature(base.difference(b)), signature(rotated.difference(b))
  end

  def test_disjoint_union_and_empty_intersection
    a = box(0, 0, 2, 2)
    b = box(0, 10, 2, 12)
    union = a.union(b)
    assert_equal 2, union.polygons.size
    assert_in_delta a.area + b.area, union.area, 1e-12
    assert a.intersection(b).empty?
    assert_equal 0.0, a.intersection(b).area
  end

  def test_shared_edge_union_merges
    a = box(0, 0, 2, 2)
    b = box(0, 2, 2, 4)
    union = a.union(b)
    assert_equal 1, union.polygons.size
    assert_equal 1, union.polygons[0].outer.size >= 4 ? 1 : 0 # sanity
    assert_in_delta a.area + b.area, union.area, 1e-9
    assert union.contains?([1, 3])
    assert union.contains?([1, 1])
    # The shared edge (the lng=2 meridian arc) must be gone from the
    # boundary: nothing outside a or b is covered.
    assert !union.contains?([1, 5])
    assert !union.contains?([3, 1])
  end

  def test_shared_edge_intersection_and_difference
    a = box(0, 0, 2, 2)
    b = box(0, 2, 2, 4)
    assert a.intersection(b).empty?
    difference = a.difference(b)
    assert_in_delta a.area, difference.area, 1e-9
    assert difference.contains?([1, 1])
    assert !difference.contains?([1, 3])
  end

  def test_point_touching_union
    a = box(0, 0, 2, 2)
    b = box(2, 2, 4, 4)
    union = a.union(b)
    assert_in_delta a.area + b.area, union.area, 1e-9
    assert union.contains?([3, 3])
    assert union.contains?([1, 1])
    assert !union.contains?([1, 3])
  end

  def test_difference_creates_hole
    a = box(0, 0, 10, 10)
    b = box(4, 4, 6, 6)
    difference = a.difference(b)
    assert_equal 1, difference.polygons.size
    assert_equal 1, difference.polygons[0].holes.size
    assert_in_delta a.area - b.area, difference.area, 1e-9
    assert difference.contains?([1, 1])
    assert !difference.contains?([5, 5])
    # Intersection with the hole region is exactly b.
    assert_in_delta b.area, a.intersection(b).area, 1e-9
    # Union restores a.
    assert_in_delta a.area, a.union(b).area, 1e-9
  end

  def test_difference_with_shared_boundary_arc
    # b sits inside a sharing the arc lat=0, lng in 4..6 (same direction).
    a = box(0, 0, 10, 10)
    b = box(0, 4, 6, 6)
    difference = a.difference(b)
    assert_in_delta a.area - b.area, difference.area, 1e-9
    assert !difference.contains?([3, 5])
    assert difference.contains?([3, 7])
    assert_in_delta a.area, a.union(b).area, 1e-9
    assert_in_delta b.area, a.intersection(b).area, 1e-9
  end

  def test_antimeridian_union_merges_across_180
    a = poly([[-10, 170], [-10, -170], [10, -170], [10, 170]])
    b = poly([[-10, -170], [-10, -150], [10, -150], [10, -170]])
    union = a.union(b)
    assert_equal 1, union.polygons.size
    assert_in_delta a.area + b.area, union.area, 1e-9
    assert union.contains?([0, -160])
    assert union.contains?([0, 180])
  end

  def test_multipolygon_chaining
    a = box(0, 0, 2, 2)
    b = box(1, 1, 3, 3)
    c = box(1.5, 1.5, 2.5, 2.5)
    union = a.union(b)
    assert_in_delta (a.area + b.area - a.intersection(b).area), union.area, 1e-9
    both = union.intersection(c)
    # c is fully inside b, so (a | b) & c == c.
    assert_in_delta c.area, both.area, 1e-9
  end

  def test_t_junction_vertex_on_edge
    # A corner of b sits exactly on the right edge of a.
    a = box(0, 0, 4, 4)
    b = box(2, 4, 6, 6)
    union = a.union(b)
    assert_equal 1, union.polygons.size
    assert_in_delta a.area + b.area, union.area, 1e-9
    assert union.contains?([4, 5])
    assert union.contains?([1, 1])
    assert !union.contains?([5, 1])
    assert a.intersection(b).empty?
  end

  def test_partially_shared_arc
    # Bottom edges of a and b overlap over lng 2..4 (same direction).
    a = box(0, 0, 2, 4)
    b = box(0, 2, 2, 6)
    union = a.union(b)
    assert_equal 1, union.polygons.size
    assert_in_delta a.area + b.area - a.intersection(b).area, union.area, 1e-9
    intersection = a.intersection(b)
    assert intersection.contains?([1, 3])
    assert !intersection.contains?([1, 5])
    assert !intersection.contains?([3, 3])
    difference = a.difference(b)
    assert_in_delta a.area - intersection.area, difference.area, 1e-9
    assert difference.contains?([1, 1])
    assert !difference.contains?([1, 3])
  end

  def test_boolean_with_holed_input
    a = poly(
      [[0, 0], [0, 10], [10, 10], [10, 0]],
      [[4, 4], [4, 6], [6, 6], [6, 4]]
    )
    # Horizontal strip covering the hole entirely.
    b = box(3, -1, 7, 11)
    intersection = a.intersection(b)
    assert_equal 1, intersection.polygons.size
    assert_equal 1, intersection.polygons[0].holes.size
    assert_in_delta a.area + b.area, a.union(b).area + intersection.area, 1e-9
    assert_in_delta a.area, a.difference(b).area + intersection.area, 1e-9
    assert intersection.contains?([5, 1])
    assert !intersection.contains?([5, 5])
  end

  def test_module_level_functions
    a = box(0, 0, 2, 2)
    b = box(1, 1, 3, 3)
    assert_in_delta a.union(b).area, Spherical.union(a, b).area, 1e-12
    assert_in_delta a.intersection(b).area, Spherical.intersection(a, b).area, 1e-12
    assert_in_delta a.difference(b).area, Spherical.difference(a, b).area, 1e-12
  end

  def test_chained_union_forms_hole
    # Four boxes sharing meridian edges, surrounding a central void.
    left = box(0, 0, 10, 4)
    right = box(0, 6, 10, 10)
    top = box(6, 4, 10, 6)
    bottom = box(0, 4, 4, 6)
    union = left.union(right).union(bottom).union(top)
    assert_equal 1, union.polygons.size
    assert_equal 1, union.polygons[0].holes.size
    total = left.area + right.area + top.area + bottom.area
    assert_in_delta total, union.area, 1e-9
    assert !union.contains?([5, 5])   # the void
    assert union.contains?([5, 2])
    assert union.contains?([2, 5])
    assert union.contains?([8, 5])
  end

  # --- invalid input reporting -------------------------------------------

  def test_self_intersecting_ring_is_reported
    assert_raise(Spherical::SelfIntersectionError) do
      poly([[0, 0], [2, 2], [2, 0], [0, 2]])
    end
  end

  def test_repeated_vertex_ring_is_reported
    assert_raise(Spherical::SelfIntersectionError) do
      poly([[0, 0], [0, 2], [2, 2], [0, 0], [0, -2], [2, -2]])
    end
  end

  def test_antipodal_edge_is_reported
    assert_raise(Spherical::AntipodalEdgeError) do
      poly([[0, 0], [0, 180], [45, 90]])
    end
  end

  def test_too_few_vertices_is_reported
    assert_raise(Spherical::InvalidPolygonError) do
      poly([[0, 0], [0, 1]])
    end
  end

  def test_hole_outside_outer_is_reported
    assert_raise(Spherical::RingNestingError) do
      poly([[0, 0], [0, 10], [10, 10], [10, 0]],
        [[20, 20], [20, 21], [21, 21], [21, 20]])
    end
  end

  def test_hole_inside_hole_is_reported
    assert_raise(Spherical::RingNestingError) do
      poly([[0, 0], [0, 10], [10, 10], [10, 0]],
        [[1, 1], [1, 9], [9, 9], [9, 1]],
        [[2, 2], [2, 3], [3, 3], [3, 2]])
    end
  end

  def test_polygon_with_hole_area_and_contains
    region = poly(
      [[0, 0], [0, 10], [10, 10], [10, 0]],
      [[4, 4], [6, 4], [6, 6], [4, 6]] # clockwise hole, auto-oriented anyway
    )
    outer_only = box(0, 0, 10, 10)
    hole_only = box(4, 4, 6, 6)
    assert_in_delta outer_only.area - hole_only.area, region.area, 1e-9
    assert region.contains?([1, 1])
    assert !region.contains?([5, 5])
  end
end
