# encoding: utf-8
require File.join(File.dirname(__FILE__), 'helper')

class SphericalDelaunayTest < Test::Unit::TestCase #:nodoc: all
  TOLERANCE = 1e-9
  FOUR_PI = 4.0 * Math::PI

  def setup
    @tetrahedron = [
      [1.0, 1.0, 1.0],
      [1.0, -1.0, -1.0],
      [-1.0, 1.0, -1.0],
      [-1.0, -1.0, 1.0]
    ]
    @octahedron = [
      [1.0, 0.0, 0.0], [-1.0, 0.0, 0.0],
      [0.0, 1.0, 0.0], [0.0, -1.0, 0.0],
      [0.0, 0.0, 1.0], [0.0, 0.0, -1.0]
    ]
    phi = (1.0 + Math.sqrt(5.0)) / 2.0
    @icosahedron = [
      [0.0, 1.0, phi], [0.0, 1.0, -phi], [0.0, -1.0, phi], [0.0, -1.0, -phi],
      [1.0, phi, 0.0], [1.0, -phi, 0.0], [-1.0, phi, 0.0], [-1.0, -phi, 0.0],
      [phi, 0.0, 1.0], [-phi, 0.0, 1.0], [phi, 0.0, -1.0], [-phi, 0.0, -1.0]
    ]
  end

  # -- regular symmetric point sets -----------------------------------------

  def test_tetrahedron
    triangulation = Geokit::SphericalDelaunay.triangulate(@tetrahedron)
    assert_equal 4, triangulation.vertex_count
    assert_equal 4, triangulation.face_count
    assert_equal 6, triangulation.edge_count
    assert_valid_sphere_triangulation(triangulation)
    # all faces of a regular tetrahedron share the same circumradius
    radii = triangulation.faces.map(&:circumradius).uniq
    assert_equal 1, radii.size
    assert_empty triangulation.degeneracies
  end

  def test_octahedron
    triangulation = Geokit::SphericalDelaunay.triangulate(@octahedron)
    assert_equal 6, triangulation.vertex_count
    assert_equal 8, triangulation.face_count
    assert_equal 12, triangulation.edge_count
    assert_valid_sphere_triangulation(triangulation)
    assert_empty triangulation.degeneracies
  end

  def test_icosahedron
    triangulation = Geokit::SphericalDelaunay.triangulate(@icosahedron)
    assert_equal 12, triangulation.vertex_count
    assert_equal 20, triangulation.face_count
    assert_equal 30, triangulation.edge_count
    assert_valid_sphere_triangulation(triangulation)
    assert_empty triangulation.degeneracies
  end

  # -- general point sets ----------------------------------------------------

  def test_random_points_empty_circumcap
    points = random_unit_vectors(40, 1234)
    triangulation = Geokit::SphericalDelaunay.triangulate(points)
    assert_valid_sphere_triangulation(triangulation)
    assert_empty_circumcaps(triangulation)
    assert_empty triangulation.degeneracies
  end

  def test_random_points_larger_set
    points = random_unit_vectors(80, 987_654)
    triangulation = Geokit::SphericalDelaunay.triangulate(points)
    assert_valid_sphere_triangulation(triangulation)
    assert_empty_circumcaps(triangulation)
  end

  def test_lat_lng_input
    stations = [
      Geokit::LatLng.new(35.264389682754654, 45.0),
      Geokit::LatLng.new(-35.264389682754654, -45.0),
      Geokit::LatLng.new(-35.264389682754654, 135.0),
      Geokit::LatLng.new(35.264389682754654, -135.0)
    ]
    triangulation = Geokit::SphericalDelaunay.triangulate(stations)
    assert_equal 4, triangulation.face_count
    assert_valid_sphere_triangulation(triangulation)
    assert_equal((0..3).to_a, triangulation.index_map.keys.sort)
  end

  def test_array_lat_lng_input
    stations = [
      [35.264389682754654, 45.0],
      [-35.264389682754654, -45.0],
      [-35.264389682754654, 135.0],
      [35.264389682754654, -135.0]
    ]
    triangulation = Geokit::SphericalDelaunay.triangulate(stations)
    assert_equal 4, triangulation.face_count
    assert_valid_sphere_triangulation(triangulation)
  end

  # -- identity mapping and duplicates ---------------------------------------

  def test_duplicate_stations_are_merged_and_mapped
    points = @octahedron.dup
    points << [1.0, 0.0, 0.0]          # exact duplicate of vertex 0
    points << [1.0, 1e-12, 0.0]        # near duplicate within tolerance
    triangulation = Geokit::SphericalDelaunay.triangulate(points)
    assert_equal 6, triangulation.vertex_count
    assert_equal 8, triangulation.face_count
    assert_equal 8, triangulation.index_map.size
    # all three copies of the +x station map to the same canonical station
    assert_equal triangulation.index_map[0], triangulation.index_map[6]
    assert_equal triangulation.index_map[0], triangulation.index_map[7]
    duplicates = triangulation.degeneracies.select { |d| d.kind == :duplicate_stations }
    assert_equal 1, duplicates.size
    assert_equal [0, 6, 7], duplicates.first.data
    assert_valid_sphere_triangulation(triangulation)
  end

  # -- order independence ----------------------------------------------------

  def test_input_order_does_not_change_topology
    points = random_unit_vectors(30, 42)
    reference = topology_signature(Geokit::SphericalDelaunay.triangulate(points))
    5.times do |seed|
      shuffled = points.shuffle(random: Random.new(seed + 100))
      signature = topology_signature(Geokit::SphericalDelaunay.triangulate(shuffled))
      assert_equal reference, signature, "topology changed under permutation #{seed}"
    end
  end

  def test_cocircular_stations_resolve_deterministically
    points = cocircular_set
    reference = topology_signature(Geokit::SphericalDelaunay.triangulate(points))
    cocircular = Geokit::SphericalDelaunay.triangulate(points).degeneracies.select do |d|
      d.kind == :cocircular_stations
    end
    assert_equal 1, cocircular.size
    5.times do |seed|
      shuffled = points.shuffle(random: Random.new(seed + 7))
      signature = topology_signature(Geokit::SphericalDelaunay.triangulate(shuffled))
      assert_equal reference, signature, "cocircular topology changed under permutation #{seed}"
    end
  end

  # -- antimeridian ----------------------------------------------------------

  def test_no_artificial_cut_at_180_meridian
    stations = [
      [10.0, 179.0],   # just east of the antimeridian
      [10.0, -179.0],  # just west of the antimeridian: 2 degrees away, not 358
      [10.0, 0.0],
      [-40.0, 120.0],
      [-40.0, -120.0],
      [60.0, 90.0]
    ]
    triangulation = Geokit::SphericalDelaunay.triangulate(stations)
    assert_valid_sphere_triangulation(triangulation)
    east = triangulation.index_map[0]
    west = triangulation.index_map[1]
    connected = triangulation.edges.any? do |edge|
      edge.vertices == [east, west].sort
    end
    assert connected, 'nearest neighbours across the 180 degree meridian must be connected'
  end

  # -- hemisphere coverage precondition --------------------------------------

  def test_closed_hemisphere_is_not_applicable
    northern = [
      [0.0, 0.0], [0.0, 90.0], [0.0, 180.0], [0.0, -90.0],
      [45.0, 45.0], [45.0, -135.0]
    ]
    assert_raise(Geokit::SphericalDelaunay::NotApplicableError) do
      Geokit::SphericalDelaunay.triangulate(northern)
    end
  end

  def test_open_hemisphere_is_not_applicable
    tropical = [
      [10.0, 0.0], [10.0, 90.0], [10.0, 180.0], [10.0, -90.0],
      [30.0, 45.0], [-5.0, -45.0]
    ]
    assert_raise(Geokit::SphericalDelaunay::NotApplicableError) do
      Geokit::SphericalDelaunay.triangulate(tropical)
    end
  end

  def test_too_few_stations_is_not_applicable
    assert_raise(Geokit::SphericalDelaunay::NotApplicableError) do
      Geokit::SphericalDelaunay.triangulate([[0.0, 0.0], [0.0, 120.0], [0.0, -120.0]])
    end
  end

  # -- degeneracy reporting --------------------------------------------------

  def test_near_cocircular_data_is_reported
    points = cocircular_set
    # nudge one cocircular station by ~5e-9 rad: outside the decision
    # tolerance but inside the unreliability band
    points[1] = [30.0 + 3e-7, 90.0]
    triangulation = Geokit::SphericalDelaunay.triangulate(points)
    assert_valid_sphere_triangulation(triangulation)
    kinds = triangulation.degeneracies.map(&:kind)
    assert kinds.include?(:near_degenerate),
      "expected a reported degeneracy, got #{kinds.inspect}"
  end

  def test_invalid_input_raises_argument_error
    assert_raise(ArgumentError) do
      Geokit::SphericalDelaunay.triangulate([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    end
  end

  private

  def random_unit_vectors(count, seed)
    random = Random.new(seed)
    Array.new(count) do
      z = random.rand * 2.0 - 1.0
      theta = random.rand * 2.0 * Math::PI
      r = Math.sqrt(1.0 - z * z)
      [r * Math.cos(theta), r * Math.sin(theta), z]
    end
  end

  # 4 stations exactly on the lat=30 circle plus 4 more below it, so the
  # lat=30 plane is a supporting plane with 4 cocircular stations.
  def cocircular_set
    [
      [30.0, 0.0], [30.0, 90.0], [30.0, 180.0], [30.0, -90.0],
      [0.0, 45.0], [-30.0, 0.0], [-30.0, 120.0], [-30.0, -120.0]
    ]
  end

  def topology_signature(triangulation)
    triangulation.faces.map { |face| face.vertices.sort }.sort
  end

  def assert_valid_sphere_triangulation(triangulation)
    assert_equal 2, triangulation.euler_characteristic,
      "Euler characteristic V-E+F must be 2 (V=#{triangulation.vertex_count}, E=#{triangulation.edge_count}, F=#{triangulation.face_count})"
    assert_in_delta FOUR_PI, triangulation.total_area, 1e-6,
      'faces must cover the whole sphere'
    # every edge has exactly two oppositely directed half-edges
    triangulation.edges.each do |edge|
      assert_equal 2, edge.faces.size, "edge #{edge.vertices.inspect} must border exactly 2 faces"
      from_a, to_a, = edge.half_edges[0]
      from_b, to_b, = edge.half_edges[1]
      assert_equal [from_a, to_a], [to_b, from_b],
        "edge #{edge.vertices.inspect} must be traversed in opposite directions by its two faces"
    end
  end

  def assert_empty_circumcaps(triangulation)
    triangulation.faces.each do |face|
      cos_radius = Math.cos(face.circumradius)
      triangulation.stations.each_with_index do |station, index|
        next if face.vertices.include?(index)

        dot = Geokit::SphericalDelaunay.dot(station, face.circumcenter)
        assert dot <= cos_radius + TOLERANCE,
          "station #{index} lies inside the circumcap of face #{face.vertices.inspect}"
      end
    end
  end
end
