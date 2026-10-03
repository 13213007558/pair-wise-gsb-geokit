require File.join(File.dirname(__FILE__), 'helper')

class SphericalDelaunayTest < Test::Unit::TestCase #:nodoc: all
  SD = Geokit::SphericalDelaunay
  FOUR_PI = 4.0 * Math::PI

  def regular_tetrahedron
    [[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]].map do |point|
      normalize(point)
    end
  end

  def octahedron
    [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]]
  end

  def cube
    signs = [1, -1]
    signs.product(signs, signs).map { |point| normalize(point) }
  end

  def icosahedron
    phi = (1.0 + Math.sqrt(5.0)) / 2.0
    vertices = []
    [1, -1].each do |a|
      [phi, -phi].each do |b|
        vertices << [0, a, b] << [a, b, 0] << [b, 0, a]
      end
    end
    vertices.map { |point| normalize(point) }
  end

  def random_sphere(count, seed)
    random = Random.new(seed)
    Array.new(count) do
      z = random.rand * 2.0 - 1.0
      angle = random.rand * 2.0 * Math::PI
      radius = Math.sqrt(1.0 - z * z)
      [radius * Math.cos(angle), radius * Math.sin(angle), z]
    end
  end

  def normalize(point)
    length = Math.sqrt(point.sum { |value| value * value })
    point.map { |value| value / length }
  end

  def face_normal(result, face)
    a = result.vertices[face[0]]
    b = result.vertices[face[1]]
    c = result.vertices[face[2]]
    [
      (b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1]),
      (b[2] - a[2]) * (c[0] - a[0]) - (b[0] - a[0]) * (c[2] - a[2]),
      (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
    ]
  end

  def unit_normal(result, face)
    normal = face_normal(result, face)
    length = Math.sqrt(normal.sum { |value| value * value })
    normal.map { |value| value / length }
  end

  # Shared geometric contract asserted by every positive test case.
  def assert_spherical_complex(result, point_count)
    assert_equal(point_count, result.vertex_count)
    assert_equal(2, result.euler_characteristic)
    assert_in_delta(FOUR_PI, result.total_area, 1e-7)

    incidence = Hash.new(0)
    result.faces.each do |face|
      normal = unit_normal(result, face)
      anchor = result.vertices[face[0]]
      # Outward orientation: the facet normal points away from the origin.
      assert_operator(normal.zip(anchor).sum { |x, y| x * y }, :>, 0.0)
      [[0, 1], [1, 2], [2, 0]].each do |left, right|
        incidence[[face[left], face[right]].sort] += 1
      end
    end

    # Every undirected edge has exactly one face on each side.
    assert(result.edges.all? { |edge| incidence[edge] == 2 })
    assert_equal(incidence.keys.sort, result.edges)

    # Adjacency is symmetric and points across the opposite edge.
    result.faces.each_with_index do |face, face_index|
      3.times do |position|
        neighbour = result.opposite_face(face_index, position)
        assert(result.adjacency[neighbour].include?(face_index))
        shared = face.values_at((position + 1) % 3, (position + 2) % 3).sort
        neighbour_face = result.faces[neighbour]
        assert(neighbour_face.include?(shared[0]))
        assert(neighbour_face.include?(shared[1]))
      end
    end
  end

  # Empty spherical circum-cap: every other vertex lies on or below the
  # supporting plane of each face.
  def assert_empty_circumcaps(result, delta = 1e-8)
    result.faces.each do |face|
      normal = unit_normal(result, face)
      anchor = result.vertices[face[0]]
      support = normal.zip(anchor).sum { |x, y| x * y }
      result.vertices.each do |point|
        value = normal.zip(point).sum { |x, y| x * y }
        assert_in_delta(0.0, [value - support, 0.0].max, delta)
      end
    end
  end

  def test_tetrahedron_counts_euler_and_area
    result = SD.triangulate(regular_tetrahedron)
    assert_equal(4, result.face_count)
    assert_equal(6, result.edge_count)
    assert_spherical_complex(result, 4)
    assert_empty_circumcaps(result)
  end

  def test_octahedron
    result = SD.triangulate(octahedron)
    assert_equal(8, result.face_count)
    assert_equal(12, result.edge_count)
    assert_spherical_complex(result, 6)
    assert_empty_circumcaps(result)
    assert(result.degeneracies.empty?)
  end

  def test_icosahedron
    result = SD.triangulate(icosahedron)
    assert_equal(20, result.face_count)
    assert_equal(30, result.edge_count)
    assert_spherical_complex(result, 12)
    assert_empty_circumcaps(result)
    # A regular icosahedron has no four co-circular vertices.
    assert(result.degeneracies.empty?)
  end

  def test_cube_cocircular_fans_are_deterministic
    result = SD.triangulate(cube)
    assert_equal(12, result.face_count)
    assert_equal(18, result.edge_count)
    assert_spherical_complex(result, 8)
    assert_empty_circumcaps(result)
    kinds = result.degeneracies.map(&:kind)
    assert(kinds.all? { |kind| kind == :cocircular })
    assert_equal(6, kinds.length)
    assert(!result.unreliable)
  end

  def test_random_point_sets
    [[40, 1], [80, 2], [120, 7]].each do |count, seed|
      result = SD.triangulate(random_sphere(count, seed))
      assert_spherical_complex(result, count)
      assert_empty_circumcaps(result)
      assert(!result.unreliable)
    end
  end

  def test_topology_is_independent_of_input_order
    points = random_sphere(60, 11)
    baseline = SD.triangulate(points)
    [3, 17, 4096].each do |seed|
      permutation = (0...points.length).to_a
      random = Random.new(seed)
      (points.length - 1).downto(1) do |i|
        j = random.rand(i + 1)
        permutation[i], permutation[j] = permutation[j], permutation[i]
      end
      shuffled_result = SD.triangulate(permutation.map { |i| points[i] })
      forward = shuffled_result.original_to_vertex.each_with_index.map do
        |canonical, input_index|
        [canonical, baseline.original_to_vertex[permutation[input_index]]]
      end.to_h
      baseline_faces = baseline.faces.map(&:sort).sort
      shuffled_faces = shuffled_result.faces.map do |face|
        face.map { |vertex| forward[vertex] }.sort
      end.sort
      assert_equal(baseline_faces, shuffled_faces)
    end
  end

  def test_duplicate_identity_mapping
    tetrahedron = regular_tetrahedron
    result = SD.triangulate(tetrahedron + [tetrahedron[0], tetrahedron[2]])
    assert_equal(4, result.vertex_count)
    assert_equal(4, result.face_count)
    assert_equal((0..5).to_a, result.original_to_vertex.each_index.to_a)
    # Original indices 0 and 4 both map to the same canonical vertex.
    assert_equal(result.original_to_vertex[0], result.original_to_vertex[4])
    assert_equal(result.original_to_vertex[2], result.original_to_vertex[5])
    grouped = result.vertex_original_indices.map(&:sort).sort
    assert_include(grouped, [0, 4])
    assert_include(grouped, [2, 5])
    result.faces.each do |face|
      originals = result.face_original_indices(result.faces.index(face))
      assert_equal(3, originals.length)
    end
  end

  def test_too_few_points_is_inapplicable
    assert_raise(SD::InapplicableError) do
      SD.triangulate([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    end
  end

  def test_closed_hemisphere_is_inapplicable
    northern = [[1, 0, 0.1], [-1, 0, 0.1], [0, 1, 0.1], [0, -1, 0.1],
                [0, 0, 1.0], [0.1, 0, 0.05]].map { |point| normalize(point) }
    assert_raise(SD::InapplicableError) { SD.triangulate(northern) }
    equator = [[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]]
    assert_raise(SD::InapplicableError) { SD.triangulate(equator) }
  end

  def test_near_duplicate_is_reported
    tetrahedron = regular_tetrahedron
    near = tetrahedron.dup
    perturbed = tetrahedron[0].dup
    perturbed[0] += 5e-11
    near << normalize(perturbed)
    assert_raise(SD::DegenerateDataError) do
      SD.triangulate(near, strict: true)
    end
    begin
      SD.triangulate(near)
      flunk('expected DegenerateDataError for an undecidable near duplicate')
    rescue SD::DegenerateDataError => error
      assert(error.degeneracies.any? { |d| d.kind == :near_duplicate })
    end
  end

  def test_dateline_and_poles_are_seams
    stations = [[0, 180], [0, -180], [60, 30], [60, 150], [-60, 60],
                [-60, -60], [30, -90], [-30, 120]]
    result = SD.triangulate(stations)
    assert_equal(7, result.vertex_count)
    assert_spherical_complex(result, 7)
    assert_equal([:duplicates], result.degeneracies.map(&:kind).uniq)

    rotated = stations.map { |lat, lng| [lat, lng + 137.0] }
    rotated_result = SD.triangulate(rotated)
    assert_equal(result.face_count, rotated_result.face_count)
    assert_in_delta(FOUR_PI, rotated_result.total_area, 1e-7)
  end

  def test_latlng_objects_and_pairs_are_accepted
    stations = [Geokit::LatLng.new(10, 20), [-30, 80], [60, -100],
                [20, 170], [-70, -20], [5, -90]]
    result = SD.triangulate(stations)
    assert_equal(6, result.vertex_count)
    assert_spherical_complex(result, 6)
    assert_empty_circumcaps(result)
  end

  def test_zero_vector_is_rejected
    assert_raise(ArgumentError) { SD.to_unit_vector([0, 0, 0]) }
    assert_raise(ArgumentError) { SD.to_unit_vector(nil) }
  end

  def test_directed_edges_are_consistent_with_orientation
    result = SD.triangulate(regular_tetrahedron)
    result.faces.each_with_index do |face, face_index|
      directed = result.directed_edges(face_index)
      expected = [[face[0], face[1]], [face[1], face[2]], [face[2], face[0]]]
      assert_equal(expected, directed)
      directed.each do |edge|
        reversed = [edge[1], edge[0]]
        neighbours = result.faces.each_index.reject { |i| i == face_index }
        assert(neighbours.any? { |i| result.directed_edges(i).include?(
          reversed) })
      end
    end
  end
end
