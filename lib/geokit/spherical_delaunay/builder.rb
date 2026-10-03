# frozen_string_literal: true

class Geokit::SphericalDelaunay::Builder # :nodoc: all
  FLOAT_EPSILON = 1e-13
  class HullFailure < StandardError; end

  def initialize(points, tolerance, strict)
    @input = points
    @tolerance = tolerance
    @strict = strict
    @degeneracies = []
  end

  def build
    prepare_vertices
    detect_near_duplicates
    build_geometry
    verify_full_hemisphere
    result = assemble_result
    fail_if_unreliable
    result
  end

  private

  def prepare_vertices
    raise ArgumentError, 'points must respond to each' unless
      @input.respond_to?(:each)
    vectors = @input.each_with_index.map do |point, index|
      [Geokit::SphericalDelaunay.to_unit_vector(point), index]
    end
    raise Geokit::SphericalDelaunay::InapplicableError,
          'at least four points are required' if vectors.length < 4

    groups = {}
    vectors.each do |vector, original_index|
      (groups[vector] ||= []) << original_index
    end
    @vectors = groups.keys.sort { |left, right| lex_compare(left, right) }
    @vertex_original_indices = @vectors.map { |vector| groups[vector].sort }
    if @vertex_original_indices.sum(&:length) != vectors.length ||
       @vectors.length < vectors.length
      add_degeneracy(:duplicates,
                     "#{vectors.length - @vectors.length} exact duplicate " +
                     'input(s) merged into canonical vertices', [])
    end
    raise Geokit::SphericalDelaunay::InapplicableError,
          'at least four distinct points are required' if @vectors.length < 4

    @original_to_vertex = Array.new(vectors.length)
    @vertex_original_indices.each_with_index do |originals, vertex|
      originals.each { |original| @original_to_vertex[original] = vertex }
    end
  end

  def lex_compare(left, right)
    3.times do |axis|
      difference = left[axis] <=> right[axis]
      return difference unless difference.zero?
    end
    0
  end

  def detect_near_duplicates
    limit = 2.0 * Math.sin(@tolerance / 2.0)
    count = @vectors.length
    count.times do |i|
      (i + 1).upto(count - 1) do |j|
        next unless chord_distance(@vectors[i], @vectors[j]) <= limit
        add_unreliable(:near_duplicate,
                       "vertices #{i} and #{j} are within the tolerance of " +
                       'being identical', [i, j])
      end
    end
  end

  def chord_distance(left, right)
    squared = 3.times.map do |axis|
      difference = left[axis] - right[axis]
      difference * difference
    end.inject(0.0, :+)
    Math.sqrt(squared)
  end

  # Face growth (gift-wrap) over exact supporting planes. Growth depends
  # only on the geometry and canonical vertex indices, never on input order.
  def build_geometry
    recover_supporting_planes
    assemble_manifold
  rescue HullFailure
    diagnose_failure
  rescue Geokit::SphericalDelaunay::DegenerateDataError
    diagnose_failure
  end

  # Distinguish "cannot surround the sphere" from "too degenerate to decide".
  def diagnose_failure
    distance = distance_to_hull_from_origin
    unless distance <= @tolerance
      raise Geokit::SphericalDelaunay::InapplicableError,
            'points lie in a closed hemisphere; a full-sphere ' +
            'triangulation is not applicable'
    end
    if distance > 0.0
      add_unreliable(:closed_hemisphere_boundary,
                     "the convex hull approaches the origin within the " +
                     "tolerance (distance #{distance}); the input is on the " +
                     'verge of lying in a closed hemisphere', [])
    end
    flagged = @degeneracies.select do |degeneracy|
      Geokit::SphericalDelaunay::UNRELIABLE_KINDS.include?(degeneracy.kind)
    end
    unless flagged.empty?
      raise Geokit::SphericalDelaunay::DegenerateDataError.new(
        flagged.map(&:message).join('; '), flagged)
    end
    raise Geokit::SphericalDelaunay::DegenerateDataError.new(
      'could not grow a closed hull; the data is too close to a degenerate ' +
      'configuration to decide reliably', @degeneracies)
  end

  # Full-sphere applicability via the distance from the origin to the
  # convex hull of the points. The origin must be strictly inside the hull;
  # a positive distance means the whole set fits in an open hemisphere, and
  # a distance within the tolerance means it effectively meets a great
  # circle (closed hemisphere). Either case makes triangulation inapplicable.
  def origin_strictly_inside_hull?
    distance_to_hull_from_origin <= @tolerance
  end

  # Minimum distance from the origin to any point in the convex hull. This is
  # a convex quadratic program over the simplex; a deterministic cyclic
  # coordinate descent on barycentric weights converges to the optimum.
  def distance_to_hull_from_origin
    count = @vectors.length
    weights = Array.new(count, 0.0)
    @vectors.each_with_index do |point, index|
      weights[index] = 1.0
      residual = hull_point(weights)
      return 0.0 if norm3(residual) <= @tolerance
    end
    weights = Array.new(count, 1.0 / count)
    best = norm3(hull_point(weights))
    10_000.times do
      moved = false
      count.times do |from|
        next if weights[from] <= 0.0
        count.times do |to|
          next if from == to
          shift = optimal_shift(from, to, weights)
          next if shift <= 1e-14
          weights[from] -= shift
          weights[to] += shift
          moved = true
        end
      end
      current = norm3(hull_point(weights))
      best = [best, current].min
      break if !moved || current <= @tolerance
    end
    [best, norm3(hull_point(weights))].min
  end

  def hull_point(weights)
    point = [0.0, 0.0, 0.0]
    weights.each_with_index do |weight, index|
      3.times { |axis| point[axis] += weight * @vectors[index][axis] }
    end
    point
  end

  # Optimal non-negative transfer of weight from one vertex to another while
  # reducing the norm of the current convex combination.
  def optimal_shift(from, to, weights)
    current = hull_point(weights)
    direction = diff_vector(@vectors[to], @vectors[from])
    numerator = dot(current, direction)
    denominator = dot(direction, direction)
    return 0.0 if numerator >= 0.0 || denominator.zero?
    [-numerator / denominator, weights[from]].min
  end

  def norm3(vector)
    Math.sqrt(dot(vector, vector))
  end

  def recover_supporting_planes
    @planes = {}
    add_seed_plane
    gift_wrap_planes
    @plane_records = @planes.values
  end

  # A lexicographically extreme vertex is necessarily a hull vertex; scan
  # the triples through it for the plane with the largest support distance.
  def add_seed_plane
    pivot = (0...@vectors.length).max_by do |vertex|
      [@vectors[vertex][0], @vectors[vertex][1], @vectors[vertex][2]]
    end
    candidates = collect_seed_candidates(pivot)
    if candidates.empty?
      raise Geokit::SphericalDelaunay::InapplicableError,
            'no supporting plane found; points may all lie in a closed '
            'hemisphere or on one great circle'
    end
    support, second, third = candidates.max_by { |entry|
      [entry[0], -entry[1], -entry[2]] }
    best = [pivot, second, third]
    if support.nil? || support <= @tolerance
      raise Geokit::SphericalDelaunay::InapplicableError,
            'points may all lie in a closed hemisphere or on one great '
            'circle'
    end
    add_plane_from_triple(best)
  end

  def collect_seed_candidates(pivot)
    @vectors.each_index.flat_map do |second|
      next [] if second == pivot
      @vectors.each_index.flat_map do |third|
        next [] if third == pivot || third == second
        triple = [pivot, second, third]
        next [] unless non_degenerate_triple?(triple)
        normal, support = normalised_plane(triple)
        next [] unless supporting_plane?(normal, support, triple)
        [[support, second, third]]
      end
    end
  end

  def non_degenerate_triple?(triple)
    a, b, c = triple.map { |vertex| @vectors[vertex] }
    normal = cross(diff_vector(b, a), diff_vector(c, a))
    dot(normal, normal) > 0.0
  end

  def gift_wrap_planes
    loop do
      frontier = boundary_edges
      break if frontier.empty?
      progressed = false
      frontier.each do |edge|
        next if add_plane_across_edge(edge) == :already
        progressed = true
      end
      raise HullFailure unless progressed
    end
  end

  def boundary_edges
    incidence = Hash.new(0)
    @planes.each_value do |record|
      polygon_edges(record[:vertices]).each { |edge| incidence[edge] += 1 }
    end
    incidence.select { |_, count| count == 1 }.keys
  end

  def polygon_edges(vertices)
    edges = vertices.each_cons(2).map { |a, b| [a, b].sort }
    edges << [vertices.first, vertices.last].sort
  end

  # Grow the unique new supporting plane incident to a boundary edge.
  def add_plane_across_edge(edge)
    candidates = supporting_triples_across(edge)
    return :already if candidates.empty?
    groups = []
    candidates.each do |candidate|
      group = groups.find { |entry|
        same_plane?(candidate[:normal], candidate[:support],
                    entry.first[:normal], entry.first[:support]) }
      if group
        group << candidate
      else
        groups << [candidate]
      end
    end
    if groups.length > 1
      supports = groups.map { |group| group.first[:support] }
      if supports.max - supports.min <= @tolerance
        add_unreliable(:near_cocircular,
                       "edge #{edge.inspect} has adjacent planes within the " +
                       'tolerance of each other', edge)
      end
    end
    chosen = groups.min_by { |group|
      [-group.first[:support]]
    }.first
    add_plane_from_triple(chosen[:triple])
    :added
  end

  def supporting_triples_across(edge)
    results = []
    seen = []
    @vectors.each_index do |candidate|
      next if edge.include?(candidate)
      triple = [edge[0], edge[1], candidate]
      next unless non_degenerate_triple?(triple)
      normal, support = normalised_plane(triple)
      next if support <= @tolerance
      next unless supporting_plane?(normal, support, triple)
      next if find_plane(normal, support, exact_plane_signature(triple))
      next if seen.any? { |entry|
        same_plane?(normal, support, entry[:normal], entry[:support]) }
      seen << { normal: normal, support: support }
      results << { triple: triple, normal: normal, support: support }
    end
    results
  end

  def add_plane_from_triple(face)
    normal, support = normalised_plane(face)
    exact = exact_plane_signature(face)
    existing = find_plane(normal, support, exact)
    return @planes[existing] if existing
    signature = { exact: exact }
    on_plane = []
    near_plane = []
    @vectors.each_index do |vertex|
      residual = float_plane_residual(normal, support, vertex)
      if residual.abs <= [@tolerance, 1e-8].max
        if residual.abs <= FLOAT_EPSILON ||
             exact_plane_value(face, vertex).zero?
          # Points agreeing to floating-point round-off are treated as
          # exactly coplanar (e.g. a latitude ring built from sin/cos).
          on_plane << vertex
        elsif residual.abs <= @tolerance
          near_plane << vertex
        end
      end
    end
    return unless supporting_plane?(normal, support, on_plane)
    on_plane = convex_polygon(face, on_plane)
    return if on_plane.length < 3
    near_plane.each do |vertex|
      add_unreliable(:near_cocircular,
                     "vertex #{vertex} is within the tolerance of the " +
                     'supporting circle of a face',
                     [vertex] + on_plane.first(3))
    end
    if on_plane.length > 3
      add_degeneracy(:cocircular,
                     "#{on_plane.length} co-circular vertices fanned " +
                     'deterministically', on_plane)
    end
    @planes[signature] = {
      normal: normal,
      support: support,
      vertices: on_plane,
      signature: signature
    }
  end

  # True when two outward oriented planes agree to round-off.
  def same_plane?(left_normal, left_support, right_normal, right_support)
    normal_gap = 3.times.map do |axis|
      (left_normal[axis] - right_normal[axis]).abs
    end.max
    normal_gap <= FLOAT_EPSILON &&
      (left_support - right_support).abs <= FLOAT_EPSILON
  end

  def find_plane(normal, support, exact)
    @planes.each_value do |record|
      return record[:signature] if record[:signature][:exact] == exact
      if same_plane?(normal, support, record[:normal], record[:support])
        return record[:signature]
      end
    end
    nil
  end

  def supporting_plane?(normal, support, on_plane)
    return false if support <= 0.0
    @vectors.each_index do |vertex|
      next if on_plane.include?(vertex)
      value = dot(normal, @vectors[vertex])
      return false if value > support + [@tolerance, 1e-8].max
    end
    true
  end

  def float_plane_residual(normal, support, vertex)
    dot(normal, @vectors[vertex]) - support
  end

  def exact_plane_signature(face)
    rows = face.map { |vertex| @vectors[vertex].map(&:to_r) }
    edge1 = vector_subtract(rows[1], rows[0])
    edge2 = vector_subtract(rows[2], rows[0])
    rational_normal = cross(edge1, edge2)
    support = dot(rational_normal, rows[0])
    coefficients = rational_normal.map(&:to_r) << -support
    common_denominator = coefficients.map(&:denominator).inject(1, :lcm)
    integers = coefficients.map { |value| (value * common_denominator).to_i }
    first_nonzero = integers.find { |value| !value.zero? }
    integers = integers.map { |value| -value } if first_nonzero < 0
    divisor = integers.map(&:abs).inject(:gcd)
    integers.map { |value| value / divisor }
  end

  def normalised_plane(face)
    edge1 = diff_vector(@vectors[face[1]], @vectors[face[0]])
    edge2 = diff_vector(@vectors[face[2]], @vectors[face[0]])
    normal = cross(edge1, edge2)
    length = Math.sqrt(dot(normal, normal))
    normal = scale_vector(normal, 1.0 / length)
    support = dot(normal, @vectors[face[0]])
    if support < 0.0
      normal = scale_vector(normal, -1.0)
      support = -support
    end
    [normal, support]
  end

  def exact_plane_value(face, vertex)
    a, b, c, p = [face[0], face[1], face[2], vertex].map do |index|
      @vectors[index].map(&:to_r)
    end
    edge1 = vector_subtract(b, a)
    edge2 = vector_subtract(c, a)
    dot(cross(edge1, edge2), vector_subtract(p, a))
  end

  # Vertices of the convex polygon on a plane, ordered counter-clockwise as
  # seen from outside the sphere.
  def convex_polygon(face, vertices)
    normal, = normalised_plane(face)
    axis = most_orthogonal_axis(normal)
    first = basis_vector(normal, axis)
    second = cross(normal, first)
    points = vertices.map do |vertex|
      vector = @vectors[vertex]
      [dot(vector, first), dot(vector, second), vertex]
    end
    anchor = points.min_by { |x, y, vertex| [x, y, vertex] }
    # second = normal x first makes (first, second) right-handed when seen
    # from the side normal points to (the outside), so increasing atan2(y,x)
    # is counter-clockwise order. Every point exactly on the sphere/plane
    # intersection lies on the same circle, so angular order traces the
    # convex polygon boundary and no general hull construction is needed.
    ordered = points.sort_by { |x, y, _| Math.atan2(y, x) }
    offset = ordered.index(anchor)
    (ordered[offset..-1] + ordered[0...offset]).map(&:last)
  end

  def most_orthogonal_axis(normal)
    components = normal.map(&:abs)
    components.index(components.min)
  end

  def basis_vector(normal, axis)
    candidate = [0.0, 0.0, 0.0]
    candidate[axis] = 1.0
    first = diff_vector(candidate,
                        scale_vector(normal, dot(candidate, normal)))
    scale_vector(first, 1.0 / Math.sqrt(dot(first, first)))
  end

  # Rebuild a canonical triangle fan for every supporting plane and connect
  # the manifold using edge incidence.
  def assemble_manifold
    records = @plane_records.sort_by do |record|
      record[:vertices].min
    end
    @faces = records.flat_map do |record|
      fan_triangles(record[:vertices])
    end
    raise HullFailure if @faces.empty?
    unless all_vertices_covered?
      ensure_vertex_coverage
      raise HullFailure
    end
    undirected = Hash.new { |hash, key| hash[key] = [] }
    @faces.each_with_index do |face, face_index|
      face_edges(face).each do |edge|
        undirected[edge.sort] << face_index
      end
    end
    bad = undirected.select { |_, incident| incident.length != 2 }
    unless bad.empty?
      raise Geokit::SphericalDelaunay::DegenerateDataError.new(
        "recovered hull is not a closed manifold; bad edges #{bad.keys.inspect}",
        @degeneracies)
    end
    @adjacency = @faces.each_with_index.map do |face, face_index|
      3.times.map do |position|
        edge = [face[(position + 1) % 3], face[(position + 2) % 3]].sort
        undirected[edge].reject { |index| index == face_index }.first
      end
    end
    @undirected_edges = undirected.keys.sort
  end

  def directed_edges(face)
    [[face[0], face[1]], [face[1], face[2]], [face[2], face[0]]]
  end

  def face_edges(face)
    directed_edges(face)
  end

  def fan_triangles(ordered)
    (1...ordered.length - 1).map do |index|
      [ordered[0], ordered[index], ordered[index + 1]]
    end
  end

  def ensure_vertex_coverage
    covered = @faces.flatten.uniq
    missing = (0...@vectors.length).to_a - covered
    missing.each do |vertex|
      add_unreliable(:unclassified_vertex,
                     "vertex #{vertex} could not be placed on a supporting " +
                     'plane', [vertex])
    end
  end

  def all_vertices_covered?
    covered = @faces.flatten
    (0...@vectors.length).all? { |vertex| covered.include?(vertex) }
  end

  # The origin is inside the hull exactly when every outward supporting
  # plane has strictly positive distance from it. A non-positive distance
  # means the points fit in a closed hemisphere.
  def verify_full_hemisphere
    minimum = Float::INFINITY
    @plane_records.each do |record|
      exact_support = exact_support(record)
      minimum = exact_support.to_f if exact_support.to_f < minimum
      next unless exact_support <= 0r
      raise Geokit::SphericalDelaunay::InapplicableError,
            'points lie in a closed hemisphere; a full-sphere triangulation '
            'is not applicable'
    end
    if minimum <= @tolerance
      raise Geokit::SphericalDelaunay::InapplicableError,
            "a supporting great circle is within #{@tolerance} of the origin; "
            'the input is effectively confined to a closed hemisphere'
    end
  end

  def exact_support(record)
    face = record[:vertices].first(3)
    rows = face.map { |vertex| @vectors[vertex].map(&:to_r) }
    edge1 = vector_subtract(rows[1], rows[0])
    edge2 = vector_subtract(rows[2], rows[0])
    rational_normal = cross(edge1, edge2)
    first_nonzero = rational_normal.find { |value| !value.zero? }
    rational_normal = rational_normal.map { |value| value / first_nonzero }
    dot(rational_normal, rows[0]).abs
  end

  def assemble_result
    area = @faces.map { |face| spherical_triangle_area(face) }
                 .inject(0.0, :+)
    euler = @vectors.length - @undirected_edges.length + @faces.length
    unless euler == 2
      raise Geokit::SphericalDelaunay::DegenerateDataError.new(
        "Euler characteristic is #{euler}, expected 2", @degeneracies)
    end
    unless (area - 4.0 * Math::PI).abs <= [100.0 * @tolerance, 1e-6].max
      add_unreliable(:area,
                     "total spherical area #{area} differs from 4*PI", [])
    end
    unreliable = @degeneracies.any? do |degeneracy|
      Geokit::SphericalDelaunay::UNRELIABLE_KINDS.include?(degeneracy.kind)
    end
    mappings = {
      vertex_original_indices: @vertex_original_indices,
      original_to_vertex: @original_to_vertex
    }
    Geokit::SphericalDelaunay::Result.new(
      @vectors, @faces, @adjacency, @undirected_edges, mappings, area,
      @degeneracies, unreliable)
  end

  def fail_if_unreliable
    return unless @strict
    flagged = @degeneracies.select do |degeneracy|
      Geokit::SphericalDelaunay::UNRELIABLE_KINDS.include?(degeneracy.kind)
    end
    return if flagged.empty?
    raise Geokit::SphericalDelaunay::DegenerateDataError.new(
      flagged.map(&:message).join('; '), flagged)
  end

  def spherical_triangle_area(face)
    points = face.map { |vertex| @vectors[vertex] }
    angles = 3.times.map { |index| interior_angle(points, index) }
    excess = angles.inject(0.0, :+) - Math::PI
    excess = 0.0 if excess.abs < 1e-12
    excess
  end

  def interior_angle(points, index)
    origin = points[index]
    left = points[(index + 1) % 3]
    right = points[(index + 2) % 3]
    toward_left = diff_vector(left, scale_vector(origin, dot(origin, left)))
    toward_right = diff_vector(right, scale_vector(origin, dot(origin, right)))
    denominator = Math.sqrt(dot(toward_left, toward_left) *
                            dot(toward_right, toward_right))
    cosine = dot(toward_left, toward_right) / denominator
    cosine = [[cosine, -1.0].max, 1.0].min
    Math.acos(cosine)
  end

  def add_degeneracy(kind, message, vertices)
    @degeneracies << Geokit::SphericalDelaunay::Degeneracy.new(
      kind, message, vertices)
  end

  def add_unreliable(kind, message, vertices)
    add_degeneracy(kind, message, vertices)
  end

  def diff_vector(left, right)
    [left[0] - right[0], left[1] - right[1], left[2] - right[2]]
  end

  def vector_subtract(left, right)
    [left[0] - right[0], left[1] - right[1], left[2] - right[2]]
  end

  def scale_vector(vector, factor)
    vector.map { |value| value * factor }
  end

  def dot(left, right)
    left[0] * right[0] + left[1] * right[1] + left[2] * right[2]
  end

  def cross(left, right)
    [left[1] * right[2] - left[2] * right[1],
     left[2] * right[0] - left[0] * right[2],
     left[0] * right[1] - left[1] * right[0]]
  end
end
