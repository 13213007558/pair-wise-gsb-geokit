require File.join(File.dirname(__FILE__), 'helper')

class MapMatchingTest < Test::Unit::TestCase #:nodoc: all
  MM = Geokit::MapMatching

  # A tiny in-memory directed-graph router used to answer road-network
  # distances for the matcher. Topology is described per edge id as
  # {from:, to:, oneway:}; travel direction follows the edge geometry.
  class GraphRouter
    def initialize(edges, topology)
      @edges = edges
      @topology = topology
      @adj = Hash.new { |hash, key| hash[key] = [] }
      topology.each do |id, topo|
        length = edges[id].length_m
        @adj[topo[:from]] << [topo[:to], length]
        @adj[topo[:to]] << [topo[:from], length] unless topo[:oneway]
      end
    end

    def route_distance(from_candidate, to_candidate)
      from_edge = from_candidate.edge
      to_edge = to_candidate.edge
      if from_edge.id == to_edge.id
        delta = (to_candidate.fraction - from_candidate.fraction) * from_edge.length_m
        return delta.abs unless @topology[from_edge.id][:oneway]
        return delta if delta >= 0
      end
      tail = (1.0 - from_candidate.fraction) * from_edge.length_m
      head = to_candidate.fraction * to_edge.length_m
      mid = shortest_path(@topology[from_edge.id][:to], @topology[to_edge.id][:from])
      mid.nil? ? nil : tail + mid + head
    end

    private

    def shortest_path(from, to)
      dist = { from => 0.0 }
      visited = {}
      until dist.empty?
        node, d = dist.min_by { |_, v| v }
        dist.delete(node)
        next if visited[node]
        visited[node] = true
        return d if node == to
        @adj[node].each do |neighbor, weight|
          next if visited[neighbor]
          nd = d + weight
          dist[neighbor] = nd if dist[neighbor].nil? || nd < dist[neighbor]
        end
      end
      nil
    end
  end

  def edge(id, geometry)
    MM::Edge.new(id, geometry)
  end

  # Small two-hop network: A(0,0) --e1--> B(0,0.001) --e2--> C(0.001,0.001),
  # plus a diagonal shortcut e3 from A to C. All edges bidirectional.
  def small_network
    edges = {
      e1: edge(:e1, [[0, 0], [0, 0.001]]),
      e2: edge(:e2, [[0, 0.001], [0.001, 0.001]]),
      e3: edge(:e3, [[0, 0], [0.001, 0.001]])
    }
    topology = {
      e1: { from: :A, to: :B },
      e2: { from: :B, to: :C },
      e3: { from: :A, to: :C }
    }
    [edges, GraphRouter.new(edges, topology)]
  end

  def test_edge_length_and_projection
    e = edge(:e, [[0, 0], [0, 0.001]])
    assert_in_delta 111.30, e.length_m, 0.01

    candidate = e.project([0.00005, 0.0005])
    assert_in_delta 0.0, candidate.point.lat, 1e-9
    assert_in_delta 0.0005, candidate.point.lng, 1e-9
    assert_in_delta 5.56, candidate.distance, 0.01
    assert_in_delta 0.5, candidate.fraction, 1e-6
    assert_same e, candidate.edge
  end

  def test_projection_clamps_to_endpoints
    e = edge(:e, [[0, 0], [0, 0.001]])
    before_start = e.project([0.00001, -0.0005])
    assert_in_delta 0.0, before_start.fraction, 1e-9
    after_end = e.project([0.0, 0.0015])
    assert_in_delta 1.0, after_end.fraction, 1e-9
  end

  def test_probability_models
    assert_in_delta(-Math.log(10.0 * Math.sqrt(2.0 * Math::PI)),
                    MM.emission_log_probability(0.0, 10.0), 1e-12)
    assert MM.emission_log_probability(5.0, 10.0) >
           MM.emission_log_probability(50.0, 10.0)

    assert_in_delta(-Math.log(20.0),
                    MM.transition_log_probability(100.0, 100.0, 20.0), 1e-12)
    assert_in_delta(-Math.log(20.0) - 2.5,
                    MM.transition_log_probability(150.0, 100.0, 20.0), 1e-12)
  end

  def test_matcher_requires_explicit_valid_parameters
    router = ->(_a, _b) { 0.0 }
    assert_raise(ArgumentError) { MM::Matcher.new(router: router, sigma: 0, beta: 10) }
    assert_raise(ArgumentError) { MM::Matcher.new(router: router, sigma: 5, beta: -1) }
    assert_raise(ArgumentError) { MM::Matcher.new(router: router, sigma: 5, beta: 10, max_time_gap: 0) }
    assert_raise(ArgumentError) { MM::Matcher.new(router: Object.new, sigma: 5, beta: 10) }
  end

  # Enumerates every candidate path through the tiny network, scores each
  # one independently with the public probability functions and checks
  # that the Viterbi matcher returns the global optimum.
  def test_viterbi_matches_brute_force_enumeration
    edges, router = small_network
    sigma = 8.0
    beta = 30.0
    observations = [
      [[0.00002, 0.0001], 0],
      [[-0.00002, 0.00098], 30],
      [[0.00098, 0.00102], 60]
    ]
    lists = [
      [edges[:e1], edges[:e3]],
      [edges[:e1], edges[:e2]],
      [edges[:e2], edges[:e3]]
    ].each_with_index.map do |edge_list, i|
      edge_list.map { |e| e.project(observations[i][0]) }
    end

    linear = observations.each_cons(2).map do |a, b|
      MM.distance_m(a[0], b[0])
    end

    best_score = -Float::INFINITY
    best_path = nil
    lists[0].product(lists[1], lists[2]).each do |path|
      next unless path.each_cons(2).all? { |a, b| !router.route_distance(a, b).nil? }
      score = path.sum { |c| MM.emission_log_probability(c.distance, sigma) }
      path.each_cons(2).with_index do |(a, b), i|
        route_m = router.route_distance(a, b)
        score += MM.transition_log_probability(route_m, linear[i], beta)
      end
      if score > best_score
        best_score = score
        best_path = path
      end
    end

    matcher = MM::Matcher.new(router: router, sigma: sigma, beta: beta)
    result = matcher.match(observations, lists)

    assert result.success?
    assert_equal 1, result.segments.size
    segment = result.segments.first
    assert_equal 3, segment.size
    assert_equal best_path.map { |c| c.edge.id },
                 segment.edges.map(&:id)
    assert_in_delta best_score, segment.score, 1e-6
    assert_in_delta best_score, result.score, 1e-6
    assert_equal [0, 1, 2], segment.matches.map(&:index)
  end

  # Two parallel one-way roads 11 m apart carrying opposite directions,
  # connected at both ends so the wrong-direction road is reachable only
  # via a long loop. Noisy observations (sigma 8 m) drift across both
  # roads; the matcher must stay on the road that matches the travel
  # direction instead of jumping between the parallel roads.
  def test_parallel_roads_do_not_jump
    edges = {
      r1: edge(:r1, [[0, 0], [0, 0.004]]),
      r2: edge(:r2, [[0.0001, 0.004], [0.0001, 0]]),
      c1: edge(:c1, [[0.0001, 0], [0, 0]]),
      c2: edge(:c2, [[0, 0.004], [0.0001, 0.004]])
    }
    topology = {
      r1: { from: :P, to: :Q, oneway: true },
      r2: { from: :R, to: :S, oneway: true },
      c1: { from: :S, to: :P, oneway: true },
      c2: { from: :Q, to: :R, oneway: true }
    }
    router = GraphRouter.new(edges, topology)

    noise = [0.0, 0.00005, -0.00003, 0.00006, -0.00005, 0.00002]
    observations = noise.each_with_index.map do |lat_offset, i|
      [[lat_offset, i * 0.0008], i * 10]
    end
    lists = observations.map do |(point, _)|
      [edges[:r1].project(point), edges[:r2].project(point)]
    end

    matcher = MM::Matcher.new(router: router, sigma: 8.0, beta: 20.0)
    result = matcher.match(observations, lists)

    assert_equal 1, result.segments.size
    segment = result.segments.first
    assert_equal 6, segment.size
    assert_equal [:r1] * 6, segment.edges.map(&:id)
    segment.matches.each do |match|
      assert_in_delta 0.0, match.point.lat, 1e-9
    end
  end

  def test_observation_without_candidates_is_unmatched_and_splits
    edges, router = small_network
    observations = [
      [[0.0, 0.0002], 0],
      [[0.0, 0.0005], 10],
      [[0.0, 0.0008], 20]
    ]
    lists = [
      [edges[:e1].project(observations[0][0])],
      [],
      [edges[:e1].project(observations[2][0])]
    ]

    matcher = MM::Matcher.new(router: router, sigma: 5.0, beta: 20.0)
    result = matcher.match(observations, lists)

    assert_equal 2, result.segments.size
    assert_equal [1, 1], result.segments.map(&:size)
    assert_not_nil result.matches[0]
    assert_nil result.matches[1]
    assert_not_nil result.matches[2]
  end

  def test_duplicate_timestamps_start_a_new_segment
    e = edge(:e, [[0, 0], [0, 0.004]])
    router = GraphRouter.new({ e: e }, e: { from: :A, to: :B })
    observations = [
      [[0, 0.0005], 0],
      [[0, 0.0010], 10],
      [[0, 0.0015], 10], # duplicate timestamp
      [[0, 0.0020], 20]
    ]
    lists = observations.map { |(point, _)| [e.project(point)] }

    matcher = MM::Matcher.new(router: router, sigma: 5.0, beta: 20.0)
    result = matcher.match(observations, lists)

    assert_equal 2, result.segments.size
    assert_equal [0, 1], result.segments[0].matches.map(&:index)
    assert_equal [2, 3], result.segments[1].matches.map(&:index)
    assert result.matches.all?
  end

  def test_unreachable_edges_are_not_connected_through_the_network
    edges = {
      ex: edge(:ex, [[0, 0], [0, 0.001]]),
      ey: edge(:ey, [[0.01, 0.01], [0.01, 0.011]])
    }
    # Two disconnected components: no route exists between ex and ey.
    router = GraphRouter.new(edges, ex: { from: :A, to: :B }, ey: { from: :C, to: :D })
    observations = [
      [[0, 0.0003], 0],
      [[0, 0.0007], 10],
      [[0.01, 0.0103], 20],
      [[0.01, 0.0107], 30]
    ]
    lists = [
      [edges[:ex].project(observations[0][0])],
      [edges[:ex].project(observations[1][0])],
      [edges[:ey].project(observations[2][0])],
      [edges[:ey].project(observations[3][0])]
    ]

    matcher = MM::Matcher.new(router: router, sigma: 5.0, beta: 20.0)
    result = matcher.match(observations, lists)

    assert_equal 2, result.segments.size
    assert_equal %i[ex ex], result.segments[0].edges.map(&:id)
    assert_equal %i[ey ey], result.segments[1].edges.map(&:id)
  end

  def test_long_signal_loss_splits_only_when_max_time_gap_set
    e = edge(:e, [[0, 0], [0, 0.01]])
    router = GraphRouter.new({ e: e }, e: { from: :A, to: :B })
    observations = [
      [[0, 0.001], 0],
      [[0, 0.002], 10],
      [[0, 0.008], 10_000] # reappears after a long signal loss
    ]
    lists = observations.map { |(point, _)| [e.project(point)] }

    with_gap = MM::Matcher.new(router: router, sigma: 5.0, beta: 50.0,
                               max_time_gap: 600)
    result = with_gap.match(observations, lists)
    assert_equal 2, result.segments.size
    assert_equal [2, 1], result.segments.map(&:size)

    without_gap = MM::Matcher.new(router: router, sigma: 5.0, beta: 50.0)
    result = without_gap.match(observations, lists)
    assert_equal 1, result.segments.size
    assert_equal 3, result.segments.first.size
  end

  def test_one_way_street_direction_is_respected
    e = edge(:ow, [[0, 0], [0, 0.002]])
    router = GraphRouter.new({ ow: e }, ow: { from: :A, to: :B, oneway: true })
    matcher = MM::Matcher.new(router: router, sigma: 5.0, beta: 20.0)

    forward = [
      [[0, 0.0005], 0],
      [[0, 0.0010], 10],
      [[0, 0.0015], 20]
    ]
    lists = forward.map { |(point, _)| [e.project(point)] }
    result = matcher.match(forward, lists)
    assert_equal 1, result.segments.size
    assert_equal 3, result.segments.first.size

    backward = forward.reverse.each_with_index.map do |(point, _), i|
      [point, i * 10]
    end
    lists = backward.map { |(point, _)| [e.project(point)] }
    result = matcher.match(backward, lists)
    # Driving against the one-way is unreachable and there is no way
    # back, so every observation becomes its own segment instead of a
    # single impossible route.
    assert_equal 3, result.segments.size
    assert_equal [1, 1, 1], result.segments.map(&:size)
  end

  def test_observation_input_forms_and_empty_input
    e = edge(:e, [[0, 0], [0, 0.001]])
    router = GraphRouter.new({ e: e }, e: { from: :A, to: :B })
    matcher = MM::Matcher.new(router: router, sigma: 5.0, beta: 20.0)

    observations = [
      MM::Observation.new([0, 0.0002], 0),
      { point: [0, 0.0005], time: 10 },
      [[0, 0.0008], 20]
    ]
    lists = observations.map do |o|
      point = o.is_a?(MM::Observation) ? o.point : (o.is_a?(Hash) ? o[:point] : o[0])
      [e.project(point)]
    end
    result = matcher.match(observations, lists)
    assert_equal 1, result.segments.size
    assert_equal 3, result.segments.first.size

    result = matcher.match([], [])
    assert_equal false, result.success?
    assert_equal [], result.segments

    assert_raise(ArgumentError) { matcher.match([[[0, 0], 0]], []) }
  end
end
