require File.join(File.dirname(__FILE__), 'helper')

# A tiny directed-graph road network used as the matcher's network
# distance oracle. Nodes are lat/lng points; segments are registered
# between two nodes. route_distance follows the contract the matcher
# expects: Numeric distance in meters, or nil when unreachable.
class TinyNetwork
  SegmentInfo = Struct.new(:from_node, :to_node, :length)

  def initialize
    @nodes = {}
    @info = {}
    @adjacency = Hash.new { |hash, key| hash[key] = [] }
  end

  def add_node(name, lat, lng)
    @nodes[name] = Geokit::LatLng.new(lat, lng)
  end

  def add_segment(id, from_node, to_node, one_way: false)
    polyline = [@nodes[from_node], @nodes[to_node]]
    segment = Geokit::MapMatching::Segment.new(id, polyline, one_way)
    length = polyline.each_cons(2).sum do |a, b|
      a.distance_to(b, units: :meters)
    end
    @info[segment.object_id] = SegmentInfo.new(from_node, to_node, length)
    @adjacency[from_node] << [to_node, length]
    @adjacency[to_node] << [from_node, length] unless one_way
    segment
  end

  def route_distance(from_cand, to_cand)
    from = @info.fetch(from_cand.segment.object_id)
    to = @info.fetch(to_cand.segment.object_id)
    if from_cand.segment.equal?(to_cand.segment)
      delta = to_cand.offset - from_cand.offset
      return delta.abs unless from_cand.segment.one_way

      return delta >= 0 ? delta : nil
    end

    exits = [[from.to_node, from.length - from_cand.offset]]
    exits << [from.from_node, from_cand.offset] unless from_cand.segment.one_way
    entries = [[to.from_node, to_cand.offset]]
    entries << [to.to_node, to.length - to_cand.offset] unless to_cand.segment.one_way

    best = nil
    exits.each do |exit_node, exit_cost|
      entries.each do |entry_node, entry_cost|
        between = dijkstra(exit_node, entry_node)
        next if between.nil?

        total = exit_cost + between + entry_cost
        best = total if best.nil? || total < best
      end
    end
    best
  end

  private

  def dijkstra(start, goal)
    distances = { start => 0.0 }
    visited = {}
    until distances.empty?
      node, dist = distances.min_by { |_, d| d }
      distances.delete(node)
      visited[node] = dist
      return dist if node == goal

      @adjacency[node].each do |neighbour, weight|
        next if visited.key?(neighbour)

        candidate = dist + weight
        if !distances.key?(neighbour) || candidate < distances[neighbour]
          distances[neighbour] = candidate
        end
      end
    end
    nil
  end
end

class MapMatchingTest < Test::Unit::TestCase
  Observation = Geokit::MapMatching::Observation
  Candidate = Geokit::MapMatching::Candidate
  Matcher = Geokit::MapMatching::Matcher

  METERS_PER_LAT = Geokit::LatLng.units_per_latitude_degree(:meters)

  def meters_per_lng(lat)
    Geokit::LatLng.units_per_longitude_degree(lat, :meters)
  end

  # Independent re-implementation of the matcher's tangent-plane
  # projection, used to cross-check scores by enumeration.
  def project_for(point, segment)
    per_lat = METERS_PER_LAT
    per_lng = meters_per_lng(point.lat)
    xy = segment.polyline.map do |p|
      p = Geokit::LatLng.normalize(p)
      [(p.lng - point.lng) * per_lng, (p.lat - point.lat) * per_lat]
    end
    best = nil
    cumulative = 0.0
    xy.each_cons(2) do |(ax, ay), (bx, by)|
      dx = bx - ax
      dy = by - ay
      length_sq = dx * dx + dy * dy
      t = length_sq.zero? ? 0.0 : -(ax * dx + ay * dy) / length_sq
      t = [[t, 0.0].max, 1.0].min
      px = ax + t * dx
      py = ay + t * dy
      along = cumulative + Math.sqrt(length_sq) * t
      dist_sq = px * px + py * py
      best = [dist_sq, along] if best.nil? || dist_sq < best[0]
      cumulative += Math.sqrt(length_sq)
    end
    [Math.sqrt(best[0]), best[1]]
  end

  def emission(distance, sigma)
    -Math.log(sigma * Math.sqrt(2.0 * Math::PI)) - 0.5 * (distance / sigma)**2
  end

  def transition(route, displacement, beta)
    -Math.log(beta) - (route - displacement).abs / beta
  end

  # Main road along the equator, parallel side road ~55.6m north,
  # connected at the west ends.
  def build_connected_network
    network = TinyNetwork.new
    network.add_node(:a, 0.0, 0.0)
    network.add_node(:b, 0.0, 0.002)
    network.add_node(:c, 0.0005, 0.0)
    network.add_node(:d, 0.0005, 0.002)
    main = network.add_segment('main', :a, :b)
    side = network.add_segment('side', :c, :d)
    network.add_segment('link', :a, :c)
    [network, main, side]
  end

  def candidate_for(point, segment)
    distance, offset = project_for(point, segment)
    Candidate.new(segment: segment, point: point, distance: distance,
                  offset: offset, fraction: 0.0)
  end

  def test_matches_best_path_by_brute_force_enumeration
    network, main, side = build_connected_network
    sigma = 15.0
    beta = 40.0
    observations = [
      Observation.new(0.0, Geokit::LatLng.new(0.0001, 0.0005)),
      Observation.new(10.0, Geokit::LatLng.new(-0.0001, 0.0010)),
      Observation.new(20.0, Geokit::LatLng.new(0.00008, 0.0015))
    ]
    segments = [main, side]
    candidates = observations.map { segments }

    matcher = Matcher.new(network: network, sigma: sigma, beta: beta)
    result = matcher.match(observations, candidates)

    # Enumerate every candidate path and score it independently.
    scored_paths = segments.product(segments, segments).map do |choice|
      score = 0.0
      choice.each_with_index do |segment, i|
        distance, = project_for(observations[i].point, segment)
        score += emission(distance, sigma)
      end
      choice.each_cons(2).with_index do |(from_seg, to_seg), i|
        from_cand = candidate_for(observations[i].point, from_seg)
        to_cand = candidate_for(observations[i + 1].point, to_seg)
        route = network.route_distance(from_cand, to_cand)
        displacement = observations[i].point.distance_to(
          observations[i + 1].point, units: :meters
        )
        score += transition(route, displacement, beta)
      end
      [choice.map(&:id), score]
    end
    best_ids, best_score = scored_paths.max_by { |_, score| score }

    assert_equal 1, result.segments.length
    assert_equal best_ids, result.points.map { |p| p.segment.id }
    assert_in_delta best_score, result.segments.first.score, 1e-6
  end

  def test_projection_snaps_onto_segment
    network, main, = build_connected_network
    matcher = Matcher.new(network: network, sigma: 15.0, beta: 40.0)
    observation = Observation.new(0.0, Geokit::LatLng.new(0.0002, 0.001))
    result = matcher.match([observation], [[main]])

    point = result.points.first
    assert_equal 'main', point.segment.id
    assert_in_delta 0.0, point.point.lat, 1e-9 # snapped onto the road line
    assert_in_delta 0.001, point.point.lng, 1e-9
    assert_in_delta 0.0002 * METERS_PER_LAT, point.distance, 0.5
  end

  def test_parallel_roads_do_not_jump_back_and_forth
    network = TinyNetwork.new
    network.add_node(:a, 0.0, 0.0)
    network.add_node(:b, 0.0, 0.004)
    network.add_node(:c, -0.0003, 0.0)
    network.add_node(:d, -0.0003, 0.004)
    north = network.add_segment('north', :a, :b, one_way: true)
    south = network.add_segment('south', :c, :d, one_way: true)
    network.add_segment('west_link', :c, :a, one_way: true)
    network.add_segment('east_link', :b, :c, one_way: true)

    # True trajectory drives east along the north road. Points 3 and 4
    # are noisy enough to sit closer to the parallel south road.
    observations = (1..6).map do |k|
      lat = (3..4).cover?(k) ? -0.000225 : 0.0
      Observation.new(k * 5.0, Geokit::LatLng.new(lat, 0.0004 * k))
    end
    candidates = observations.map { [north, south] }

    # Sanity check the trap: the noisy middle points really are closer
    # to the south road, so greedy nearest-segment snapping would jump.
    noisy = observations[2].point
    south_distance, = project_for(noisy, south)
    north_distance, = project_for(noisy, north)
    assert south_distance < north_distance

    matcher = Matcher.new(network: network, sigma: 12.0, beta: 50.0)
    result = matcher.match(observations, candidates)

    assert_equal 1, result.segments.length
    assert_equal [], result.breaks
    assert_equal(['north'] * 6, result.points.map { |p| p.segment.id })
  end

  def test_long_signal_gap_splits_trajectory
    network, main, = build_connected_network
    matcher = Matcher.new(network: network, sigma: 15.0, beta: 40.0,
                          max_time_gap: 60.0)
    times = [0.0, 10.0, 20.0, 500.0, 510.0, 520.0]
    observations = times.each_with_index.map do |time, i|
      Observation.new(time, Geokit::LatLng.new(0.0, 0.0002 * (i + 1)))
    end
    candidates = observations.map { [main] }

    result = matcher.match(observations, candidates)

    assert_equal 2, result.segments.length
    assert_equal [3, 3], result.segments.map { |seg| seg.points.length }
    assert_equal [:time_gap], result.breaks.map(&:kind)
    assert_equal 500.0, result.breaks.first.observation.time
  end

  def test_observation_without_candidates_breaks_trajectory
    network, main, = build_connected_network
    matcher = Matcher.new(network: network, sigma: 15.0, beta: 40.0)
    times = [0.0, 10.0, 20.0, 30.0, 40.0]
    observations = times.each_with_index.map do |time, i|
      Observation.new(time, Geokit::LatLng.new(0.0, 0.0002 * (i + 1)))
    end
    candidates = observations.map { [main] }
    candidates[2] = []

    result = matcher.match(observations, candidates)

    assert_equal 2, result.segments.length
    assert_equal [2, 2], result.segments.map { |seg| seg.points.length }
    assert_equal [:no_candidates], result.breaks.map(&:kind)
    assert_equal 20.0, result.breaks.first.observation.time
    assert_equal [0.0, 10.0, 30.0, 40.0],
                 result.points.map { |p| p.observation.time }
  end

  def test_one_way_road_against_direction_is_unreachable
    network = TinyNetwork.new
    network.add_node(:a, 0.0, 0.0)
    network.add_node(:b, 0.0, 0.002)
    road = network.add_segment('one_way', :a, :b, one_way: true)

    matcher = Matcher.new(network: network, sigma: 15.0, beta: 40.0)
    # Driving west on an eastbound one-way road.
    observations = (1..3).map do |k|
      Observation.new(k * 10.0, Geokit::LatLng.new(0.0, 0.002 - 0.0004 * k))
    end
    candidates = observations.map { [road] }

    result = matcher.match(observations, candidates)

    # Each observation stays snapped to the road, but the matcher must
    # not invent a connected route against the one-way direction.
    assert_equal 3, result.segments.length
    assert_equal [1, 1, 1], result.segments.map { |seg| seg.points.length }
    assert_equal [:unreachable, :unreachable], result.breaks.map(&:kind)
    assert_equal(['one_way'] * 3, result.points.map { |p| p.segment.id })
  end

  def test_duplicate_timestamps_are_rejected
    network, main, = build_connected_network
    matcher = Matcher.new(network: network, sigma: 15.0, beta: 40.0)
    observations = [
      Observation.new(10.0, Geokit::LatLng.new(0.0, 0.0005)),
      Observation.new(10.0, Geokit::LatLng.new(0.0, 0.0010))
    ]
    assert_raise(Geokit::MapMatching::DuplicateTimestampError) do
      matcher.match(observations, [[main], [main]])
    end
  end

  def test_unsorted_observations_are_matched_in_time_order
    network, main, = build_connected_network
    matcher = Matcher.new(network: network, sigma: 15.0, beta: 40.0)
    observations = [
      Observation.new(20.0, Geokit::LatLng.new(0.0, 0.0015)),
      Observation.new(0.0, Geokit::LatLng.new(0.0, 0.0005)),
      Observation.new(10.0, Geokit::LatLng.new(0.0, 0.0010))
    ]
    candidates = observations.map { [main] }
    result = matcher.match(observations, candidates)
    assert_equal [0.0, 10.0, 20.0], result.points.map { |p| p.observation.time }
  end

  def test_mismatched_observations_and_candidates_raise
    network, main, = build_connected_network
    matcher = Matcher.new(network: network, sigma: 15.0, beta: 40.0)
    observations = [Observation.new(0.0, Geokit::LatLng.new(0.0, 0.0))]
    assert_raise(ArgumentError) { matcher.match(observations, [[main], [main]]) }
  end

  def test_invalid_parameters_raise
    network, = build_connected_network
    assert_raise(ArgumentError) do
      Matcher.new(network: network, sigma: 0.0, beta: 40.0)
    end
    assert_raise(ArgumentError) do
      Matcher.new(network: network, sigma: 15.0, beta: -1.0)
    end
    assert_raise(ArgumentError) do
      Matcher.new(network: Object.new, sigma: 15.0, beta: 40.0)
    end
  end

  def test_lambda_network_is_accepted
    network, main, = build_connected_network
    router = ->(from, to) { network.route_distance(from, to) }
    matcher = Matcher.new(network: router, sigma: 15.0, beta: 40.0)
    observation = Observation.new(0.0, Geokit::LatLng.new(0.0, 0.0005))
    result = matcher.match([observation], [[main]])
    assert result.matched?
    assert_equal 'main', result.points.first.segment.id
  end
end

