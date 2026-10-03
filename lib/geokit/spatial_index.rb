module Geokit
  # In-memory spatial index for points on a sphere. Point IDs are stable and
  # several IDs may share the same coordinate.
  class SpatialIndex
    # A point returned by a spatial query.
    Result = Struct.new(:id, :location, :distance) do
      alias_method :latlng, :location
      alias_method :point, :location
    end

    class DuplicateIdError < ArgumentError; end
    class MissingIdError < KeyError; end

    Node = Struct.new(
      :vector, :location, :ids, :left, :right, :min, :max, :axis
    )

    attr_reader :size, :distance_calculation_count

    # Accepts an optional enumerable of +[id, point]+ pairs, or a Hash mapping
    # IDs to points. Points may be LatLng objects or [lat, lng] arrays.
    def initialize(entries = nil, options = {})
      @units = options[:units] || Geokit.default_units
      clear
      bulk_load(entries) if entries
    end

    def add(id, point_or_lat, lng = nil)
      vector, location = prepare_point(point_or_lat, lng)
      old_vector = @locations[id]

      if old_vector
        return self if old_vector == vector
        remove_point(id, old_vector)
      else
        @size += 1
      end

      node = @point_nodes[vector]
      if node
        node.ids << id
      else
        node = Node.new(vector, location, [id], nil, nil, nil, nil, 0)
        @root = insert_node(@root, node, 0)
      end

      @locations[id] = vector
      self
    end
    alias_method :insert, :add
    alias_method :upsert, :add

    def remove(id)
      vector = @locations.delete(id)
      return false unless vector

      remove_point(id, vector)
      @size -= 1
      true
    end
    alias_method :delete, :remove

    def remove!(id)
      remove(id) or raise MissingIdError, "index does not contain #{id.inspect}"
    end

    def include?(id)
      @locations.key?(id)
    end

    def location_for(id)
      vector = @locations.fetch(id) do
        raise MissingIdError, "index does not contain #{id.inspect}"
      end
      @point_nodes[vector].location
    end

    def bulk_load(entries)
      clear
      groups = {}

      each_entry(entries) do |id, point_or_lat, lng|
        if @locations.key?(id)
          raise DuplicateIdError, "duplicate point ID #{id.inspect}"
        end

        vector, location = prepare_point(point_or_lat, lng)
        @locations[id] = vector
        group = groups[vector]
        if group
          group.ids << id
        else
          groups[vector] = Node.new(vector, location, [id], nil, nil,
                                   nil, nil, 0)
        end
      end

      @root = build_tree(groups.values, 0)
      @size = @locations.size
      self
    end

    def clear
      @root = nil
      @locations = {}
      @point_nodes = {}
      @size = 0
      @distance_calculation_count = 0
      self
    end

    # Returns every point within +radius+ of +origin+, ordered by distance and
    # then by ID. Radius and returned distances use +units+ (defaulting to the
    # index's units).
    def within_radius(origin, radius, options = {})
      units = get_units(options)
      raise ArgumentError, 'radius must be non-negative' unless valid_radius?(radius)

      _, origin_location = prepare_point(origin)
      radius_angle = radius.to_f / earth_radius(units)
      results = []
      @distance_calculation_count = 0
      radius_search(@root, vector_for(origin_location), origin_location,
                    radius_angle, radius.to_f, units, results)
      sort_results(results)
    end
    alias_method :radius_query, :within_radius

    # Returns up to +k+ nearest points, ordered by distance and then by ID.
    # Ties are never resolved arbitrarily.
    def nearest(origin, k = 1, options = {})
      units = get_units(options)
      raise ArgumentError, 'k must be a positive Integer' unless valid_k?(k)

      _, origin_location = prepare_point(origin)
      best = []
      @distance_calculation_count = 0
      nearest_search(@root, vector_for(origin_location), origin_location,
                     k, units, best)
      sort_results(best)
    end
    alias_method :knn, :nearest
    alias_method :k_nearest, :nearest

    private

    # The slack is deliberately conservative. It can make the index examine a
    # few additional points, but never allows an approximate or missed answer.
    PRUNING_EPSILON = 1e-12

    def prepare_point(point_or_lat, lng = nil)
      location = if lng
                   Geokit::LatLng.new(point_or_lat, lng)
                 else
                   Geokit::LatLng.normalize(point_or_lat)
                 end

      lat = location.lat.to_f
      lng_value = location.lng.to_f
      unless finite?(lat) && finite?(lng_value)
        raise ArgumentError, 'latitude and longitude must be finite'
      end
      unless lat >= -90.0 && lat <= 90.0
        raise ArgumentError, 'latitude must be between -90 and 90'
      end

      lng_value = 0.0 if lat.abs == 90.0
      lng_value = ((lng_value + 180.0) % 360.0) - 180.0
      location = Geokit::LatLng.new(lat, lng_value)
      [vector_for(location), location]
    end

    def vector_for(location)
      lat = deg2rad(location.lat)
      lng = deg2rad(location.lng)
      cos_lat = Math.cos(lat)
      [
        cos_lat * Math.cos(lng),
        cos_lat * Math.sin(lng),
        Math.sin(lat)
      ]
    end

    def build_tree(nodes, depth)
      return nil if nodes.empty?

      axis = depth % 3
      nodes.sort! { |left, right| left.vector[axis] <=> right.vector[axis] }
      middle = nodes.size / 2
      node = nodes[middle]
      node.axis = axis
      node.left = build_tree(nodes[0...middle], depth + 1)
      node.right = build_tree(nodes[(middle + 1)..-1], depth + 1)
      refresh_bounds(node)
      @point_nodes[node.vector] = node
      node
    end

    def insert_node(parent, node, depth)
      unless parent
        node.axis = depth % 3
        node.min = node.vector.dup
        node.max = node.vector.dup
        @point_nodes[node.vector] = node
        return node
      end

      axis = parent.axis
      if node.vector[axis] < parent.vector[axis]
        parent.left = insert_node(parent.left, node, depth + 1)
      else
        parent.right = insert_node(parent.right, node, depth + 1)
      end
      refresh_bounds(parent)
      parent
    end

    def remove_point(id, vector)
      node = @point_nodes[vector]
      node.ids.delete(id)
      return if node.ids.any?

      @point_nodes.delete(vector)
      @root = remove_node(@root, vector, 0)
    end

    def remove_node(node, vector, depth)
      return nil unless node

      axis = depth % 3
      if node.vector == vector
        if node.left.nil? && node.right.nil?
          return nil
        elsif node.right
          successor = minimum_node(node.right, axis, depth + 1)
          replace_node(node, successor)
          node.right = remove_node(
            node.right, successor.vector, depth + 1
          )
        else
          predecessor = maximum_node(node.left, axis, depth + 1)
          replace_node(node, predecessor)
          node.left = remove_node(
            node.left, predecessor.vector, depth + 1
          )
        end
        refresh_bounds(node)
        return node
      end

      if vector[axis] < node.vector[axis]
        node.left = remove_node(node.left, vector, depth + 1)
      elsif vector[axis] > node.vector[axis]
        node.right = remove_node(node.right, vector, depth + 1)
      elsif contains_point?(node.left, vector)
        node.left = remove_node(node.left, vector, depth + 1)
      else
        node.right = remove_node(node.right, vector, depth + 1)
      end

      refresh_bounds(node)
      node
    end

    def replace_node(node, replacement)
      node.vector = replacement.vector
      node.location = replacement.location
      node.ids = replacement.ids
      @point_nodes[node.vector] = node
    end

    def minimum_node(node, axis, depth)
      return node unless node
      return node if depth % 3 == axis && node.left.nil?

      candidates = [node]
      candidates << minimum_node(node.left, axis, depth + 1) if node.left
      candidates << minimum_node(node.right, axis, depth + 1) if node.right
      candidates.compact.min_by { |candidate| candidate.vector[axis] }
    end

    def maximum_node(node, axis, depth)
      return node unless node
      return node if depth % 3 == axis && node.right.nil?

      candidates = [node]
      candidates << maximum_node(node.left, axis, depth + 1) if node.left
      candidates << maximum_node(node.right, axis, depth + 1) if node.right
      candidates.compact.max_by { |candidate| candidate.vector[axis] }
    end

    def contains_point?(node, vector)
      return false unless node
      return true if node.vector == vector
      contains_point?(node.left, vector) || contains_point?(node.right, vector)
    end

    def refresh_bounds(node)
      node.min = node.vector.dup
      node.max = node.vector.dup
      [node.left, node.right].each do |child|
        next unless child
        3.times do |axis|
          node.min[axis] = [node.min[axis], child.min[axis]].min
          node.max[axis] = [node.max[axis], child.max[axis]].max
        end
      end
      node
    end

    def radius_search(
      node, origin_vector, origin, radius_angle, radius_distance, units, results
    )
      return unless node
      return if prunable?(node, origin_vector, radius_angle)

      distance = sphere_distance(origin, node.location, units)
      if distance <= radius_distance
        node.ids.each { |id| results << Result.new(id, node.location, distance) }
      end

      radius_search(
        node.left, origin_vector, origin, radius_angle,
        radius_distance, units, results
      )
      radius_search(
        node.right, origin_vector, origin, radius_angle,
        radius_distance, units, results
      )
    end

    def nearest_search(node, origin_vector, origin, k, units, best)
      return unless node
      cutoff = best.size == k ? best.last.distance / earth_radius(units) : nil
      return if prunable?(node, origin_vector, cutoff)

      distance = sphere_distance(origin, node.location, units)
      node.ids.each do |id|
        add_candidate(Result.new(id, node.location, distance), k, best)
      end

      children = [node.left, node.right].compact
      children.sort_by! { |child| bound_score(child, origin_vector) }
      children.each do |child|
        nearest_search(child, origin_vector, origin, k, units, best)
      end
    end

    def add_candidate(candidate, k, best)
      if best.size < k
        insert_candidate(candidate, best)
      elsif candidate_before?(candidate, best.last)
        best.pop
        insert_candidate(candidate, best)
      end
    end

    def insert_candidate(candidate, best)
      index = best.index { |existing| candidate_before?(candidate, existing) }
      if index
        best.insert(index, candidate)
      else
        best << candidate
      end
    end

    def candidate_before?(left, right)
      left.distance < right.distance ||
        (left.distance == right.distance && left.id < right.id)
    end

    def prunable?(node, origin_vector, max_angle)
      return false if max_angle.nil? || max_angle >= Math::PI

      # The existing haversine-style implementation returns zero when acos
      # receives a value just outside [-1, 1]. Keep such branches exact.
      return false if numerical_zero_distance_possible?(node, origin_vector)

      distance_squared = bound_distance_squared(node, origin_vector)
      chord_limit = 2.0 * Math.sin(max_angle / 2.0)
      distance_squared > (chord_limit * chord_limit) + PRUNING_EPSILON
    end

    def numerical_zero_distance_possible?(node, origin_vector)
      minimum_dot = 0.0
      maximum_dot = 0.0
      3.times do |axis|
        coordinate = origin_vector[axis]
        if coordinate >= 0.0
          minimum_dot += coordinate * node.min[axis]
          maximum_dot += coordinate * node.max[axis]
        else
          minimum_dot += coordinate * node.max[axis]
          maximum_dot += coordinate * node.min[axis]
        end
      end
      minimum_dot < -1.0 + PRUNING_EPSILON ||
        maximum_dot > 1.0 - PRUNING_EPSILON
    end

    def bound_score(node, origin_vector)
      return 0.0 if numerical_zero_distance_possible?(node, origin_vector)
      bound_distance_squared(node, origin_vector)
    end

    def bound_distance_squared(node, origin_vector)
      distance_squared = 0.0
      3.times do |axis|
        closest = if origin_vector[axis] < node.min[axis]
                    node.min[axis]
                  elsif origin_vector[axis] > node.max[axis]
                    node.max[axis]
                  else
                    origin_vector[axis]
                  end
        difference = closest - origin_vector[axis]
        distance_squared += difference * difference
      end
      distance_squared
    end

    def sphere_distance(origin, location, units)
      @distance_calculation_count += 1
      Geokit::LatLng.distance_between(
        origin, location, units: units, formula: :sphere
      )
    end

    def sort_results(results)
      results.sort_by! { |result| [result.distance, result.id] }
    end

    def each_entry(entries)
      case entries
      when Hash
        entries.each { |id, point| yield id, point, nil }
      when Enumerable
        entries.each do |entry|
          unless entry.is_a?(Array) && entry.size >= 2
            raise ArgumentError, 'entries must be [id, point] pairs'
          end
          id, point, lng = entry
          yield id, point, lng
        end
      else
        raise ArgumentError, 'entries must be a Hash or Enumerable'
      end
    end

    def valid_radius?(radius)
      radius.is_a?(Numeric) && finite?(radius) && radius >= 0
    end

    def valid_k?(k)
      k.is_a?(Integer) && k > 0
    end

    def finite?(value)
      value.is_a?(Numeric) && value.finite?
    end

    def get_units(options)
      units = options[:units] || @units
      Geokit::LatLng.send(:get_units!, options.merge(units: units))
    end

    def earth_radius(units)
      Geokit::LatLng.units_sphere_multiplier(units)
    end

    def deg2rad(degrees)
      Geokit::LatLng.deg2rad(degrees)
    end
  end
end
