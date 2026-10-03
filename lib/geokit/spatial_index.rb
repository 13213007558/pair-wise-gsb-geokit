module Geokit
  # In-memory spatial index for stable identifiers on a spherical Earth.
  class SpatialIndex
    RAD_PER_DEGREE = Math::PI / 180.0
    ANGULAR_EPSILON = 1e-12
    DEFAULT_LAT_BINS = 16
    DEFAULT_LNG_BINS = 32

    Cell = Struct.new(:points, :vectors, :min_lat, :max_lat, :min_lng,
                      :max_lng)

    attr_reader :lat_bins, :lng_bins, :distance_computations, :units

    def self.build(points = {}, options = {})
      index = new(options.merge(preallocate: points.size))
      index.bulk_insert(points)
      index
    end

    def initialize(options = {})
      @units = Geokit::LatLng.get_units!(options)
      @earth_radius = Geokit::LatLng.units_sphere_multiplier(@units)
      @distance_computations = 0
      @manual_size = options.key?(:lat_bins) || options.key?(:lng_bins)
      set_bin_sizes(options[:lat_bins] || DEFAULT_LAT_BINS,
                    options[:lng_bins] || DEFAULT_LNG_BINS)
      prepare_containers(options[:preallocate].to_i)
    end

    def size
      @points.size
    end

    def empty?
      @points.empty?
    end

    def [](id)
      @points[id]
    end

    def include?(id)
      @points.key?(id)
    end

    def insert(id, point)
      normalized = normalize_point(point)
      remove_point(id, @points[id]) if @points.key?(id)
      add_point(id, normalized, vector_for(normalized))
      auto_rebuild
      self
    end
    alias_method :upsert, :insert
    alias_method :update, :insert
    alias_method :add, :insert
    alias_method :[]=, :insert

    def bulk_insert(points)
      points.each do |id, point|
        normalized = normalize_point(point)
        remove_point(id, @points[id]) if @points.key?(id)
        add_point(id, normalized, vector_for(normalized))
      end
      rebuild unless @manual_size
      self
    end
    alias_method :bulk_load, :bulk_insert
    alias_method :bulk_add, :bulk_insert

    def delete(id)
      point = @points.delete(id)
      remove_point(id, point) if point
      auto_rebuild
      point
    end
    alias_method :remove, :delete

    def rebuild(lat_bins = nil, lng_bins = nil)
      if lat_bins || lng_bins
        @manual_size = true
        set_bin_sizes(lat_bins || @lat_bins, lng_bins || @lng_bins)
      elsif !@manual_size
        set_bin_sizes_for_count(@points.size)
      end

      saved_points = @points
      saved_vectors = @vectors
      prepare_containers(saved_points.size)
      saved_points.each { |id, point| add_point(id, point, saved_vectors[id]) }
      self
    end

    def within_radius(center, radius, options = {})
      query(center, options) do |normalized_center|
        angular_radius = angular_distance(radius, options)
        center_vector = vector_for(normalized_center)
        results = []

        each_candidate_cell(center_vector, angular_radius) do |cell|
          cell.points.each_key do |id|
            angular = angular_distance_between(center_vector, cell.vectors[id])
            next if angular > angular_radius + radius_epsilon(angular_radius)

            results << [id, angular * @earth_radius]
          end
        end

        format_results(sort_pairs(results), options)
      end
    end
    alias_method :radius_query, :within_radius
    alias_method :radius_search, :within_radius

    def nearest(center, count = 1, options = {})
      raise ArgumentError, 'count must be positive' unless count.is_a?(Integer) && count.positive?

      query(center, options) do |normalized_center|
        center_vector = vector_for(normalized_center)
        cells_by_lower_bound = @cells.reject { |cell| cell.points.empty? }
                               .sort_by { |cell| lower_angular_distance(center_vector, cell) }
        heap = []

        cells_by_lower_bound.each do |cell|
          lower = lower_angular_distance(center_vector, cell)
          cutoff = heap.first&.first
          break if heap.size == count && lower > cutoff + ANGULAR_EPSILON

          cell.points.each_key do |id|
            angular = angular_distance_between(center_vector, cell.vectors[id])
            next if heap.size == count && heap_entry_worse?(angular, id, heap.first)

            heap << [angular, id]
            heap.sort! { |left, right| compare_heap_entries(left, right) }
            heap.shift if heap.size > count
          end
        end

        ordered = heap.map { |angular, id| [id, angular * @earth_radius] }
        format_results(sort_pairs(ordered), options)
      end
    end
    alias_method :k_nearest, :nearest
    alias_method :knn, :nearest

    def reset_distance_computations
      @distance_computations = 0
    end

    private

    def prepare_containers(capacity = 0)
      @points = {}
      @vectors = {}
      @cells = Array.new(@lat_bins * @lng_bins) do
        Cell.new({}, {}, nil, nil, nil, nil)
      end
    end

    def set_bin_sizes(lat_bins, lng_bins)
      @lat_bins = positive_integer(lat_bins, 'lat_bins')
      @lng_bins = positive_integer(lng_bins, 'lng_bins')
    end

    def set_bin_sizes_for_count(count)
      lat_bins = [1, Math.sqrt(count / 8.0).ceil].max
      set_bin_sizes(lat_bins, lat_bins * 2)
    end

    def positive_integer(value, name)
      value = value.to_i unless value.is_a?(Integer)
      raise ArgumentError, "#{name} must be a positive integer" unless value.positive?

      value
    end

    def normalize_point(point)
      latlng = Geokit::LatLng.normalize(point)
      lat = latlng.lat.to_f
      lng = normalize_longitude(latlng.lng.to_f)
      raise ArgumentError, 'coordinates must be finite' unless finite?(lat) && finite?(lng)
      raise ArgumentError, 'latitude must be between -90 and 90' unless lat.between?(-90.0, 90.0)

      Geokit::LatLng.new(lat, lng)
    end

    def normalize_longitude(lng)
      normalized = ((lng + 180.0) % 360.0) - 180.0
      normalized == 180.0 ? -180.0 : normalized
    end

    def finite?(number)
      number.finite? && !number.nan?
    end

    def add_point(id, point, vector)
      @points[id] = point
      @vectors[id] = vector
      cell = cell_for(point)
      cell.points[id] = point
      cell.vectors[id] = vector
      update_cell_bounds(cell, point)
    end

    def remove_point(id, point)
      cell = cell_for(point)
      cell.points.delete(id)
      cell.vectors.delete(id)
      @vectors.delete(id)

      if cell.points.empty?
        reset_cell_bounds(cell)
      else
        rebuild_cell_bounds(cell)
      end
    end

    def cell_for(point)
      lat_index = (((point.lat + 90.0) / 180.0) * @lat_bins).floor
      lat_index = @lat_bins - 1 if lat_index == @lat_bins
      lng_index = (((point.lng + 180.0) / 360.0) * @lng_bins).floor
      lng_index = @lng_bins - 1 if lng_index == @lng_bins
      @cells[(lat_index * @lng_bins) + lng_index]
    end

    def update_cell_bounds(cell, point)
      cell.min_lat = point.lat if cell.min_lat.nil? || point.lat < cell.min_lat
      cell.max_lat = point.lat if cell.max_lat.nil? || point.lat > cell.max_lat
      cell.min_lng = point.lng if cell.min_lng.nil? || point.lng < cell.min_lng
      cell.max_lng = point.lng if cell.max_lng.nil? || point.lng > cell.max_lng
    end

    def rebuild_cell_bounds(cell)
      reset_cell_bounds(cell)
      cell.points.each { |id, point| update_cell_bounds(cell, point) }
    end

    def reset_cell_bounds(cell)
      cell.min_lat = cell.max_lat = nil
      cell.min_lng = cell.max_lng = nil
    end

    def auto_rebuild
      return if @manual_size || @points.empty?

      target_lat_bins = [1, Math.sqrt(@points.size / 8.0).ceil].max
      target_cells = target_lat_bins * target_lat_bins * 2
      rebuild if @cells.size * 4 < target_cells
    end

    def query(center, options)
      old_units = @units
      old_radius = @earth_radius
      @units = options_units(options)
      @earth_radius = sphere_multiplier(@units)
      yield normalize_point(center)
    ensure
      @units = old_units
      @earth_radius = old_radius
    end

    def options_units(options)
      if options.key?(:units)
        Geokit::LatLng.get_units!(options)
      else
        @units
      end
    end

    def sphere_multiplier(units)
      Geokit::LatLng.units_sphere_multiplier(units)
    end

    def angular_distance(distance, options)
      distance = Float(distance)
      raise ArgumentError, 'radius must be finite and non-negative' unless distance.finite? && distance >= 0.0

      distance / sphere_multiplier(options_units(options))
    end

    def radius_epsilon(angular_radius)
      ANGULAR_EPSILON * [1.0, angular_radius.abs].max
    end

    def each_candidate_cell(center_vector, angular_radius)
      @cells.each do |cell|
        next if cell.points.empty?
        next if angular_radius < Math::PI - ANGULAR_EPSILON &&
                lower_angular_distance(center_vector, cell) > angular_radius + radius_epsilon(angular_radius)

        yield cell
      end
    end

    def lower_angular_distance(center_vector, cell)
      lat0 = Math.asin(clamp(center_vector[2], -1.0, 1.0))
      lng0 = Math.atan2(center_vector[1], center_vector[0])
      min_lat = cell.min_lat * RAD_PER_DEGREE
      max_lat = cell.max_lat * RAD_PER_DEGREE
      min_lng = cell.min_lng * RAD_PER_DEGREE
      max_lng = cell.max_lng * RAD_PER_DEGREE
      dots = [
        dot_at(min_lat, min_lng, lat0, lng0),
        dot_at(min_lat, max_lng, lat0, lng0),
        dot_at(max_lat, min_lng, lat0, lng0),
        dot_at(max_lat, max_lng, lat0, lng0),
        longitude_edge_dot(min_lat, min_lng, max_lng, lat0, lng0),
        longitude_edge_dot(max_lat, min_lng, max_lng, lat0, lng0),
        latitude_edge_dot(min_lng, min_lat, max_lat, lat0, lng0),
        latitude_edge_dot(max_lng, min_lat, max_lat, lat0, lng0),
        interior_dot(min_lat, max_lat, min_lng, max_lng, lat0, lng0)
      ]
      Math.acos(clamp(dots.compact.max, -1.0, 1.0))
    end

    def dot_at(lat, lng, lat0, lng0)
      Math.sin(lat0) * Math.sin(lat) +
        Math.cos(lat0) * Math.cos(lat) * Math.cos(lng - lng0)
    end

    def longitude_edge_dot(lat, min_lng, max_lng, lat0, lng0)
      cos_lat = Math.cos(lat)
      return Math.sin(lat0) * Math.sin(lat) if cos_lat.abs < 1e-15

      phase = wrap_angle(lng0)
      candidate_phase(phase, min_lng, max_lng).map do |lng|
        dot_at(lat, lng, lat0, lng0)
      end.max
    end

    def latitude_edge_dot(lng, min_lat, max_lat, lat0, lng0)
      a = Math.cos(lat0) * Math.cos(lng - lng0)
      b = Math.sin(lat0)
      cos_plus_sin_max(a, b, min_lat, max_lat)
    end

    def interior_dot(min_lat, max_lat, min_lng, max_lng, lat0, lng0)
      lat = lat0
      return nil unless lat.between?(min_lat, max_lat)

      lng = lng0
      return nil unless lng.between?(min_lng, max_lng)

      dot_at(lat, lng, lat0, lng0)
    end

    def candidate_phase(phase, minimum, maximum)
      candidates = [minimum, maximum]
      if phase.between?(minimum, maximum)
        candidates << phase
      else
        shifted = phase > maximum ? phase - (2.0 * Math::PI) : phase + (2.0 * Math::PI)
        candidates << shifted if shifted.between?(minimum, maximum)
      end
      candidates
    end

    def cos_plus_sin_max(a, b, minimum, maximum)
      values = [a * Math.cos(minimum) + b * Math.sin(minimum),
                a * Math.cos(maximum) + b * Math.sin(maximum)]
      stationary = wrap_angle(Math.atan2(b, a))
      [stationary, stationary - Math::PI, stationary + Math::PI].each do |candidate|
        next unless candidate.between?(minimum, maximum)

        values << a * Math.cos(candidate) + b * Math.sin(candidate)
      end
      values.max
    end

    def wrap_angle(angle)
      (angle + Math::PI) % (2.0 * Math::PI) - Math::PI
    end

    def angular_distance_between(center, vector)
      @distance_computations += 1
      dot = center[0] * vector[0] + center[1] * vector[1] + center[2] * vector[2]
      Math.acos(clamp(dot, -1.0, 1.0))
    end

    def clamp(value, minimum, maximum)
      return minimum if value < minimum
      return maximum if value > maximum

      value
    end

    def vector_for(point)
      lat = point.lat * RAD_PER_DEGREE
      lng = point.lng * RAD_PER_DEGREE
      cos_lat = Math.cos(lat)
      [cos_lat * Math.cos(lng), cos_lat * Math.sin(lng), Math.sin(lat)]
    end

    def format_results(results, options)
      options[:with_distances] ? results : results.map(&:first)
    end

    def sort_pairs(pairs)
      pairs.sort do |left, right|
        distance = left.last - right.last
        tolerance = ANGULAR_EPSILON * @earth_radius
        next compare_ids(left.first, right.first) if distance.abs <= tolerance

        distance <=> 0
      end
    end

    def heap_entry_worse?(angular, id, entry)
      angular > entry.first + ANGULAR_EPSILON ||
        (angular >= entry.first - ANGULAR_EPSILON && compare_ids(id, entry.last) > 0)
    end

    def compare_heap_entries(left, right)
      distance_comparison = right.first <=> left.first
      return distance_comparison unless distance_comparison.zero?

      compare_ids(right.last, left.last)
    end

    def compare_ids(left, right)
      class_comparison = left.class.name.to_s <=> right.class.name.to_s
      return class_comparison unless class_comparison.zero?

      begin
        left <=> right
      rescue ArgumentError
        left.object_id <=> right.object_id
      end
    end

  end
end
