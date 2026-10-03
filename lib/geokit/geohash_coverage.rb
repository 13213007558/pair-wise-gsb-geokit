# frozen_string_literal: true

module Geokit
  # Build antichain sets of Geohash prefixes covering spherical polygons.
  module GeohashCoverage
    ALPHABET = '0123456789bcdefghjkmnpqrstuvwxyz'.freeze
    MAX_PRECISION = 12
    DEFAULT_SEARCH_BUDGET = 20_000
    EARTH_RADIUS_IN_METERS =
      Geokit::Mappable::ClassMethods::EARTH_RADIUS_IN_METERS

    def self.deg2rad(value)
      value.to_f * Math::PI / 180.0
    end

    def self.rad2deg(value)
      value.to_f * 180.0 / Math::PI
    end

    Cover = Struct.new(
      :prefixes,
      :cell_count,
      :cell_area,
      :polygon_area,
      :overcoverage_area,
      :optimal,
      :cell_budget_reached,
      :cell_budget_limited,
      :max_precision_reached,
      :search_budget_reached,
      :nodes_evaluated
    ) do
      alias optimal? optimal
      alias cell_budget_reached? cell_budget_reached
      alias cell_budget_limited? cell_budget_limited
      alias max_precision_reached? max_precision_reached
      alias search_budget_reached? search_budget_reached
    end

    def self.plan(geometry, max_precision:, max_cells:, holes: [],
                  search_budget: DEFAULT_SEARCH_BUDGET)
      Planner.new(
        geometry,
        holes: holes,
        max_precision: max_precision,
        max_cells: max_cells,
        search_budget: search_budget
      ).call
    end

    def self.all_covers(geometry, max_precision:, holes: [])
      exterior = Ring.new(geometry_for(geometry))
      hole_rings = holes.map { |ring| Ring.new(geometry_for(ring)) }
      analyzer = Analyzer.new(exterior, hole_rings, nil)
      root = Node.root
      covers_for(root, max_precision, analyzer)
    end

    def self.bounds_for_prefix(prefix)
      Node.new(prefix).to_h_bounds
    end

    def self.geometry_for(geometry)
      return geometry.points if geometry.is_a?(Geokit::Polygon)

      geometry
    end

    private_class_method def self.covers_for(node, max_precision, analyzer)
      state = analyzer.state_for(node)
      return [[]] if state == :outside

      take = [[node.prefix]]
      return take if node.depth >= max_precision

      child_covers = node.children.map do |child|
        covers_for(child, max_precision, analyzer)
      end
      splits = child_covers.shift.product(*child_covers).map do |parts|
        parts.flatten.reject(&:empty?)
      end
      take + splits
    end

    module Geometry
      module_function

      def vector(lat, lng)
        lat_r = GeohashCoverage.deg2rad(lat.to_f)
        lng_r = GeohashCoverage.deg2rad(normalize_longitude(lng.to_f))
        cos_lat = Math.cos(lat_r)
        [
          cos_lat * Math.cos(lng_r),
          cos_lat * Math.sin(lng_r),
          Math.sin(lat_r)
        ]
      end

      def normalize_longitude(lng)
        ((lng + 180.0) % 360.0) - 180.0
      end

      def dot(first, second)
        first[0] * second[0] + first[1] * second[1] + first[2] * second[2]
      end

      def cross(first, second)
        [
          first[1] * second[2] - first[2] * second[1],
          first[2] * second[0] - first[0] * second[2],
          first[0] * second[1] - first[1] * second[0]
        ]
      end

      def angle(first, second)
        Math.acos(clamp(dot(first, second), -1.0, 1.0))
      end

      def clamp(value, low, high)
        [low, [value, high].min].max
      end

      def point_on_arc?(point, start, finish, epsilon = 1e-7)
        edge_angle = angle(start, finish)
        return angle(point, start) <= epsilon if edge_angle < 1e-10

        triple = dot(point, cross(start, finish)).abs
        along = angle(start, point) + angle(point, finish) - edge_angle
        triple <= 1e-8 && along.abs <= epsilon
      end

      def point_in_ring(point, ring)
        attempts = [point, jitter(point, 1e-8), jitter(point, -1.3e-8),
                    jitter(point, 2.7e-8)]
        attempts.each do |candidate|
          result = winding_result(candidate, ring)
          return result unless result == :degenerate
        end
        :unknown
      end

      def winding_result(point, ring)
        ring.edges.each do |start, finish|
          return :boundary if point_on_arc?(point, start, finish)
        end

        winding = 0.0
        ring.edges.each do |start, finish|
          start_tangent = tangent(point, start)
          finish_tangent = tangent(point, finish)
          return :degenerate if start_tangent.nil? || finish_tangent.nil?

          winding += Math.atan2(
            dot(point, cross(start_tangent, finish_tangent)),
            dot(start_tangent, finish_tangent)
          )
        end
        winding.abs > Math::PI ? :inside : :outside
      end

      def tangent(origin, point)
        projection = dot(origin, point)
        vector = point.each_with_index.map { |value, index|
          value - origin[index] * projection
        }
        length = Math.sqrt(dot(vector, vector))
        return if length < 1e-12

        vector.map { |value| value / length }
      end

      def jitter(point, amount)
        lat = Math.asin(clamp(point[2] + amount, -1.0, 1.0))
        lng = Math.atan2(point[1], point[0]) + amount
        [
          Math.cos(lat) * Math.cos(lng),
          Math.cos(lat) * Math.sin(lng),
          Math.sin(lat)
        ]
      end

      def ring_area(ring)
        reference_points.each do |reference|
          next if ring.edges.any? { |start, finish|
            point_on_arc?(reference, start, finish, 1e-9)
          }

          area = signed_ring_area(ring, reference).abs
          return area if area > 1e-12
        end

        raise ArgumentError, 'unable to calculate spherical polygon area'
      end

      def signed_ring_area(ring, reference)
        ring.edges.inject(0.0) do |total, (start, finish)|
          denominator = 1.0 + dot(reference, start) + dot(start, finish) +
                        dot(finish, reference)
          return 0.0 if denominator.abs < 1e-11

          numerator = dot(reference, cross(start, finish))
          total + 2.0 * Math.atan2(numerator, denominator)
        end
      end

      def reference_points
        [
          vector(0, 0),
          vector(37, -51),
          vector(-63, 127),
          vector(82, 73),
          vector(-82, -17)
        ]
      end
    end

    class Ring
      attr_reader :lat_lngs, :vectors

      def initialize(points)
        @lat_lngs = normalize_points(points)
        @vectors = @lat_lngs.map do |point|
          Geometry.vector(point.lat, point.lng)
        end
        validate!
      end

      def edges
        @edges ||= vectors[0...-1].zip(vectors[1...vectors.length])
      end

      private

      def normalize_points(points)
        raise ArgumentError, 'a polygon ring requires points' if points.nil?
        normalized = points.map { |point| Geokit::LatLng.normalize(point) }
        normalized << normalized.first unless normalized.empty? ||
          normalized.first == normalized.last
        normalized
      end

      def validate!
        if lat_lngs.length < 4
          raise ArgumentError, 'a polygon ring requires at least 3 points'
        end
        lat_lngs.each do |point|
          unless point.lat.finite? && point.lng.finite?
            raise ArgumentError, 'polygon coordinates must be finite'
          end
          unless point.lat.between?(-90.0, 90.0)
            raise ArgumentError, 'latitude must be between -90 and 90'
          end
        end
        edges.each do |start, finish|
          if (Geometry.angle(start, finish) - Math::PI).abs < 1e-9
            raise ArgumentError,
                  'antipodal adjacent polygon points require an '
                  'intermediate point'
          end
        end
      end
    end

    class Node
      attr_reader :prefix, :south, :north, :west, :east, :depth, :area, :state
      attr_writer :state

      def self.root
        new('')
      end

      def initialize(prefix)
        @prefix = prefix
        @depth = prefix.length
        @state = nil
        @south, @north, @west, @east = bounds(prefix)
        @area = calculate_area
      end

      def children
        ALPHABET.each_char.map { |character| Node.new(prefix + character) }
      end

      def to_h_bounds
        {
          south: south,
          north: north,
          west: west,
          east: east
        }
      end

      private

      def bounds(prefix)
        south = -90.0
        north = 90.0
        west = -180.0
        east = 180.0

        prefix.each_char do |character|
          value = ALPHABET.index(character)
          5.times do |bit|
            bit_value = value & (1 << (4 - bit))
            if bit.even?
              middle = (west + east) / 2.0
              if bit_value != 0
                west = middle
              else
                east = middle
              end
            else
              middle = (south + north) / 2.0
              if bit_value != 0
                south = middle
              else
                north = middle
              end
            end
          end
        end

        [south, north, west, east]
      end

      def calculate_area
        longitude_span = GeohashCoverage.deg2rad(east - west)
        latitude_factor = Math.sin(GeohashCoverage.deg2rad(north)) -
                          Math.sin(GeohashCoverage.deg2rad(south))
        longitude_span * latitude_factor * EARTH_RADIUS_IN_METERS**2
      end
    end

    class Analyzer
      attr_reader :nodes_evaluated
      attr_writer :truncated

      def initialize(exterior, holes, search_budget)
        @rings = [exterior] + holes
        @search_budget = search_budget
        @remaining = search_budget
        @nodes_evaluated = 0
        @truncated = false
        @states = {}
      end

      def truncated?
        @truncated
      end

      def state_for(node)
        return @states[node.prefix] if @states.key?(node.prefix)

        if @remaining && @remaining <= 0
          @truncated = true
          return @states[node.prefix] = :unknown
        end

        @remaining -= 1 if @remaining
        @nodes_evaluated += 1
        center = Geometry.vector((node.south + node.north) / 2.0,
                                 (node.west + node.east) / 2.0)
        inside = polygon_contains?(center)
        intersected = @rings.any? { |ring| ring_intersects?(ring, node) }
        state =
          if inside == :inside || inside == :boundary || inside == :unknown
            intersected || inside != :inside ? :mixed : :inside
          else
            intersected ? :mixed : :outside
          end
        node.state = state
        @states[node.prefix] = state
      end

      private

      def polygon_contains?(point)
        exterior = @rings.first
        unless Geometry.point_in_ring(point, exterior) == :inside
          return :outside
        end

        @rings[1..-1].each do |hole|
          result = Geometry.point_in_ring(point, hole)
          return :boundary if result == :boundary
          return :outside if result == :inside
          return :unknown if result == :unknown
        end
        :inside
      end

      def ring_intersects?(ring, node)
        ring.edges.any? { |start, finish|
          edge_intersects_node?(start, finish, node)
        }
      end

      def edge_intersects_node?(start, finish, node)
        omega = Geometry.angle(start, finish)
        return false if omega < 1e-10

        times = [0.0, 1.0]
        times += latitude_roots(start, finish, omega, node)
        unless node.prefix.empty?
          times += longitude_roots(start, finish, omega, node.west)
          times += longitude_roots(start, finish, omega, node.east)
        end
        times = normalize_times(times)
        points = times.map { |time| arc_point(start, finish, omega, time) }
        times.each_cons(2) do |first, second|
          points << arc_point(start, finish, omega, (first + second) / 2.0)
        end
        points.any? { |point| point_in_node_interior?(point, node) }
      end

      def latitude_roots(start, finish, omega, node)
        [node.south, node.north].flat_map do |latitude|
          target = Math.sin(GeohashCoverage.deg2rad(latitude))
          roots(start[2], finish[2], target, omega).map { |value|
            value / omega
          }
        end
      end

      def longitude_roots(start, finish, omega, longitude)
        longitude_r = GeohashCoverage.deg2rad(longitude)
        sin_longitude = Math.sin(longitude_r)
        cos_longitude = Math.cos(longitude_r)
        coefficient_start = start[0] * sin_longitude - start[1] * cos_longitude
        coefficient_finish = finish[0] * sin_longitude -
                             finish[1] * cos_longitude
        roots(coefficient_start, coefficient_finish, 0.0, omega)
          .map { |value| value / omega }
      end

      def roots(start_value, finish_value, target, omega)
        cos_omega = Math.cos(omega)
        sin_omega = Math.sin(omega)
        sine_coefficient =
          (finish_value - start_value * cos_omega) / sin_omega
        radius = Math.sqrt(
          sine_coefficient**2 + start_value**2
        )
        return [] if radius < 1e-14

        ratio = Geometry.clamp(target / radius, -1.0, 1.0)
        phase = Math.atan2(start_value, sine_coefficient)
        alpha = Math.asin(ratio)
        [alpha - phase, Math::PI - alpha - phase].select { |value|
          value.between?(-1e-10, omega + 1e-10)
        }.map { |value| Geometry.clamp(value, 0.0, omega) }
      end

      def arc_point(start, finish, omega, time)
        angle = omega * time
        sin_angle = Math.sin(angle)
        cos_angle = Math.cos(angle)
        sin_omega = Math.sin(omega)
        vector = [
          (Math.sin(omega - angle) * start[0] + sin_angle * finish[0]) /
            sin_omega,
          (Math.sin(omega - angle) * start[1] + sin_angle * finish[1]) /
            sin_omega,
          (Math.sin(omega - angle) * start[2] + sin_angle * finish[2]) /
            sin_omega
        ]
        length = Math.sqrt(Geometry.dot(vector, vector))
        vector.map { |value| value / length }
      end

      def normalize_times(times)
        times.sort!
        values = []
        times.each do |time|
          time = Geometry.clamp(time, 0.0, 1.0)
          values << time if values.empty? || (values.last - time).abs > 1e-9
        end
        values
      end

      def point_in_node_interior?(point, node)
        latitude_rad = Math.asin(Geometry.clamp(point[2], -1.0, 1.0))
        latitude = GeohashCoverage.rad2deg(latitude_rad)
        return false unless latitude > node.south + 1e-10 &&
                            latitude < node.north - 1e-10
        return true if node.prefix.empty?

        equatorial = Math.sqrt(point[0]**2 + point[1]**2)
        return false if equatorial < 1e-12

        longitude = GeohashCoverage.rad2deg(Math.atan2(point[1], point[0]))
        longitude > node.west + 1e-10 && longitude < node.east - 1e-10
      end
    end

    class Planner
      Option = Struct.new(:count, :area, :exact, :refinable, :selection)

      def initialize(geometry, holes:, max_precision:, max_cells:,
                     search_budget:)
        @exterior = Ring.new(GeohashCoverage.geometry_for(geometry))
        @holes = holes.map { |ring|
          Ring.new(GeohashCoverage.geometry_for(ring))
        }
        @max_precision = validate_precision(max_precision)
        @max_cells = validate_positive_integer(max_cells, :max_cells)
        @search_budget = validate_positive_integer(
          search_budget, :search_budget
        )
        @analyzer = Analyzer.new(@exterior, @holes, @search_budget)
        @nodes = {}
      end

      def call
        root = node('')
        option = best_option(root, @max_cells)
        polygon_area = calculate_polygon_area
        cell_area = option.selection.sum(&:area)
        prefixes = option.selection.map(&:prefix).sort
        overcoverage = [cell_area - polygon_area, 0.0].max
        cell_count = option.selection.length
        precision_reached = option.selection.any? do |selected|
          selected.depth == @max_precision && selected.state == :mixed
        end

        Cover.new(
          prefixes,
          cell_count,
          cell_area,
          polygon_area,
          overcoverage,
          !@analyzer.truncated? && option.exact,
          cell_count == @max_cells,
          cell_count == @max_cells && option.refinable,
          precision_reached,
          @analyzer.truncated?,
          @analyzer.nodes_evaluated
        )
      end

      private

      def best_option(node, budget)
        options = options_for(node, budget)
        options.min_by(&:area)
      end

      def options_for(node, budget)
        state = @analyzer.state_for(node)
        return [empty_option] if state == :outside

        if state == :inside || node.depth == @max_precision
          return [take_option(node, state)]
        end
        if state == :unknown
          return [take_option(node, :unknown)]
        end

        options = [take_option(node, state)]
        split = split_options(node, budget)
        prune(options + split)
      end

      def split_options(node, budget)
        current = [empty_option]
        GeohashCoverage::ALPHABET.each_char do |character|
          child = node(node.prefix + character)
          child_options = options_for(child, budget)
          combined = []
          current.each do |first|
            child_options.each do |second|
              next if first.count + second.count > budget

              combined << Option.new(
                first.count + second.count,
                first.area + second.area,
                first.exact && second.exact,
                first.refinable || second.refinable,
                first.selection + second.selection
              )
            end
          end
          current = prune(combined)
          if current.empty?
            return [take_option(node, @analyzer.state_for(node))]
          end
        end
        current
      end

      def take_option(node, state)
        Option.new(
          1,
          node.area,
          state != :unknown,
          state == :mixed && node.depth < @max_precision,
          [node]
        )
      end

      def empty_option
        Option.new(0, 0.0, true, false, [])
      end

      def prune(options)
        best = {}
        options.each do |option|
          existing = best[option.count]
          if existing.nil? || option.area < existing.area
            best[option.count] = option
          end
        end

        retained = []
        best.keys.sort.each do |count|
          option = best[count]
          if retained.empty? || option.area < retained.last.area - 1e-9
            retained << option
          end
        end
        retained
      end

      def node(prefix)
        @nodes[prefix] ||= Node.new(prefix)
      end

      def calculate_polygon_area
        exterior_area = Geometry.ring_area(@exterior)
        hole_area = @holes.sum { |hole| Geometry.ring_area(hole) }
        (exterior_area - hole_area) * EARTH_RADIUS_IN_METERS**2
      end

      def validate_precision(value)
        precision = validate_positive_integer(value, :max_precision)
        unless precision <= MAX_PRECISION
          raise ArgumentError, "max_precision must be <= #{MAX_PRECISION}"
        end

        precision
      end

      def validate_positive_integer(value, name)
        integer =
          begin
            Integer(value)
          rescue ArgumentError, TypeError
            raise ArgumentError, "#{name} must be an integer"
          end
        minimum = name == :max_precision ? 0 : 1
        if integer < minimum
          raise ArgumentError,
                "#{name} must be greater than or equal to #{minimum}"
        end
        integer
      end
    end
  end
end
