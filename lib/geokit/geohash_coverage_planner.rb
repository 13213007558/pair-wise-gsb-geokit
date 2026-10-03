require 'geokit/geohash'
require 'geokit/spherical_polygon'

module Geokit
  # Builds an antichain of Geohash prefixes covering a spherical polygon.
  class GeohashCoveragePlanner
    DEFAULT_SEARCH_BUDGET = 100_000
    EARTH_RADIUS = 6_376_772.71

    Covering = Struct.new(:prefixes, :cell_count, :covered_area,
                          :overcoverage_area, :polygon_area,
                          :minimum_cell_count_at_max_precision,
                          :cell_budget_limited, :search_budget,
                          :nodes_searched, :search_budget_exhausted,
                          :optimal, keyword_init: true) do
      alias_method :to_a, :prefixes

      def optimal?
        optimal
      end

      def cell_budget_limited?
        cell_budget_limited
      end

      def search_budget_exhausted?
        search_budget_exhausted
      end

      def antichain?
        prefixes.none? do |prefix|
          prefixes.any? { |other| other != prefix && other.start_with?(prefix) }
        end
      end
    end

    class << self
      def plan(polygon, max_precision:, max_cells:, search_budget: DEFAULT_SEARCH_BUDGET,
               earth_radius: EARTH_RADIUS)
        new(max_precision: max_precision, max_cells: max_cells,
            search_budget: search_budget, earth_radius: earth_radius).plan(polygon)
      end

      def exhaustive_plan(polygon, max_precision:, max_cells: nil,
                          earth_radius: EARTH_RADIUS)
        new(max_precision: max_precision, max_cells: max_cells || Float::INFINITY,
            search_budget: nil, earth_radius: earth_radius)
          .exhaustive_plan(polygon)
      end

      def each_covering(polygon, max_precision:, max_cells: nil, &block)
        new(max_precision: max_precision, max_cells: max_cells || Float::INFINITY,
            search_budget: nil).each_covering(polygon, &block)
      end
    end

    def initialize(max_precision:, max_cells:, search_budget: DEFAULT_SEARCH_BUDGET,
                   earth_radius: EARTH_RADIUS)
      Geohash.validate_precision(max_precision)
      @max_precision = Integer(max_precision)
      @max_cells = max_cells.nil? || max_cells.infinite? ? Float::INFINITY : Integer(max_cells)
      raise ArgumentError, 'max_cells must be positive' unless @max_cells.positive?
      raise ArgumentError, 'search_budget must be positive' if search_budget && Integer(search_budget) <= 0

      @search_budget = search_budget.nil? ? Float::INFINITY : Integer(search_budget)
      @earth_radius = earth_radius.to_f
      @nodes_searched = 0
      @search_budget_exhausted = false
    end

    def plan(input_polygon)
      polygon = SphericalPolygon.new(input_polygon, @earth_radius)
      root = build_adaptive_node('', 0, polygon)
      explored_minimum = minimum_cell_count(root)
      choice = best_choice(root, @max_cells)
      raise 'failed to construct a covering' unless choice

      exact_search = !@search_budget_exhausted
      build_covering(polygon, choice, explored_minimum, exact_search, false)
    end

    def exhaustive_plan(input_polygon)
      polygon = SphericalPolygon.new(input_polygon, @earth_radius)
      root = build_exhaustive_node('', 0, polygon)
      choice = best_choice(root, @max_cells)
      raise 'requested cell budget cannot represent a covering' unless choice

      minimum = minimum_cell_count(root)
      build_covering(polygon, choice, minimum, true, true)
    end

    def each_covering(input_polygon, &block)
      return enum_for(__method__, input_polygon) unless block_given?

      polygon = SphericalPolygon.new(input_polygon, @earth_radius)
      root = build_exhaustive_node('', 0, polygon)
      enumerate_node(root, @max_cells) { |prefixes| block.call(prefixes) }
    end

    private

    Node = Struct.new(:prefix, :bounds, :depth, :kind, :cell_area,
                      :intersection_area, :children, keyword_init: true)
    Choice = Struct.new(:count, :cost, :prefixes, keyword_init: true)

    def build_adaptive_node(prefix, depth, polygon)
      node = measure_node(prefix, depth, polygon)
      return node if node.kind != :partial || depth == @max_precision

      if search_remaining < Geohash::BASE32.length
        @search_budget_exhausted = true
        return node
      end

      node.children = build_children(prefix, depth, polygon, :adaptive)
      node
    end

    def build_exhaustive_node(prefix, depth, polygon)
      node = measure_node(prefix, depth, polygon)
      return node if node.kind == :outside || depth == @max_precision

      node.children = build_children(prefix, depth, polygon, :exhaustive)
      node
    end

    def build_children(prefix, depth, polygon, mode)
      Geohash::BASE32.each_char.map do |character|
        child_prefix = prefix + character
        if mode == :adaptive
          build_adaptive_node(child_prefix, depth + 1, polygon)
        else
          build_exhaustive_node(child_prefix, depth + 1, polygon)
        end
      end
    end

    def measure_node(prefix, depth, polygon)
      @nodes_searched += 1
      bounds = Geohash.bounds(prefix)
      area = rectangle_area(bounds)
      intersection = polygon.intersection_area(bounds)
      tolerance = [area * 1e-12, 1e-2].max
      full = area - intersection <= tolerance
      kind = if full
               :inside
             elsif intersection <= tolerance
               :outside
             else
               :partial
             end

      Node.new(prefix: prefix, bounds: bounds, depth: depth, kind: kind,
               cell_area: area, intersection_area: intersection, children: [])
    end

    def rectangle_area(bounds)
      width = (bounds[:east] - bounds[:west]) * Math::PI / 180.0
      @earth_radius**2 * width *
        (Math.sin(bounds[:north] * Math::PI / 180.0) -
         Math.sin(bounds[:south] * Math::PI / 180.0))
    end

    def minimum_cell_count(node)
      case node.kind
      when :outside
        0
      when :inside
        1
      else
        node.children.empty? ? 1 : node.children.sum { |child| minimum_cell_count(child) }
      end
    end

    def best_choice(node, limit)
      choices = node_choices(node, limit)
      choices.min_by { |choice| [choice.cost, choice.count, choice.prefixes] }
    end

    def node_choices(node, limit)
      return [] if node.kind == :outside || limit <= 0

      self_choice = Choice.new(count: 1, cost: overcoverage_cost(node),
                               prefixes: [node.prefix])
      return [self_choice] if node.kind == :inside || node.children.empty?

      combined = combine_children(node.children, limit)
      [self_choice, *combined].select { |choice| choice.count <= limit }
    end

    def combine_children(children, limit)
      combinations = [Choice.new(count: 0, cost: 0.0, prefixes: [])]

      children.each do |child|
        child_choices = node_choices(child, limit)
        child_minimum = minimum_cell_count(child)
        next_start = []
        combinations.each do |left|
          available = limit - left.count
          if child_minimum.zero?
            next_start << left
          else
            child_options = child_choices.select { |choice| choice.count <= available }
            next if child_options.empty?

            child_options.each do |right|
              count = left.count + right.count
              next_start << Choice.new(
                count: count,
                cost: left.cost + right.cost,
                prefixes: (left.prefixes + right.prefixes).sort
              )
            end
          end
        end
        combinations = pareto(next_start)
      end

      combinations
    end

    def pareto(choices)
      choices.group_by(&:count).values.map do |group|
        group.min_by { |choice| [choice.cost, choice.prefixes] }
      end.reject do |choice|
        choices.any? do |other|
          other.count < choice.count && other.cost <= choice.cost &&
            (other.cost < choice.cost || lex_before?(other.prefixes, choice.prefixes))
        end
      end
    end

    def lex_before?(left, right)
      (left <=> right) == -1
    end

    def enumerate_node(node, limit)
      return yield([]) if node.kind == :outside
      return if limit <= 0

      yield([node.prefix])
      return if node.depth == @max_precision || node.children.empty?

      enumerate_children(node.children, 0, limit, []) do |prefixes|
        yield(prefixes)
      end
    end

    def enumerate_children(children, index, remaining, selected)
      return yield(selected.sort) if index == children.length

      enumerate_node(children[index], remaining) do |child_prefixes|
        next if child_prefixes.length > remaining

        enumerate_children(children, index + 1, remaining - child_prefixes.length,
                           selected + child_prefixes) do |prefixes|
          yield(prefixes)
        end
      end
    end

    def build_covering(polygon, choice, minimum, exact_search, exhaustive)
      covered_area = choice.prefixes.sum { |prefix| rectangle_area(Geohash.bounds(prefix)) }
      overcoverage = [covered_area - polygon.area, 0.0].max
      limited = minimum > @max_cells
      search_exhausted = !exhaustive && @search_budget_exhausted

      Covering.new(
        prefixes: choice.prefixes,
        cell_count: choice.count,
        covered_area: covered_area,
        overcoverage_area: overcoverage,
        polygon_area: polygon.area,
        minimum_cell_count_at_max_precision:
          exact_search || limited ? minimum : nil,
        cell_budget_limited: exact_search ? limited : (limited ? true : nil),
        search_budget: @search_budget.infinite? ? nil : @search_budget,
        nodes_searched: @nodes_searched,
        search_budget_exhausted: search_exhausted,
        optimal: !search_exhausted && exact_search
      )
    end

    def overcoverage_cost(node)
      return 0.0 if node.kind == :inside

      [node.cell_area - node.intersection_area, 0.0].max
    end

    def search_remaining
      @search_budget - @nodes_searched
    end
  end
end
