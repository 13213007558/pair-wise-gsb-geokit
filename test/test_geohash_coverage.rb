# frozen_string_literal: true

require File.join(File.dirname(__FILE__), 'helper')

class GeohashCoverageTest < Test::Unit::TestCase
  EARTH_RADIUS = Geokit::Mappable::ClassMethods::EARTH_RADIUS_IN_METERS

  def rectangle(south:, west:, north:, east:)
    [
      [south, west],
      [south, east],
      [north, east],
      [north, west],
      [south, west]
    ].map { |lat, lng| Geokit::LatLng.new(lat, lng) }
  end

  def polar_triangle
    [
      [80, 0],
      [80, 120],
      [80, -120],
      [80, 0]
    ].map { |lat, lng| Geokit::LatLng.new(lat, lng) }
  end

  def antichain?(prefixes)
    prefixes.none? do |prefix|
      prefixes.any? do |other|
        prefix != other && other.start_with?(prefix)
      end
    end
  end

  def prefix_area(prefix)
    Geokit::GeohashCoverage::Node.new(prefix).area
  end

  def test_returns_spherical_areas_for_a_square
    polygon = rectangle(south: 40, west: -10, north: 50, east: 10)
    result = Geokit::GeohashCoverage.plan(
      polygon, max_precision: 4, max_cells: 1000
    )
    flat_area = EARTH_RADIUS**2 *
                (20.0 * Math::PI / 180.0) *
                (10.0 * Math::PI / 180.0) *
                Math.cos(45.0 * Math::PI / 180.0)

    assert_in_delta flat_area, result.polygon_area, flat_area * 0.01
    assert_equal result.cell_area - result.polygon_area,
                 result.overcoverage_area, 1.0
    assert result.optimal?
    assert !result.search_budget_reached?
    assert antichain?(result.prefixes)
    assert result.prefixes.map(&:length).uniq.length > 1
  end

  def test_uses_standard_geohash_character_bit_order
    assert_equal(
      { south: 45.0, north: 90.0, west: 0.0, east: 45.0 },
      Geokit::GeohashCoverage.bounds_for_prefix('u')
    )
    assert_equal(
      { south: 45.0, north: 90.0, west: 45.0, east: 90.0 },
      Geokit::GeohashCoverage.bounds_for_prefix('v')
    )
  end

  def test_small_tree_enumeration_matches_exact_planner
    polygon = rectangle(south: 40, west: -10, north: 60, east: 10)
    covers = Geokit::GeohashCoverage.all_covers(
      polygon, max_precision: 2
    )
    result = Geokit::GeohashCoverage.plan(
      polygon, max_precision: 2, max_cells: 1000
    )
    best_area = covers
                .select { |prefixes| prefixes.length <= 1000 }
                .map { |prefixes| prefixes.sum { |prefix| prefix_area(prefix) } }
                .min

    assert result.optimal?
    assert !result.search_budget_reached?
    assert_in_delta best_area, result.cell_area, 1.0
    best_cover = covers
                 .select { |prefixes| prefixes.length <= 1000 }
                 .min_by { |prefixes|
                   prefixes.sum { |prefix| prefix_area(prefix) }
                 }
                 .sort
    assert_equal best_cover, result.prefixes.sort
  end

  def test_supports_holes
    polygon = rectangle(south: 38, west: -12, north: 52, east: 12)
    hole = rectangle(south: 42, west: -6, north: 48, east: 6)
    result = Geokit::GeohashCoverage.plan(
      polygon,
      holes: [hole],
      max_precision: 3,
      max_cells: 500
    )

    assert result.optimal?
    assert antichain?(result.prefixes)
    assert result.polygon_area.positive?
    assert_in_delta result.cell_area - result.polygon_area,
                    result.overcoverage_area, 1.0
    assert result.overcoverage_area > 100_000_000_000
  end

  def test_supports_antimeridian
    polygon = rectangle(south: 40, west: 170, north: 50, east: -170)
    reference = rectangle(south: 40, west: -10, north: 50, east: 10)
    result = Geokit::GeohashCoverage.plan(
      polygon, max_precision: 3, max_cells: 500
    )
    reference_result = Geokit::GeohashCoverage.plan(
      reference, max_precision: 3, max_cells: 500
    )
    bounds = result.prefixes.map do |prefix|
      Geokit::GeohashCoverage.bounds_for_prefix(prefix)
    end

    assert result.optimal?
    assert antichain?(result.prefixes)
    assert bounds.any? { |bound| bound[:east] == 180.0 }
    assert bounds.any? { |bound| bound[:west] == -180.0 }
    assert_in_delta reference_result.polygon_area,
                    result.polygon_area, 1_000_000
    assert_in_delta reference_result.cell_area,
                    result.cell_area, 1.0
  end

  def test_supports_polar_region
    result = Geokit::GeohashCoverage.plan(
      polar_triangle, max_precision: 3, max_cells: 500
    )

    assert result.optimal?
    assert !result.search_budget_reached?
    assert result.polygon_area > 1_000_000_000_000
    assert result.overcoverage_area.positive?
    assert antichain?(result.prefixes)
    assert result.prefixes.length < 500
  end

  def test_cell_budget_conservatively_expands_to_root
    polygon = rectangle(south: 40, west: -10, north: 50, east: 10)
    result = Geokit::GeohashCoverage.plan(
      polygon, max_precision: 3, max_cells: 1
    )

    assert_equal [''], result.prefixes
    assert_equal 1, result.cell_count
    assert result.optimal?
    assert result.cell_budget_reached?
    assert result.cell_budget_limited?
    assert !result.search_budget_reached?
    assert !result.max_precision_reached?
    assert_in_delta 4 * Math::PI * EARTH_RADIUS**2,
                    result.cell_area, 1.0
  end

  def test_reports_max_precision_boundary
    polygon = rectangle(south: 40, west: -10, north: 50, east: 10)
    result = Geokit::GeohashCoverage.plan(
      polygon, max_precision: 1, max_cells: 32
    )

    assert result.optimal?
    assert result.max_precision_reached?
    assert !result.cell_budget_limited?
    assert result.overcoverage_area.positive?
  end

  def test_search_budget_shortfall_is_separate_from_cell_budget
    polygon = rectangle(south: 40, west: -10, north: 50, east: 10)
    result = Geokit::GeohashCoverage.plan(
      polygon,
      max_precision: 3,
      max_cells: 100,
      search_budget: 1
    )

    assert_equal [''], result.prefixes
    assert result.search_budget_reached?
    assert !result.optimal?
    assert !result.cell_budget_limited?
    assert result.overcoverage_area.finite?
  end

  def test_partial_search_still_returns_a_valid_cover
    polygon = rectangle(south: 40, west: -10, north: 50, east: 10)
    result = Geokit::GeohashCoverage.plan(
      polygon,
      max_precision: 6,
      max_cells: 500,
      search_budget: 1000
    )

    assert result.search_budget_reached?
    assert !result.optimal?
    assert !result.cell_budget_limited?
    assert result.cell_count > 1
    assert result.cell_count <= 500
    assert antichain?(result.prefixes)
    assert result.overcoverage_area.finite?
  end
end
