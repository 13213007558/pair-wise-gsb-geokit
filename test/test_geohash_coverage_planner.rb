require 'test/unit'
require 'geokit'

class GeohashCoveragePlannerTest < Test::Unit::TestCase
  SQUARE = [[0, 0], [0, 20], [20, 20], [20, 0]].freeze
  MERIDIAN = [[10, 170], [-10, 170], [-10, -170], [10, -170]].freeze
  POLAR_CAP = [[80, -180], [80, 180], [90, 180], [90, -180]].freeze
  DONUT = [
    [[0, 0], [0, 30], [30, 30], [30, 0]],
    [[10, 10], [20, 10], [20, 20], [10, 20]]
  ].freeze

  def test_returns_adaptive_antichain_and_reports_real_overcoverage
    result = Geokit::GeohashCoveragePlanner.plan(
      SQUARE, max_precision: 2, max_cells: 100, search_budget: nil
    )

    assert_equal %w[s0 s1 s2 s3 s4 s5 s6 s7], result.prefixes
    assert result.antichain?
    assert result.optimal?
    assert_equal 8, result.cell_count
    assert_in_delta 1_256_163_639_289.86, result.overcoverage_area, 1
    assert_equal result.polygon_area + result.overcoverage_area,
                 covered_area(result)
  end

  def test_conservatively_enlarges_when_cell_budget_is_tight
    result = Geokit::GeohashCoveragePlanner.plan(
      SQUARE, max_precision: 2, max_cells: 1, search_budget: nil
    )

    assert_equal ['s'], result.prefixes
    assert result.optimal?
    assert result.cell_budget_limited?
    assert_equal 8, result.minimum_cell_count_at_max_precision
    assert result.covered_area >= result.polygon_area
  end

  def test_exhaustive_enumeration_lists_all_legal_coverings
    covers = Geokit::GeohashCoveragePlanner.each_covering(
      SQUARE, max_precision: 2
    ).to_a
    expected = [
      [''],
      ['s'],
      %w[s0 s1 s2 s3 s4 s5 s6 s7]
    ]

    assert_equal expected, covers
    covers.each { |prefixes| assert_antichain(prefixes) }
  end

  def test_exhaustive_planner_agrees_with_dynamic_programming
    planned = Geokit::GeohashCoveragePlanner.plan(
      SQUARE, max_precision: 2, max_cells: 8, search_budget: nil
    )
    exhaustive = Geokit::GeohashCoveragePlanner.exhaustive_plan(
      SQUARE, max_precision: 2, max_cells: 8
    )

    assert_equal planned.prefixes, exhaustive.prefixes
    assert_in_delta planned.overcoverage_area,
                    exhaustive.overcoverage_area, 1e-6
    assert exhaustive.optimal?
  end

  def test_holes_reduce_polygon_area_but_stay_covered
    result = Geokit::GeohashCoveragePlanner.plan(
      DONUT, max_precision: 2, max_cells: 100, search_budget: nil
    )
    solid_area = Geokit::SphericalPolygon.new([DONUT[0]]).area

    assert result.optimal?
    assert result.polygon_area < solid_area
    assert_equal result.polygon_area + result.overcoverage_area,
                 covered_area(result)
    assert_antichain(result.prefixes)
  end

  def test_covers_polygon_crossing_the_antimeridian
    result = Geokit::GeohashCoveragePlanner.plan(
      MERIDIAN, max_precision: 1, max_cells: 100, search_budget: nil
    )

    expected = %w[2 8 r x]
    assert_equal expected, result.prefixes
    assert result.optimal?
    assert_in_delta 85_401_410_986_588.78,
                    result.overcoverage_area, 1
  end

  def test_covers_polar_cap
    result = Geokit::GeohashCoveragePlanner.plan(
      POLAR_CAP, max_precision: 1, max_cells: 100, search_budget: nil
    )

    assert_equal %w[b c f g u v y z], result.prefixes
    assert result.optimal?
    assert_in_delta 3_881_537_226_821.76,
                    result.polygon_area, 1
  end

  def test_search_truncation_is_reported_separately
    result = Geokit::GeohashCoveragePlanner.plan(
      DONUT, max_precision: 5, max_cells: 100, search_budget: 5
    )

    assert result.search_budget_exhausted?
    refute result.optimal?
    assert_nil result.cell_budget_limited
    assert_nil result.minimum_cell_count_at_max_precision
    assert_equal 5, result.search_budget
    assert result.nodes_searched <= 5
  end

  def test_accepts_geojson_polygon_hashes
    geojson = {
      'type' => 'Polygon',
      'coordinates' => [
        [[0, 0], [20, 0], [20, 20], [0, 20], [0, 0]]
      ]
    }
    result = Geokit::GeohashCoveragePlanner.plan(
      geojson, max_precision: 2, max_cells: 100, search_budget: nil
    )

    assert result.optimal?
    assert_equal 8, result.cell_count
  end

  def test_validates_search_budget_and_cell_budget
    assert_raise(ArgumentError) do
      Geokit::GeohashCoveragePlanner.new(
        max_precision: 2, max_cells: 0
      )
    end
    assert_raise(ArgumentError) do
      Geokit::GeohashCoveragePlanner.new(
        max_precision: 2, max_cells: 10, search_budget: 0
      )
    end
  end

  private

  def covered_area(result)
    result.prefixes.sum do |prefix|
      bounds = Geokit::Geohash.bounds(prefix)
      rectangle_area(bounds)
    end
  end

  def rectangle_area(bounds)
    radius = Geokit::GeohashCoveragePlanner::EARTH_RADIUS
    radius**2 * ((bounds[:east] - bounds[:west]) * Math::PI / 180.0) *
      (Math.sin(bounds[:north] * Math::PI / 180.0) -
       Math.sin(bounds[:south] * Math::PI / 180.0))
  end

  def assert_antichain(prefixes)
    prefixes.combination(2).each do |left, right|
      refute left.start_with?(right), "#{left} contains #{right}"
      refute right.start_with?(left), "#{right} contains #{left}"
    end
  end
end
