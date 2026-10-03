require 'test/unit'
require 'geokit/geohash'

class GeohashTest < Test::Unit::TestCase
  def test_encodes_standard_prefixes
    assert_equal 's0000000d6ds', Geokit::Geohash.encode(0.0001, 0.0001, 12)
    assert_equal 'u4pruydqqvj8', Geokit::Geohash.encode(57.64911, 10.40744, 12)
  end

  def test_reports_exact_cell_bounds
    assert_equal({
      south: -90.0, west: -180.0, north: 90.0, east: 180.0
    }, Geokit::Geohash.bounds(''))
    assert_equal({
      south: 0.0, west: 0.0, north: 45.0, east: 45.0
    }, Geokit::Geohash.bounds('s'))
  end

  def test_rejects_invalid_prefixes_and_precision
    assert_raise(ArgumentError) { Geokit::Geohash.bounds('a') }
    assert_raise(ArgumentError) { Geokit::Geohash.encode(0, 0, 13) }
  end
end
