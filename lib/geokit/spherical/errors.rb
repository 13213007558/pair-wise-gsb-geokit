module Geokit
  module Spherical
    # Base class for all spherical geometry errors.
    class Error < StandardError; end

    # Raised when a polygon ring is structurally invalid (too few vertices,
    # degenerate edges, ...).
    class InvalidPolygonError < Error; end

    # Raised when an edge connects (nearly) antipodal endpoints. Such an
    # edge is exactly half a great circle and its short-arc definition is
    # ambiguous, so it is reported instead of guessed.
    class AntipodalEdgeError < InvalidPolygonError; end

    # Raised when a ring crosses or touches itself.
    class SelfIntersectionError < InvalidPolygonError; end

    # Raised when rings of one polygon relate incorrectly (hole outside the
    # outer ring, nested holes, rings crossing each other, ...).
    class RingNestingError < InvalidPolygonError; end
  end
end
