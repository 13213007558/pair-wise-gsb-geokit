require 'geokit/spherical/vector'
require 'geokit/spherical/errors'
require 'geokit/spherical/core'
require 'geokit/spherical/polygon'
require 'geokit/spherical/boolean'

module Geokit
  # Spherical polygon geometry: polygons whose edges are great-circle
  # short arcs (strictly shorter than half a great circle) on the unit
  # sphere, with boolean operations (union / intersection / difference)
  # that work unchanged across the antimeridian and around the poles.
  #
  # Ring orientation convention: the polygon interior is on the LEFT of
  # travel. Outer rings are counterclockwise, holes are clockwise. Input
  # rings are re-oriented automatically when needed; each input ring must
  # enclose less than a hemisphere so its interior is unambiguous.
  #
  # Everything is computed locally in 3D unit-vector space; no planar
  # projection and no network access is involved.
  module Spherical
    def self.union(a, b)
      Boolean.union(a, b)
    end

    def self.intersection(a, b)
      Boolean.intersection(a, b)
    end

    def self.difference(a, b)
      Boolean.difference(a, b)
    end
  end
end
