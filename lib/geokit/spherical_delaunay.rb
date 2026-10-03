# frozen_string_literal: true

module Geokit
  # Spherical Delaunay triangulation of points on a unit sphere.
  #
  # The spherical Delaunay complex of a point set on the unit sphere equals
  # the boundary complex of the Euclidean convex hull of the points lifted to
  # R3: every triangular hull facet has an empty spherical circum-cap, the
  # spherical analogue of the planar empty circumcircle test.
  #
  # The implementation builds the hull incrementally, then rebuilds every
  # supporting plane with exact rational predicates. Exact duplicates are
  # merged deterministically, exact co-circular configurations are fanned with
  # an order-independent rule, and anything that cannot be decided reliably
  # inside the requested tolerance is reported as a degeneracy.
  #
  # The input must cover the whole sphere and therefore must not be contained
  # in any closed hemisphere.
  module SphericalDelaunay
    # Base error for this module.
    class Error < StandardError; end

    # Raised when no full-sphere triangulation exists: fewer than four
    # affinely independent points, or all points lying in a closed
    # hemisphere.
    class InapplicableError < Error; end

    # Raised in strict mode when data is too close to a degenerate
    # configuration (near duplicates or near co-circular points) for the
    # topology to be guaranteed inside the requested tolerance.
    class DegenerateDataError < Error
      attr_reader :degeneracies

      def initialize(message, degeneracies)
        super(message)
        @degeneracies = degeneracies
      end
    end

    # A single condition worth reporting.
    Degeneracy = Struct.new(:kind, :message, :vertex_indices)

    # The triangulation of a point set.
    class Result
      # Canonical unit vectors used by the computation.
      attr_reader :vertices
      # Faces as [i, j, k], counter-clockwise seen from outside the sphere.
      attr_reader :faces
      # Per-face neighbours across the edge opposite each face vertex.
      attr_reader :adjacency
      # Sorted undirected edges as [i, j].
      attr_reader :edges
      # Canonical vertex index => original input indices mapped to it.
      attr_reader :vertex_original_indices
      # Original input index => canonical vertex index.
      attr_reader :original_to_vertex
      # Sum of spherical face areas (4*PI for a full sphere).
      attr_reader :total_area
      # V - E + F; always 2 for a full-sphere triangulation.
      attr_reader :euler_characteristic
      # Degeneracy records describing how the input was handled.
      attr_reader :degeneracies
      # True when a condition could not be decided inside the tolerance.
      attr_reader :unreliable

      def initialize(vertices, faces, adjacency, edges, mappings, area,
                     degeneracies, unreliable)
        @vertices = vertices
        @faces = faces
        @adjacency = adjacency
        @edges = edges
        @vertex_original_indices = mappings[:vertex_original_indices]
        @original_to_vertex = mappings[:original_to_vertex]
        @total_area = area
        @euler_characteristic = vertices.length - edges.length + faces.length
        @degeneracies = degeneracies
        @unreliable = unreliable
      end

      def vertex_count
        @vertices.length
      end

      def face_count
        @faces.length
      end

      def edge_count
        @edges.length
      end

      # The three directed edges of a face, consistent with its orientation.
      def directed_edges(face_index)
        face = @faces[face_index]
        [[face[0], face[1]], [face[1], face[2]], [face[2], face[0]]]
      end

      # Face index adjacent across the edge opposite the given face vertex.
      def opposite_face(face_index, edge_position)
        @adjacency[face_index][edge_position]
      end

      # Original station indices making up a face.
      def face_original_indices(face_index)
        @faces[face_index].map do |vertex|
          @vertex_original_indices[vertex]
        end
      end
    end

    # Default angular tolerance in radians used for predicate decisions.
    DEFAULT_TOLERANCE = 1e-9

    UNRELIABLE_KINDS = %i[
      near_duplicate
      near_cocircular
      predicate_ambiguous
      unclassified_vertex
      closed_hemisphere_boundary
      area
    ].freeze

    class << self
      # Triangulate points on the unit sphere.
      #
      # +points+ may be Geokit::LatLng instances, [lat, lng] pairs in
      # degrees or already-normalized [x, y, z] vectors.
      #
      # Options:
      # +tolerance+:: angular tolerance in radians (default 1e-9).
      # +strict+:: when true, raises DegenerateDataError for near
      #             degeneracies; otherwise the result carries the records
      #             and marks itself unreliable.
      def triangulate(points, tolerance: DEFAULT_TOLERANCE, strict: false)
        Builder.new(points, tolerance.to_f, strict).build
      end

      # Convert any accepted point representation to a unit [x, y, z] vector.
      def to_unit_vector(point)
        if point.respond_to?(:lat) && point.respond_to?(:lng)
          unit_vector_from_lat_lng(point.lat, point.lng)
        elsif point.is_a?(Array) && point.size == 2
          unit_vector_from_lat_lng(point[0], point[1])
        elsif point.is_a?(Array) && point.size == 3
          unit_vector_from_xyz(point)
        else
          raise ArgumentError,
                "expected LatLng, [lat, lng] or [x, y, z], got #{point.inspect}"
        end
      end

      private

      def unit_vector_from_lat_lng(lat, lng)
        lat = Float(lat)
        lng = Float(lng)
        raise ArgumentError, "latitude out of range: #{lat}" unless
          lat.between?(-90.0, 90.0)

        # Collapse the +/-180 meridian and the poles so a physical location
        # never splits across those seams.
        lng = lng % 360.0
        lng -= 360.0 if lng > 180.0
        lng = 0.0 if lat.abs == 90.0

        phi = lat * Math::PI / 180.0
        theta = lng * Math::PI / 180.0
        cos_phi = Math.cos(phi)
        [cos_phi * Math.cos(theta),
         cos_phi * Math.sin(theta),
         Math.sin(phi)]
      end

      def unit_vector_from_xyz(values)
        values = values.map { |value| Float(value) }
        norm = Math.sqrt(values.map { |value| value * value }.inject(0.0, :+))
        raise ArgumentError, 'zero-length vector cannot be placed on sphere' if
          norm.zero?
        values.map { |value| value / norm }
      end
    end

    require_relative 'spherical_delaunay/builder'
  end
end
