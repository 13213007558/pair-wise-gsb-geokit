module Geokit
  # Spherical Delaunay triangulation of stations on the unit sphere.
  #
  # Given a set of stations (Geokit::LatLng, [lat, lng] pairs in degrees, or
  # [x, y, z] cartesian vectors), computes the Delaunay triangulation of the
  # sphere: the convex hull of the stations lifted to 3D, whose faces are the
  # spherical Delaunay triangles. Returns directed triangular faces, edge
  # adjacency and the mapping back to the original station identities.
  #
  # The triangulation is only defined when the stations cover the whole
  # sphere, i.e. they are not contained in any closed hemisphere (equivalently
  # the origin lies strictly inside their convex hull). When this precondition
  # does not hold, a NotApplicableError is raised.
  #
  # Stability rules:
  # * Stations are canonicalized (sorted by coordinates) before triangulation,
  #   so permuting the input never changes the resulting topology.
  # * Duplicate stations (within tolerance) are merged into one canonical
  #   station; index_map keeps the identity of every input point.
  # * Cocircular stations (4 or more sharing a circumcircle) are triangulated
  #   by a deterministic fan around the canonically smallest station.
  # * All computations happen on 3D unit vectors, so nothing special happens
  #   at the 180 degree meridian or at the poles.
  #
  # Numerical degeneracies that cannot be decided reliably (near-cocircular
  # points, near-collinear triples, marginal hemisphere coverage, ...) are
  # collected in Triangulation#degeneracies instead of being silently
  # resolved.
  module SphericalDelaunay
    # Tolerance used for all geometric predicates, in units of chord /
    # projected distance on the unit sphere.
    DEFAULT_TOLERANCE = 1e-9
    # Distances within tolerance are "equal"; distances within
    # NEAR_FACTOR * tolerance are "unreliable" and reported as degenerate.
    NEAR_FACTOR = 100

    class Error < StandardError; end

    # Raised when the stations do not cover the whole sphere (they fit in a
    # closed hemisphere, or there are too few of them), so a spherical
    # Delaunay triangulation is not applicable.
    class NotApplicableError < Error; end

    # A reported numerical degeneracy. kind is a Symbol, message a human
    # readable description, data holds the involved station indices.
    Degeneracy = Struct.new(:kind, :message, :data)

    # An oriented triangular face. vertices are station indices ordered
    # counter-clockwise seen from outside the sphere. normal is the outward
    # unit normal of the supporting plane, which is also the direction of the
    # spherical circumcenter; offset is the distance from the origin to the
    # plane (cosine of the circumradius).
    Face = Struct.new(:vertices, :normal, :offset) do
      def directed_edges
        v = vertices
        [[v[0], v[1]], [v[1], v[2]], [v[2], v[0]]]
      end

      def circumcenter
        normal
      end

      def circumradius
        Math.acos([[-1.0, offset].max, 1.0].min)
      end
    end

    # An undirected edge with its adjacency: the two incident faces and the
    # two oppositely directed half-edges [from, to, face_index].
    Edge = Struct.new(:vertices, :faces, :half_edges)

    # The result of a triangulation.
    class Triangulation
      # Canonical station positions as [x, y, z] unit vectors.
      attr_reader :stations
      # station index => array of input positions merged into that station.
      attr_reader :station_sources
      # input position => station index.
      attr_reader :index_map
      attr_reader :faces, :edges, :degeneracies

      def initialize(stations, station_sources, index_map, faces, edges, degeneracies)
        @stations = stations
        @station_sources = station_sources
        @index_map = index_map
        @faces = faces
        @edges = edges
        @degeneracies = degeneracies
      end

      def vertex_count
        stations.size
      end

      def edge_count
        edges.size
      end

      def face_count
        faces.size
      end

      def euler_characteristic
        vertex_count - edge_count + face_count
      end

      # Total spherical area covered by the faces. Equals 4 * PI for a valid
      # full-sphere triangulation.
      def total_area
        faces.inject(0.0) do |sum, face|
          a, b, c = face.vertices.map { |i| stations[i] }
          sum + SphericalDelaunay.solid_angle(a, b, c)
        end
      end

      def degenerate?
        !degeneracies.empty?
      end
    end

    class << self
      # Triangulates the given stations on the unit sphere.
      #
      # points: array of Geokit::LatLng, [lat, lng] (degrees) or [x, y, z].
      # tolerance: geometric tolerance, see DEFAULT_TOLERANCE.
      #
      # Returns a Triangulation. Raises NotApplicableError when the stations
      # do not cover the whole sphere.
      def triangulate(points, tolerance: DEFAULT_TOLERANCE)
        raise ArgumentError, 'tolerance must be positive' unless tolerance > 0

        vectors = points.map { |p| to_unit_vector(p) }
        stations, index_map, station_sources, degeneracies =
          canonical_stations(vectors, tolerance)

        if stations.size < 4
          raise NotApplicableError,
            "spherical Delaunay triangulation needs at least 4 distinct stations covering the sphere, got #{stations.size}"
        end

        groups = supporting_planes(stations, tolerance, degeneracies)
        if groups.empty?
          raise NotApplicableError,
            'stations do not cover the whole sphere (no closed triangulation found)'
        end

        groups.each_value do |group|
          check_coverage(group, tolerance, degeneracies)
        end

        faces = build_faces(groups, stations, degeneracies)
        edges = build_edges(faces, degeneracies)
        triangulation = Triangulation.new(
          stations, station_sources, index_map, faces, edges, degeneracies
        )
        audit(triangulation, tolerance)
        triangulation
      end

      # Spherical excess (area) of the triangle with unit-vector vertices.
      def solid_angle(a, b, c)
        num = dot(a, cross(b, c)).abs
        den = 1.0 + dot(a, b) + dot(b, c) + dot(c, a)
        2.0 * Math.atan2(num, den)
      end

      def to_unit_vector(point)
        if point.is_a?(Array) && point.size == 3
          v = point.map { |c| Float(c) }
          length = norm(v)
          raise ArgumentError, "zero cartesian vector: #{point.inspect}" if length.zero?

          scale(v, 1.0 / length)
        else
          lat_lng = to_lat_lng(point)
          lat = lat_lng.lat.to_f * Math::PI / 180.0
          lng = lat_lng.lng.to_f * Math::PI / 180.0
          [Math.cos(lat) * Math.cos(lng), Math.cos(lat) * Math.sin(lng), Math.sin(lat)]
        end
      end

      def dot(a, b)
        a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
      end

      def cross(a, b)
        [
          a[1] * b[2] - a[2] * b[1],
          a[2] * b[0] - a[0] * b[2],
          a[0] * b[1] - a[1] * b[0]
        ]
      end

      def sub(a, b)
        [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
      end

      def add(a, b)
        [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
      end

      def scale(a, factor)
        [a[0] * factor, a[1] * factor, a[2] * factor]
      end

      def norm(a)
        Math.sqrt(dot(a, a))
      end

      def normalize(a)
        scale(a, 1.0 / norm(a))
      end

      private

      def to_lat_lng(point)
        return point if point.is_a?(Geokit::LatLng)
        return Geokit::LatLng.new(point[0], point[1]) if point.is_a?(Array) && point.size == 2
        return point.to_lat_lng if point.respond_to?(:to_lat_lng)

        raise ArgumentError, "#{point.inspect} cannot be used as a station (want LatLng, [lat, lng] or [x, y, z])"
      end

      # Merges duplicate stations and sorts the rest into a canonical
      # (coordinate) order so the result does not depend on input order.
      def canonical_stations(vectors, tolerance)
        order = (0...vectors.size).sort_by { |i| vectors[i].map { |c| c.round(10) } + [i] }
        clusters = [] # [representative vector, [input indices]]
        order.each do |i|
          vector = vectors[i]
          cluster = clusters.find { |rep, _| norm(sub(rep, vector)) <= tolerance }
          if cluster
            cluster[1] << i
          else
            clusters << [vector, [i]]
          end
        end

        stations = []
        station_sources = []
        index_map = {}
        degeneracies = []
        clusters.each_with_index do |(rep, members), station_index|
          stations << normalize(rep)
          station_sources << members.sort
          members.each { |input_index| index_map[input_index] = station_index }
          next unless members.size > 1

          degeneracies << Degeneracy.new(
            :duplicate_stations,
            "#{members.size} input stations coincide within tolerance and were merged into station #{station_index}",
            members.sort
          )
        end
        [stations, index_map, station_sources, degeneracies]
      end

      # Finds every supporting plane of the station set: a plane through 3+
      # stations with all remaining stations strictly on the origin side.
      # Returns a hash of sorted-vertex-set => {normal:, offset:, vertices:}.
      def supporting_planes(stations, tolerance, degeneracies)
        groups = {}
        near = NEAR_FACTOR * tolerance
        count = stations.size
        (0...count).to_a.combination(3) do |i, j, k|
          a = stations[i]
          normal = cross(sub(stations[j], a), sub(stations[k], a))
          length = norm(normal)
          if length < tolerance
            degeneracies << Degeneracy.new(
              :collinear_triple,
              "stations #{i}, #{j}, #{k} are nearly collinear on a great circle; triple skipped",
              [i, j, k]
            )
            next
          end
          normal = scale(normal, 1.0 / length)
          offset = dot(normal, a)
          dots = stations.map { |p| dot(normal, p) }
          if dots.max - offset <= tolerance
            # normal already points outward
          elsif offset - dots.min <= tolerance
            normal = scale(normal, -1.0)
            offset = -offset
            dots = dots.map { |d| -d }
          else
            next
          end

          on_plane = (0...count).select { |s| (dots[s] - offset).abs <= tolerance }.sort
          key = on_plane
          next if groups.key?(key)

          groups[key] = { normal: normal, offset: offset, vertices: on_plane }
          marginal = (0...count).select do |s|
            distance = (dots[s] - offset).abs
            distance > tolerance && distance <= near
          end
          marginal.each do |s|
            degeneracies << Degeneracy.new(
              :near_degenerate,
              "station #{s} lies within #{NEAR_FACTOR}x tolerance of the face plane through #{on_plane.inspect}; classification is unreliable",
              on_plane + [s]
            )
          end
        end
        groups
      end

      # A supporting plane tangent to or cutting through the origin means the
      # stations fit in a closed hemisphere: the triangulation is undefined.
      def check_coverage(group, tolerance, degeneracies)
        offset = group[:offset]
        if offset <= tolerance
          raise NotApplicableError,
            "stations are contained in a closed hemisphere (face plane through #{group[:vertices].inspect} passes within tolerance of the sphere center)"
        elsif offset <= NEAR_FACTOR * tolerance
          degeneracies << Degeneracy.new(
            :marginal_coverage,
            "face plane through #{group[:vertices].inspect} passes within #{NEAR_FACTOR}x tolerance of the sphere center; hemisphere coverage is marginal",
            group[:vertices]
          )
        end
      end

      def build_faces(groups, stations, degeneracies)
        faces = []
        groups.each_value do |group|
          vertices = group[:vertices]
          normal = group[:normal]
          if vertices.size > 3
            degeneracies << Degeneracy.new(
              :cocircular_stations,
              "#{vertices.size} stations #{vertices.inspect} share a circumcircle; resolved by canonical fan triangulation",
              vertices
            )
          end
          fan_triangles(vertices, normal, stations).each do |triangle|
            faces << Face.new(triangle, normal, group[:offset])
          end
        end
        faces.sort_by { |face| face.vertices.sort }
      end

      # Splits a (possibly > 3 vertices) cocircular polygon into triangles.
      # The polygon vertices are sorted angularly around their centroid and
      # fanned from the canonically smallest station index, which makes the
      # diagonal choice deterministic and input-order independent.
      def fan_triangles(vertices, normal, stations)
        pivot = vertices.min
        centroid = scale(vertices.inject([0.0, 0.0, 0.0]) { |sum, v| add(sum, stations[v]) }, 1.0 / vertices.size)
        e1 = normalize(sub(stations[pivot], centroid))
        e2 = cross(normal, e1)
        sorted = vertices.sort_by do |v|
          relative = sub(stations[v], centroid)
          Math.atan2(dot(relative, e2), dot(relative, e1))
        end
        sorted = sorted.rotate(sorted.index(pivot))
        (1...(sorted.size - 1)).map do |i|
          orient([sorted[0], sorted[i], sorted[i + 1]], normal, stations)
        end
      end

      # Orders a triangle counter-clockwise seen from outside (normal side).
      def orient(triangle, normal, stations)
        a, b, c = triangle.map { |v| stations[v] }
        if dot(cross(sub(b, a), sub(c, a)), normal) < 0
          [triangle[0], triangle[2], triangle[1]]
        else
          triangle
        end
      end

      def build_edges(faces, degeneracies)
        directed = Hash.new { |hash, key| hash[key] = [] }
        faces.each_with_index do |face, face_index|
          face.directed_edges.each { |edge| directed[edge] << face_index }
        end

        undirected = Hash.new { |hash, key| hash[key] = [] }
        directed.each do |(from, to), face_indices|
          key = from < to ? [from, to] : [to, from]
          undirected[key] << [from, to, face_indices]
        end

        undirected.keys.sort.map do |key|
          halves = undirected[key]
          face_indices = halves.flat_map { |_, _, fs| fs }
          unless halves.size == 2 && halves.all? { |_, _, fs| fs.size == 1 } &&
                 halves[0][1] == halves[1][0] && halves[0][0] == halves[1][1]
            degeneracies << Degeneracy.new(
              :non_manifold_edge,
              "edge #{key.inspect} is incident to #{face_indices.size} face(s) instead of 2 with opposite orientations",
              key
            )
          end
          Edge.new(key, face_indices, halves)
        end
      end

      # Whole-triangulation sanity checks; failures are reported, not raised.
      def audit(triangulation, tolerance)
        if triangulation.euler_characteristic != 2
          triangulation.degeneracies << Degeneracy.new(
            :euler_characteristic,
            "V - E + F = #{triangulation.euler_characteristic}, expected 2 for a sphere triangulation",
            nil
          )
        end
        area_error = (triangulation.total_area - 4.0 * Math::PI).abs
        if area_error > NEAR_FACTOR * tolerance
          triangulation.degeneracies << Degeneracy.new(
            :total_area,
            "faces cover #{triangulation.total_area} steradians instead of 4*pi (error #{area_error})",
            nil
          )
        end
      end
    end
  end
end

