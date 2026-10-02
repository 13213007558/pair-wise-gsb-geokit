module Geokit
  # Deterministic 3D convex hull support enumeration.
  module SphericalCapHull
    include SphericalCapVectors

    private

    def enumerate_triangle_planes(records)
      frame = affine_frame(records)
      if frame.length == 4
        convex_hull_faces(records).each do |face|
          plane = triangle_normal(records, face)
          next unless plane
          yield plane
          yield scale(plane, -1.0)
        end
      elsif frame.length == 3
        plane = frame.last
        yield plane
        yield scale(plane, -1.0)
      end
    end

    def affine_frame(records)
      first = records[0]
      second = records[1..-1].find do |record|
        norm(subtract(record[:vector], first[:vector])) > 1e-12
      end
      return [first] unless second

      first_edge = subtract(second[:vector], first[:vector])
      third = records.find do |record|
        norm(cross(first_edge, subtract(record[:vector], first[:vector]))) > 1e-12
      end
      return [first, second] unless third

      second_edge = subtract(third[:vector], first[:vector])
      normal = normalize(cross(first_edge, second_edge))
      fourth = records.find do |record|
        dot(normal, subtract(record[:vector], first[:vector])).abs > 1e-12
      end
      return [first, second, normal] unless fourth
      [first, second, third, fourth]
    end

    def convex_hull_faces(records)
      indices = initial_tetrahedron(records)
      return [] if indices.length < 4
      faces = indices.combination(3).map do |face|
        orient_face(records, face, (indices - face).first)
      end

      (records.each_index.to_a - indices).each do |point_index|
        visible = faces.select do |face|
          side(records, face, hull_vector(records[point_index])) > 1e-10
        end
        next if visible.empty?

        boundary = boundary_edges(faces, visible)
        faces -= visible
        boundary.each { |edge| faces << [edge[0], edge[1], point_index] }
      end
      faces
    end

    def initial_tetrahedron(records)
      indices = [0]
      first = hull_vector(records[0])
      second_index = records.each_index.drop(1).find do |index|
        norm(subtract(hull_vector(records[index]), first)) > 1e-9
      end
      indices << second_index unless second_index.nil?

      second_edge = subtract(hull_vector(records[second_index]), first)
      third_index = records.each_index.find do |index|
        norm(cross(second_edge,
                   subtract(hull_vector(records[index]), first))) > 1e-9
      end
      indices << third_index unless third_index.nil?

      third_edge = subtract(hull_vector(records[third_index]), first)
      normal = cross(second_edge, third_edge)
      fourth_index = records.each_index.find do |index|
        dot(normal, subtract(hull_vector(records[index]), first)).abs > 1e-9
      end
      indices << fourth_index unless fourth_index.nil?
      indices
    end

    def orient_face(records, face, opposite_index)
      normal = hull_face_normal(records, face)
      origin = hull_vector(records[face.first])
      opposite = hull_vector(records[opposite_index])
      return face if dot(normal, subtract(opposite, origin)).negative?
      [face[0], face[2], face[1]]
    end

    def boundary_edges(faces, visible)
      invisible_edges = (faces - visible).each_with_object({}) do |face, edges|
        directed_edges(face).each { |edge| edges[edge] = true }
      end
      visible.flat_map { |face| directed_edges(face) }.select do |edge|
        invisible_edges[[edge[1], edge[0]]]
      end
    end

    def directed_edges(face)
      [[face[0], face[1]], [face[1], face[2]], [face[2], face[0]]]
    end

    def triangle_normal(records, face)
      first = records[face[0]][:vector]
      second = records[face[1]][:vector]
      third = records[face[2]][:vector]
      plane = cross(subtract(second, first), subtract(third, first))
      return nil if norm(plane) < 1e-12
      normalize(plane)
    end

    def hull_face_normal(records, face)
      first = hull_vector(records[face[0]])
      second = hull_vector(records[face[1]])
      third = hull_vector(records[face[2]])
      normalize(cross(subtract(second, first), subtract(third, first)))
    end

    def side(records, face, vector)
      origin = hull_vector(records[face.first])
      dot(hull_face_normal(records, face), subtract(vector, origin))
    end

    def hull_vector(record)
      offset = record[:vector].each_with_index.map do |value, index|
        hash = value.hash ^ (index * 7919)
        ((hash.abs % 1009) - 504) / 504.0
      end
      normalize(add(record[:vector], scale(offset, 1e-10)))
    end
  end
end
