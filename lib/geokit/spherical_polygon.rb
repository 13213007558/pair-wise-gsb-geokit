module Geokit
  # Spherical polygon represented as an outer ring and zero or more holes.
  # Longitudes may cross the antimeridian; rings are unwrapped during use.
  class SphericalPolygon
    RADIAN_SCALE = Math::PI / 180.0

    attr_reader :rings, :radius

    def initialize(rings, radius = 6_376_772.71)
      @rings = normalize_rings(rings)
      @radius = radius.to_f
      raise ArgumentError, 'polygon must contain an outer ring' if @rings.empty?
      raise ArgumentError, 'polygon rings must each contain three points' unless @rings.all? { |ring| ring.length >= 3 }
      raise ArgumentError, 'polygon has no spherical area' if area <= 1e-7
    end

    def area
      outer_area = ring_area(@rings[0])
      hole_area = @rings[1..-1].sum { |ring| ring_area(ring) }
      [outer_area - hole_area, 0.0].max
    end

    def intersection_area(rectangle)
      clipped_outer = clipped_rings(@rings[0], rectangle)
      outer = clipped_outer.sum { |ring| ring_area(ring) }
      holes = @rings[1..-1].sum do |ring|
        clipped_rings(ring, rectangle).sum { |candidate| ring_area(candidate) }
      end
      [outer - holes, 0.0].max
    end

    private

    def normalize_rings(rings)
      source = rings.is_a?(Hash) ? rings_from_hash(rings) : Array(rings)
      multi_ring = source.any? { |item| ring_like?(item) }
      source = [source] unless multi_ring
      source.map { |ring| normalize_ring(Array(ring)) }.reject { |ring| ring.length < 3 }
    end

    def rings_from_hash(hash)
      coordinates = hash[:coordinates] || hash['coordinates']
      type = hash[:type] || hash['type']
      if coordinates && type
        unless type.to_s == 'Polygon'
          raise ArgumentError, "unsupported geometry type #{type.inspect}"
        end
        coordinates.map { |ring| ring.map { |point| [point[1], point[0]] } }
      elsif (outer = hash[:outer_ring] || hash['outer_ring'])
        holes = hash[:holes] || hash['holes'] || []
        [outer, *Array(holes)]
      else
        raise ArgumentError, 'polygon hash must contain coordinates or outer_ring'
      end
    end

    def point?(value)
      value.is_a?(Array) && value.length >= 2 &&
        value[0].is_a?(Numeric) && value[1].is_a?(Numeric)
    end

    def ring_like?(item)
      (item.is_a?(Array) && item.any? && point?(item[0])) ||
        item.respond_to?(:points)
    end

    def normalize_ring(points)
      points = points.points if points.respond_to?(:points)
      normalized = points.map do |point|
        if point.is_a?(Array)
          [point[0].to_f, point[1].to_f]
        else
          location = point.respond_to?(:lat) ? point : Geokit::LatLng.normalize(point)
          [location.lat.to_f, location.lng.to_f]
        end
      end
      normalized.each do |lat, lng|
        raise ArgumentError, 'latitude must be between -90 and 90' unless lat.between?(-90.0, 90.0)
        raise ArgumentError, 'longitude must be finite' unless lng.finite?
      end
      normalized.pop if normalized.length > 1 && normalized.first == normalized.last
      normalize_longitudes(normalized)
    end

    def normalize_longitudes(points)
      points.each do |point|
        point[1] -= 360.0 while point[1] > 180.0
        point[1] += 360.0 while point[1] < -180.0
      end
      unwrap_ring(points)
    end

    def unwrap_ring(points)
      previous = points[0]
      points.each_with_index.map do |point, index|
        delta = point[1] - previous[1]
        shifted = if index > 0 && delta < -180.0 && delta > -360.0
                    [point[0], point[1] + 360.0]
                  elsif index > 0 && delta > 180.0 && delta < 360.0
                    [point[0], point[1] - 360.0]
                  elsif index > 0 && delta.abs == 360.0 && point[0] != previous[0]
                    [point[0], point[1] + (delta.negative? ? 360.0 : -360.0)]
                  else
                    point.dup
                  end
        previous = shifted
        shifted
      end
    end

    def ring_area(ring)
      total = closed_ring(ring).each_cons(2).sum do |start_point, end_point|
        (end_point[1] - start_point[1]) * RADIAN_SCALE *
          (2.0 + Math.sin(start_point[0] * RADIAN_SCALE) +
           Math.sin(end_point[0] * RADIAN_SCALE)) / 2.0
      end
      (total.abs * @radius**2)
    end

    def clipped_rings(ring, rectangle)
      candidate = ring_for_clip(ring, rectangle)
      clipped = clip_ring(candidate, rectangle)
      clipped.length < 3 ? [] : [clipped]
    end

    def ring_for_clip(ring, rectangle)
      west = rectangle[:west]
      global_edge = ring.each_cons(2).find do |start_point, end_point|
        (end_point[1] - start_point[1]).abs == 360.0 &&
          start_point[0] == end_point[0]
      end
      if global_edge
        delta = global_edge[1][1] - global_edge[0][1]
        offset = delta.positive? ? west + 180.0 : west - 180.0
        return ring.map { |lat, lng| [lat, lng + offset] }
      end

      east = rectangle[:east]
      longitudes = ring.map { |_lat, lng| lng }
      if longitudes.any? { |value| value >= west && value <= east }
        shift = 0.0
      elsif longitudes.any? { |value| value + 360.0 >= west && value + 360.0 <= east }
        shift = 360.0
      elsif longitudes.any? { |value| value - 360.0 >= west && value - 360.0 <= east }
        shift = -360.0
      else
        shift = 0.0
      end
      ring.map { |lat, lng| [lat, lng + shift] }
    end

    def clip_ring(ring, rectangle)
      clipped = ring
      clipped = clip_axis(clipped, rectangle[:west], true) { |point, boundary| point[1] >= boundary }
      clipped = clip_axis(clipped, rectangle[:east], true) { |point, boundary| point[1] <= boundary }
      clipped = clip_axis(clipped, rectangle[:south], false) { |point, boundary| point[0] >= boundary }
      clip_axis(clipped, rectangle[:north], false) { |point, boundary| point[0] <= boundary }
    end

    def clip_axis(points, boundary, longitude)
      return [] if points.empty?

      output = []
      closed_ring(points).each_cons(2) do |start_point, end_point|
        start_inside = yield(start_point, boundary)
        end_inside = yield(end_point, boundary)
        if end_inside
          output << intersection(start_point, end_point, boundary, longitude) unless start_inside
          output << end_point
        elsif start_inside
          output << intersection(start_point, end_point, boundary, longitude)
        end
      end
      output
    end

    def closed_ring(ring)
      ring.length > 1 && ring.first == ring.last ? ring : ring + [ring.first]
    end

    def intersection(start_point, end_point, boundary, longitude)
      axis = longitude ? 1 : 0
      denominator = end_point[axis] - start_point[axis]
      ratio = denominator.zero? ? 0.0 : (boundary - start_point[axis]) / denominator
      if longitude
        [start_point[0] + (end_point[0] - start_point[0]) * ratio, boundary]
      else
        [boundary, start_point[1] + (end_point[1] - start_point[1]) * ratio]
      end
    end

  end
end
