module Geokit
  # A complex polygon made of multiple points.  End point must equal start point to close the poly.
  class Polygon
    attr_accessor :points, :holes

    # Pass in an array of Geokit::LatLng
    def initialize(points, holes = [])
      @points = points
      @holes = holes

      # A Polygon must be 'closed', the last point equal to the first point
      # Append the first point to the array to close the polygon
      @points << points[0] if points[0] != points[-1]
    end

    def contains?(point)
      last_point = @points[-1]
      oddNodes = false
      x = point.lng
      y = point.lat

      @points.each do |p|
        yi = p.lat
        xi = p.lng
        yj = last_point.lat
        xj = last_point.lng
        if yi < y && yj >= y ||
            yj < y && yi >= y
          oddNodes = !oddNodes if xi + (y - yi) / (yj - yi) * (xj - xi) < x
        end

        last_point = p
      end

      oddNodes
    end # contains?

    # A polygon is static and can not be updated with new points, as a result
    # calculate the longitude/latitude planar centroid once when requested.
    # Use spherical_center for the true first-surface-moment center.
    def centroid
      @centroid ||= calculate_centroid
    end # end centroid

    def spherical_area(radius: Geokit::SphericalPolygon::DEFAULT_RADIUS)
      to_spherical_polygon(radius: radius).area
    end

    def spherical_first_moment(radius: Geokit::SphericalPolygon::DEFAULT_RADIUS)
      to_spherical_polygon(radius: radius).first_moment
    end

    def spherical_center(epsilon: Geokit::SphericalPolygon::CENTER_EPSILON)
      to_spherical_polygon.center(epsilon: epsilon)
    end
    alias_method :area, :spherical_area
    alias_method :first_moment, :spherical_first_moment
    alias_method :center, :spherical_center

    def to_spherical_polygon(radius: Geokit::SphericalPolygon::DEFAULT_RADIUS)
      Geokit::SphericalPolygon.new(@points, @holes, radius: radius)
    end

    private

    def calculate_centroid
      centroid_lat = 0.0
      centroid_lng = 0.0
      signed_area = 0.0

      # Iterate over each element in the list but the last item as it's
      # calculated by the i+1 logic
      @points[0...-1].each_index do |i|
        x0 = @points[i].lat
        y0 = @points[i].lng
        x1 = @points[i + 1].lat
        y1 = @points[i + 1].lng
        a = (x0 * y1) - (x1 * y0)
        signed_area += a
        centroid_lat += (x0 + x1) * a
        centroid_lng += (y0 + y1) * a
      end

      signed_area *= 0.5
      centroid_lat /= (6.0 * signed_area)
      centroid_lng /= (6.0 * signed_area)

      Geokit::LatLng.new(centroid_lat, centroid_lng)
    end # end calculate_centroid
  end # class Polygon
end
