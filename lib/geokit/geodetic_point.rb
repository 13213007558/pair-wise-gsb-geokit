module Geokit
  # A geodetic (curvilinear) position: latitude and longitude in decimal
  # degrees plus ellipsoidal height in metres, optionally tagged with the
  # Datum whose reference ellipsoid the coordinates refer to.
  #
  # A GeodeticPoint is never implicitly associated with a datum: pass one
  # explicitly to the constructor or to Datum#to_ecef.
  class GeodeticPoint
    attr_reader :lat, :lng, :height, :datum

    def initialize(lat, lng, height = 0.0, datum: nil)
      @lat = Geokit.assert_numeric(lat, 'lat')
      @lng = Geokit.assert_numeric(lng, 'lng')
      @height = Geokit.assert_numeric(height, 'height')
      unless @lat.between?(-90.0, 90.0)
        raise ArgumentError, "lat must be between -90 and 90 degrees, " \
                            "got #{@lat}"
      end
      @datum = datum
    end

    # Builds a point from a datum-agnostic LatLng without guessing a datum:
    # the Datum must be supplied by the caller.
    def self.from_lat_lng(lat_lng, datum, height: 0.0)
      new(lat_lng.lat, lat_lng.lng, height, datum: datum)
    end

    # Returns a datum-agnostic LatLng. No conversion is performed and no
    # datum is attached to the returned object.
    def to_lat_lng
      Geokit::LatLng.new(lat, lng)
    end

    # Converts to Earth-Centred, Earth-Fixed Cartesian coordinates. Raises if
    # the point has no explicit datum.
    def to_ecef(datum = @datum)
      raise ArgumentError, 'a Datum is required to convert to ECEF' if datum.nil?
      datum.to_ecef(self)
    end

    def to_a
      [lat, lng, height]
    end

    def ==(other)
      other.is_a?(GeodeticPoint) && lat == other.lat && lng == other.lng &&
        height == other.height && datum == other.datum
    end
    alias eql? ==

    def hash
      [lat, lng, height, datum].hash
    end
  end

  # An Earth-Centred, Earth-Fixed Cartesian position given in metres. The X
  # axis points to the intersection of the reference meridian and equator and
  # the Z axis to the north pole; XYZ form a right-handed system.
  class EcefPoint
    attr_reader :x, :y, :z

    def initialize(x, y, z)
      @x = Geokit.assert_numeric(x, 'x')
      @y = Geokit.assert_numeric(y, 'y')
      @z = Geokit.assert_numeric(z, 'z')
    end

    # ECEF is already Cartesian; returns the point itself so it can be
    # passed transparently into chains such as
    # +point.to_ecef.transform(helmert)+.
    def to_ecef(*)
      self
    end

    # Converts back to geodetic coordinates on the given Datum.
    def to_geodetic(datum)
      datum.from_ecef(self)
    end

    # Applies a Helmert transformation. +direction+ is +:forward+ or
    # +:backward+; +epoch+ (decimal years) is mandatory when the Helmert
    # transformation carries parameter rates.
    def transform(helmert, direction: :forward, epoch: nil)
      case direction
      when :forward  then helmert.forward(self, epoch: epoch)
      when :backward then helmert.backward(self, epoch: epoch)
      else raise ArgumentError,
               "direction must be :forward or :backward, got #{direction.inspect}"
      end
    end

    def to_a
      [x, y, z]
    end

    def ==(other)
      other.is_a?(EcefPoint) && x == other.x && y == other.y && z == other.z
    end
    alias eql? ==

    def hash
      [x, y, z].hash
    end
  end
end
