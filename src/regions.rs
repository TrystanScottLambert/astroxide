use std::f64::consts::PI;

const DEG_TO_RAD: f64 = PI / 180.0;

/// Result of point location query
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PointLocation {
    /// Point is inside the spherical polygon
    Inside,
    /// Point is on the boundary of the spherical polygon
    OnBoundary,
    /// Point is outside the spherical polygon
    Outside,
    /// Point is antipodal to reference point X (location undetermined)
    Antipodal,
}

/// Relative bearing between two points on a sphere
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum RelativeBearing {
    /// Target is east of reference
    East,
    /// Target is west of reference
    West,
    /// Target is neither east nor west (same meridian)
    SameMeridian,
}

/// Represents a spherical polygon boundary
pub struct SphericalPolygon {
    /// Vertex latitudes (degrees)
    vertex_latitudes: Vec<f64>,
    /// Vertex longitudes (degrees)
    vertex_longitudes: Vec<f64>,
    /// Number of vertices
    num_vertices: usize,
    /// Interior point X latitude (degrees)
    interior_point_lat: f64,
    /// Interior point X longitude (degrees)
    interior_point_lon: f64,
    /// Transformed longitudes of vertices (in X-centered coordinate system)
    transformed_vertex_longitudes: Vec<f64>,
}

impl SphericalPolygon {
    /// Create a new spherical polygon boundary
    ///
    /// # Arguments
    /// * `vertex_latitudes` - Vector of vertex latitudes in degrees (positive north, negative south)
    /// * `vertex_longitudes` - Vector of vertex longitudes in degrees (positive east, negative west)
    /// * `interior_point_lat` - Latitude of interior point X
    /// * `interior_point_lon` - Longitude of interior point X
    ///
    /// # Returns
    /// Result containing the SphericalPolygon or an error message
    pub fn new(
        vertex_latitudes: Vec<f64>,
        vertex_longitudes: Vec<f64>,
        interior_point_lat: f64,
        interior_point_lon: f64,
    ) -> Result<Self, String> {
        let num_vertices = vertex_latitudes.len();

        if num_vertices != vertex_longitudes.len() {
            return Err(
                "vertex_latitudes and vertex_longitudes must have the same length".to_string(),
            );
        }

        if num_vertices < 3 {
            return Err("Polygon must have at least 3 vertices".to_string());
        }

        // Compute transformed longitudes for each vertex
        let transformed_vertex_longitudes: Vec<f64> = vertex_latitudes
            .iter()
            .zip(&vertex_longitudes)
            .map(|(&lat, &lon)| {
                transform_longitude(interior_point_lat, interior_point_lon, lat, lon)
            })
            .collect();

        // Check that neighboring vertices are distinct and not antipodal
        for i in 0..num_vertices {
            let prev_index = if i == 0 { num_vertices - 1 } else { i - 1 };

            // Check if vertices are distinct
            if vertex_latitudes[i] == vertex_latitudes[prev_index]
                && vertex_longitudes[i] == vertex_longitudes[prev_index]
            {
                return Err(format!(
                    "Vertices {} and {} are not distinct",
                    i, prev_index
                ));
            }

            // Check if interior point lies on great circle through vertices
            if transformed_vertex_longitudes[i] == transformed_vertex_longitudes[prev_index] {
                return Err(format!(
                    "Interior point lies on great circle through vertices {} and {}",
                    i, prev_index
                ));
            }

            // Check if vertices are antipodal
            if are_vertices_antipodal(
                vertex_latitudes[i],
                vertex_longitudes[i],
                vertex_latitudes[prev_index],
                vertex_longitudes[prev_index],
            ) {
                return Err(format!("Vertices {} and {} are antipodal", i, prev_index));
            }
        }

        Ok(SphericalPolygon {
            vertex_latitudes,
            vertex_longitudes,
            num_vertices,
            interior_point_lat,
            interior_point_lon,
            transformed_vertex_longitudes,
        })
    }

    /// Locate a point P relative to the spherical polygon boundary
    ///
    /// # Arguments
    /// * `point_lat` - Latitude of point P in degrees
    /// * `point_lon` - Longitude of point P in degrees
    ///
    /// # Returns
    /// PointLocation indicating whether P is inside, outside, on boundary, or antipodal to X
    pub fn locate_point(&self, point_lat: f64, point_lon: f64) -> PointLocation {
        // Check if P is antipodal to interior point
        if is_antipodal(
            point_lat,
            point_lon,
            self.interior_point_lat,
            self.interior_point_lon,
        ) {
            return PointLocation::Antipodal;
        }

        // Check if P coincides with interior point
        if point_lat == self.interior_point_lat && point_lon == self.interior_point_lon {
            return PointLocation::Inside;
        }

        // Transform P's longitude in interior-point-centered coordinate system
        let transformed_point_lon = transform_longitude(
            self.interior_point_lat,
            self.interior_point_lon,
            point_lat,
            point_lon,
        );

        // Count crossings of arc from interior point to P with polygon boundary
        let mut crossing_count = 0;

        for i in 0..self.num_vertices {
            let next_index = (i + 1) % self.num_vertices;

            let vertex_a = VertexView {
                lat: self.vertex_latitudes[i],
                lon: self.vertex_longitudes[i],
                transformed_lon: self.transformed_vertex_longitudes[i],
            };

            let vertex_b = VertexView {
                lat: self.vertex_latitudes[next_index],
                lon: self.vertex_longitudes[next_index],
                transformed_lon: self.transformed_vertex_longitudes[next_index],
            };

            // Test for necessary strike condition
            let has_necessary_strike =
                self.check_necessary_strike(transformed_point_lon, &vertex_a, &vertex_b);

            if has_necessary_strike {
                // Check if P lies on vertex A
                if point_lat == vertex_a.lat && point_lon == vertex_a.lon {
                    return PointLocation::OnBoundary;
                }

                // Check if arc crosses this polygon side
                match self.check_arc_crossing(point_lat, point_lon, &vertex_a, &vertex_b) {
                    ArcCrossing::OnBoundary => return PointLocation::OnBoundary,
                    ArcCrossing::Crosses => crossing_count += 1,
                    ArcCrossing::NoCrossing => {}
                }
            }
        }

        // Even number of crossings means P is inside
        if crossing_count % 2 == 0 {
            PointLocation::Inside
        } else {
            PointLocation::Outside
        }
    }

    /// Check if the arc from interior point to test point has necessary strike
    /// to potentially cross the polygon side from vertex A to vertex B
    fn check_necessary_strike(
        &self,
        transformed_point_lon: f64,
        vertex_a: &VertexView,
        vertex_b: &VertexView,
    ) -> bool {
        if transformed_point_lon == vertex_a.transformed_lon {
            return true;
        }

        let bearing_a_to_b =
            compute_relative_bearing(vertex_a.transformed_lon, vertex_b.transformed_lon);
        let bearing_a_to_point =
            compute_relative_bearing(vertex_a.transformed_lon, transformed_point_lon);
        let bearing_point_to_b =
            compute_relative_bearing(transformed_point_lon, vertex_b.transformed_lon);

        bearing_a_to_point == bearing_a_to_b && bearing_point_to_b == bearing_a_to_b
    }

    /// Check if the arc from interior point to test point crosses the polygon side
    fn check_arc_crossing(
        &self,
        point_lat: f64,
        point_lon: f64,
        vertex_a: &VertexView,
        vertex_b: &VertexView,
    ) -> ArcCrossing {
        // Transform to A-centered coordinate system
        let transformed_interior_lon = transform_longitude(
            vertex_a.lat,
            vertex_a.lon,
            self.interior_point_lat,
            self.interior_point_lon,
        );
        let transformed_b_lon_from_a =
            transform_longitude(vertex_a.lat, vertex_a.lon, vertex_b.lat, vertex_b.lon);
        let transformed_point_lon_from_a =
            transform_longitude(vertex_a.lat, vertex_a.lon, point_lat, point_lon);

        // Check if P lies on side AB
        if transformed_point_lon_from_a == transformed_b_lon_from_a {
            return ArcCrossing::OnBoundary;
        }

        // Check if interior point and P lie on opposite sides of arc AB
        let bearing_b_to_interior =
            compute_relative_bearing(transformed_b_lon_from_a, transformed_interior_lon);
        let bearing_b_to_point =
            compute_relative_bearing(transformed_b_lon_from_a, transformed_point_lon_from_a);

        if are_bearings_opposite(bearing_b_to_interior, bearing_b_to_point) {
            ArcCrossing::Crosses
        } else {
            ArcCrossing::NoCrossing
        }
    }
}

/// Helper struct to pass vertex information around
struct VertexView {
    lat: f64,
    lon: f64,
    transformed_lon: f64,
}

/// Result of checking if an arc crosses a polygon side
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ArcCrossing {
    /// Arc crosses the polygon side
    Crosses,
    /// Point lies on the polygon boundary
    OnBoundary,
    /// Arc does not cross the polygon side
    NoCrossing,
}

/// Transform longitude of point Q to coordinate system where P is north pole
///
/// Returns the 'longitude' of point Q in a coordinate system where point P
/// acts as the north pole.
fn transform_longitude(pole_lat: f64, pole_lon: f64, target_lat: f64, target_lon: f64) -> f64 {
    if pole_lat == 90.0 {
        return target_lon;
    }

    let pole_lat_rad = pole_lat * DEG_TO_RAD;
    let target_lat_rad = target_lat * DEG_TO_RAD;
    let longitude_diff_rad = (target_lon - pole_lon) * DEG_TO_RAD;

    let numerator = longitude_diff_rad.sin() * target_lat_rad.cos();
    let denominator = target_lat_rad.sin() * pole_lat_rad.cos()
        - target_lat_rad.cos() * pole_lat_rad.sin() * longitude_diff_rad.cos();

    numerator.atan2(denominator) / DEG_TO_RAD
}

/// Determine the relative bearing from reference point to target point
///
/// Returns East, West, or SameMeridian based on the angular difference
fn compute_relative_bearing(reference_lon: f64, target_lon: f64) -> RelativeBearing {
    let longitude_diff = normalize_longitude_difference(target_lon - reference_lon);

    if longitude_diff > 0.0 {
        RelativeBearing::East
    } else if longitude_diff < 0.0 {
        RelativeBearing::West
    } else {
        RelativeBearing::SameMeridian
    }
}

/// Check if two relative bearings are opposite (one East, one West)
fn are_bearings_opposite(bearing1: RelativeBearing, bearing2: RelativeBearing) -> bool {
    matches!(
        (bearing1, bearing2),
        (RelativeBearing::East, RelativeBearing::West)
            | (RelativeBearing::West, RelativeBearing::East)
    )
}

/// Normalize longitude difference to range [-180, 180]
fn normalize_longitude_difference(mut diff: f64) -> f64 {
    while diff > 180.0 {
        diff -= 360.0;
    }
    while diff < -180.0 {
        diff += 360.0;
    }
    diff
}

/// Check if two points are antipodal (180 degrees apart)
fn is_antipodal(lat1: f64, lon1: f64, lat2: f64, lon2: f64) -> bool {
    if lat1 != -lat2 {
        return false;
    }

    let longitude_diff = normalize_longitude_difference(lon1 - lon2);
    longitude_diff.abs() == 180.0
}

/// Check if two vertices are antipodal
fn are_vertices_antipodal(lat1: f64, lon1: f64, lat2: f64, lon2: f64) -> bool {
    is_antipodal(lat1, lon1, lat2, lon2)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simple_polygon() {
        // Define a simple quadrilateral
        let vertex_lats = vec![0.0, 10.0, 10.0, 0.0];
        let vertex_lons = vec![0.0, 0.0, 10.0, 10.0];
        let interior_lat = 5.0;
        let interior_lon = 5.0;

        let poly =
            SphericalPolygon::new(vertex_lats, vertex_lons, interior_lat, interior_lon).unwrap();

        // Test point inside
        assert_eq!(poly.locate_point(5.0, 5.0), PointLocation::Inside);

        // Test point outside
        assert_eq!(poly.locate_point(20.0, 20.0), PointLocation::Outside);

        // Test point on vertex
        assert_eq!(poly.locate_point(0.0, 0.0), PointLocation::OnBoundary);
    }

    #[test]
    fn test_antipodal_vertices() {
        let vertex_lats = vec![45.0, -45.0, 0.0];
        let vertex_lons = vec![0.0, 180.0, 90.0];
        let interior_lat = 0.0;
        let interior_lon = 45.0;

        let result = SphericalPolygon::new(vertex_lats, vertex_lons, interior_lat, interior_lon);
        assert!(result.is_err());
    }
}
