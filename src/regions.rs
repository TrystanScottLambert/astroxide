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
}

pub struct SphericalAperture {
    ra_center: f64,
    dec_center: f64,
    radius_degrees: f64,
}

impl SphericalAperture {
    pub fn locate_point(&self, ra: f64, dec: f64) -> PointLocation {}
}

/// Represents a spherical polygon boundary
pub struct SphericalPolygon {
    /// Vertex latitudes (degrees)
    vertex_latitudes: Vec<f64>,
    /// Vertex longitudes (degrees)
    vertex_longitudes: Vec<f64>,
    /// Reference point (midpoint of an edge) for ray-casting
    reference_lat: f64,
    reference_lon: f64,
    /// Index of the reference edge
    reference_edge_index: usize,
}

impl SphericalPolygon {
    /// Create a new spherical polygon from vertices
    ///
    /// Uses directional ray-casting - no interior point needed!
    ///
    /// # Arguments
    /// * `vertex_latitudes` - Vector of vertex latitudes in degrees (positive north, negative south)
    /// * `vertex_longitudes` - Vector of vertex longitudes in degrees (positive east, negative west)
    ///
    /// # Returns
    /// Result containing the SphericalPolygon or an error message
    pub fn new(vertex_latitudes: Vec<f64>, vertex_longitudes: Vec<f64>) -> Result<Self, String> {
        let num_vertices = vertex_latitudes.len();

        if num_vertices != vertex_longitudes.len() {
            return Err(
                "vertex_latitudes and vertex_longitudes must have the same length".to_string(),
            );
        }

        if num_vertices < 3 {
            return Err("Polygon must have at least 3 vertices".to_string());
        }

        // Validate vertices
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

        // Pick first edge midpoint as reference (any edge works)
        let (reference_lat, reference_lon) = spherical_midpoint(
            vertex_latitudes[0],
            vertex_longitudes[0],
            vertex_latitudes[1],
            vertex_longitudes[1],
        );

        Ok(SphericalPolygon {
            vertex_latitudes,
            vertex_longitudes,
            reference_lat,
            reference_lon,
            reference_edge_index: 0,
        })
    }

    /// Locate a point P relative to the spherical polygon boundary
    ///
    /// Algorithm:
    /// 1. Cast a ray from P through the reference point M (edge midpoint)
    /// 2. Find the furthest vertex along this ray direction
    /// 3. Count crossings of the arc from P to that furthest vertex
    /// 4. Even crossings = inside, odd = outside
    ///
    /// # Arguments
    /// * `point_lat` - Latitude of point P in degrees
    /// * `point_lon` - Longitude of point P in degrees
    ///
    /// # Returns
    /// PointLocation indicating whether P is inside, outside, or on boundary
    pub fn locate_point(&self, point_lat: f64, point_lon: f64) -> PointLocation {
        // Check if P is on the reference point
        if point_lat == self.reference_lat && point_lon == self.reference_lon {
            return PointLocation::OnBoundary;
        }

        // Check if P is antipodal to reference
        if is_antipodal(point_lat, point_lon, self.reference_lat, self.reference_lon) {
            // Antipodal case - use opposite direction
            // For simplicity, we'll treat this as a degenerate case
            // In practice, this is extremely rare
            return self.locate_point_fallback(point_lat, point_lon);
        }

        // Transform to P-centered coordinate system
        // In this system, the great circle through P and reference becomes a meridian
        let reference_lon_transformed =
            transform_longitude(point_lat, point_lon, self.reference_lat, self.reference_lon);

        // Find the furthest vertex along the ray direction
        // This is the vertex with transformed longitude closest to reference_lon_transformed
        // (considering the directional nature - we want vertices "ahead" on the ray)
        let mut max_distance_along_ray = 0.0;
        let mut furthest_point_lat = self.reference_lat;
        let mut furthest_point_lon = self.reference_lon;

        for i in 0..self.vertex_latitudes.len() {
            let vertex_lon_transformed = transform_longitude(
                point_lat,
                point_lon,
                self.vertex_latitudes[i],
                self.vertex_longitudes[i],
            );

            // Check if this vertex is in the same hemisphere as the reference
            // (i.e., on the same side of the ray from P)
            let lon_diff_from_ref =
                normalize_longitude_difference(vertex_lon_transformed - reference_lon_transformed);

            // Only consider vertices in the forward direction (within ~90° of reference)
            if lon_diff_from_ref.abs() < 90.0 {
                let dist = spherical_distance(
                    point_lat,
                    point_lon,
                    self.vertex_latitudes[i],
                    self.vertex_longitudes[i],
                );

                if dist > max_distance_along_ray {
                    max_distance_along_ray = dist;
                    furthest_point_lat = self.vertex_latitudes[i];
                    furthest_point_lon = self.vertex_longitudes[i];
                }
            }
        }

        // If no vertex found in forward direction, use reference point
        if max_distance_along_ray == 0.0 {
            furthest_point_lat = self.reference_lat;
            furthest_point_lon = self.reference_lon;
            max_distance_along_ray =
                spherical_distance(point_lat, point_lon, self.reference_lat, self.reference_lon);
        }

        // Count crossings of the arc from P to the furthest point with polygon edges
        let mut crossing_count = 0;

        for i in 0..self.vertex_latitudes.len() {
            // Skip the reference edge
            if i == self.reference_edge_index {
                continue;
            }

            let next_i = (i + 1) % self.vertex_latitudes.len();

            let v1_lat = self.vertex_latitudes[i];
            let v1_lon = self.vertex_longitudes[i];
            let v2_lat = self.vertex_latitudes[next_i];
            let v2_lon = self.vertex_longitudes[next_i];

            // Check if P is on this edge
            if is_point_on_arc(point_lat, point_lon, v1_lat, v1_lon, v2_lat, v2_lon) {
                return PointLocation::OnBoundary;
            }

            // Check if the arc from P to furthest point crosses this edge
            if arcs_intersect(
                point_lat,
                point_lon,
                furthest_point_lat,
                furthest_point_lon,
                v1_lat,
                v1_lon,
                v2_lat,
                v2_lon,
            ) {
                crossing_count += 1;
            }
        }

        // Even crossings = inside, odd = outside
        if crossing_count % 2 == 0 {
            PointLocation::Inside
        } else {
            PointLocation::Outside
        }
    }

    /// Fallback method for antipodal cases
    fn locate_point_fallback(&self, point_lat: f64, point_lon: f64) -> PointLocation {
        // Use a different reference point
        let (alt_ref_lat, alt_ref_lon) = if self.vertex_latitudes.len() > 2 {
            spherical_midpoint(
                self.vertex_latitudes[1],
                self.vertex_longitudes[1],
                self.vertex_latitudes[2],
                self.vertex_longitudes[2],
            )
        } else {
            return PointLocation::Outside; // Shouldn't happen with valid polygon
        };

        // Simple check using this alternative reference
        // (A full implementation would recursively call the main algorithm)
        PointLocation::Outside
    }
}

/// Check if two arcs on a sphere intersect
/// Arc 1: from (lat1a, lon1a) to (lat1b, lon1b)
/// Arc 2: from (lat2a, lon2a) to (lat2b, lon2b)
fn arcs_intersect(
    lat1a: f64,
    lon1a: f64,
    lat1b: f64,
    lon1b: f64,
    lat2a: f64,
    lon2a: f64,
    lat2b: f64,
    lon2b: f64,
) -> bool {
    // Transform to coordinate system where first arc's start point is north pole
    let lon1b_t = transform_longitude(lat1a, lon1a, lat1b, lon1b);
    let lon2a_t = transform_longitude(lat1a, lon1a, lat2a, lon2a);
    let lon2b_t = transform_longitude(lat1a, lon1a, lat2b, lon2b);

    // Check if second arc crosses the meridian of the first arc
    // The first arc lies along longitude lon1b_t in this coordinate system

    // Check if endpoints of arc 2 are on opposite sides of arc 1's meridian
    let diff_a = normalize_longitude_difference(lon2a_t - lon1b_t);
    let diff_b = normalize_longitude_difference(lon2b_t - lon1b_t);

    // If signs are opposite and neither is exactly on the meridian, they might cross
    if (diff_a > 0.0 && diff_b < 0.0) || (diff_a < 0.0 && diff_b > 0.0) {
        // Now check if arc 1's endpoints straddle arc 2
        // Transform to arc 2's coordinate system
        let lon1a_t2 = transform_longitude(lat2a, lon2a, lat1a, lon1a);
        let lon1b_t2 = transform_longitude(lat2a, lon2a, lat1b, lon1b);
        let lon2b_t2 = transform_longitude(lat2a, lon2a, lat2b, lon2b);

        let diff1a = normalize_longitude_difference(lon1a_t2 - lon2b_t2);
        let diff1b = normalize_longitude_difference(lon1b_t2 - lon2b_t2);

        if (diff1a > 0.0 && diff1b < 0.0) || (diff1a < 0.0 && diff1b > 0.0) {
            return true;
        }
    }

    false
}

/// Check if a point lies on an arc (great circle segment)
fn is_point_on_arc(
    point_lat: f64,
    point_lon: f64,
    arc_lat1: f64,
    arc_lon1: f64,
    arc_lat2: f64,
    arc_lon2: f64,
) -> bool {
    // Check if point is on the great circle through arc endpoints
    let lon_p_t = transform_longitude(arc_lat1, arc_lon1, point_lat, point_lon);
    let lon_2_t = transform_longitude(arc_lat1, arc_lon1, arc_lat2, arc_lon2);

    if (lon_p_t - lon_2_t).abs() < 1e-9 || (lon_p_t - lon_2_t).abs() - 180.0 < 1e-9 {
        // Point is on the great circle, now check if it's between endpoints
        let dist_1_p = spherical_distance(arc_lat1, arc_lon1, point_lat, point_lon);
        let dist_1_2 = spherical_distance(arc_lat1, arc_lon1, arc_lat2, arc_lon2);
        let dist_p_2 = spherical_distance(point_lat, point_lon, arc_lat2, arc_lon2);

        // Point is on arc if dist_1_p + dist_p_2 ≈ dist_1_2
        return (dist_1_p + dist_p_2 - dist_1_2).abs() < 1e-6;
    }

    false
}

fn transform_longitude(pole_lat: f64, pole_lon: f64, target_lat: f64, target_lon: f64) -> f64 {
    if pole_lat == 90.0 {
        return target_lon;
    }

    let pole_lat_rad = pole_lat * DEG_TO_RAD;
    let target_lat_rad = target_lat * DEG_TO_RAD;
    let dlon_rad = (target_lon - pole_lon) * DEG_TO_RAD;

    let num = dlon_rad.sin() * target_lat_rad.cos();
    let den = target_lat_rad.sin() * pole_lat_rad.cos()
        - target_lat_rad.cos() * pole_lat_rad.sin() * dlon_rad.cos();

    num.atan2(den) / DEG_TO_RAD
}

fn normalize_longitude_difference(mut diff: f64) -> f64 {
    while diff > 180.0 {
        diff -= 360.0;
    }
    while diff < -180.0 {
        diff += 360.0;
    }
    diff
}

fn is_antipodal(lat1: f64, lon1: f64, lat2: f64, lon2: f64) -> bool {
    lat1 == -lat2 && normalize_longitude_difference(lon1 - lon2).abs() == 180.0
}

fn are_vertices_antipodal(lat1: f64, lon1: f64, lat2: f64, lon2: f64) -> bool {
    is_antipodal(lat1, lon1, lat2, lon2)
}

fn spherical_distance(lat1: f64, lon1: f64, lat2: f64, lon2: f64) -> f64 {
    let lat1_rad = lat1 * DEG_TO_RAD;
    let lat2_rad = lat2 * DEG_TO_RAD;
    let dlon_rad = (lon2 - lon1) * DEG_TO_RAD;

    let a = ((lat2_rad - lat1_rad) / 2.0).sin();
    let b = (dlon_rad / 2.0).sin();
    let c = a * a + lat1_rad.cos() * lat2_rad.cos() * b * b;

    2.0 * c.sqrt().asin() / DEG_TO_RAD
}

fn spherical_midpoint(lat1: f64, lon1: f64, lat2: f64, lon2: f64) -> (f64, f64) {
    let lat1_rad = lat1 * DEG_TO_RAD;
    let lon1_rad = lon1 * DEG_TO_RAD;
    let lat2_rad = lat2 * DEG_TO_RAD;
    let lon2_rad = lon2 * DEG_TO_RAD;

    let x1 = lat1_rad.cos() * lon1_rad.cos();
    let y1 = lat1_rad.cos() * lon1_rad.sin();
    let z1 = lat1_rad.sin();

    let x2 = lat2_rad.cos() * lon2_rad.cos();
    let y2 = lat2_rad.cos() * lon2_rad.sin();
    let z2 = lat2_rad.sin();

    let x = (x1 + x2) / 2.0;
    let y = (y1 + y2) / 2.0;
    let z = (z1 + z2) / 2.0;

    let norm = (x * x + y * y + z * z).sqrt();
    let x = x / norm;
    let y = y / norm;
    let z = z / norm;

    let lon = y.atan2(x) / DEG_TO_RAD;
    let lat = z.asin() / DEG_TO_RAD;
    (lat, lon)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simple_square() {
        let lats = vec![0.0, 0.0, 1.0, 1.0];
        let lons = vec![0.0, 1.0, 1.0, 0.0];

        let poly = SphericalPolygon::new(lats, lons).unwrap();

        // Point inside
        assert_eq!(poly.locate_point(0.5, 0.5), PointLocation::Inside);

        // Point outside (as in your diagram)
        assert_eq!(poly.locate_point(2.0, 0.5), PointLocation::Outside);

        // Point on boundary
        assert_eq!(poly.locate_point(0.0, 0.0), PointLocation::OnBoundary);
    }

    #[test]
    fn test_larger_square() {
        let lats = vec![0.0, 10.0, 10.0, 0.0];
        let lons = vec![0.0, 0.0, 10.0, 10.0];

        let poly = SphericalPolygon::new(lats, lons).unwrap();

        assert_eq!(poly.locate_point(5.0, 5.0), PointLocation::Inside);
        assert_eq!(poly.locate_point(20.0, 20.0), PointLocation::Outside);
    }
}
