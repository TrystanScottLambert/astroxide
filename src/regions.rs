use crate::spherical_trig::{
    Point, angular_separation, build_kd_tree, convert_equitorial_to_cartesian, find_idx_within,
};

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
    pub fn locate_point(&self, ra: f64, dec: f64) -> PointLocation {
        let sep = angular_separation(self.ra_center, self.dec_center, ra, dec);
        if sep < self.radius_degrees {
            PointLocation::Inside
        } else if sep > self.radius_degrees {
            PointLocation::Outside
        } else {
            PointLocation::OnBoundary
        }
    }

    pub fn new(ra_center: f64, dec_center: f64, radius_degrees: f64) -> Self {
        SphericalAperture {
            ra_center,
            dec_center,
            radius_degrees,
        }
    }

    pub fn locate_points(&self, ras: &[f64], decs: &[f64]) -> Vec<PointLocation> {
        let tree = build_kd_tree(ras, decs);
        let point = Point {
            ra_deg: self.ra_center,
            dec_deg: self.dec_center,
        };
        let idx = find_idx_within(&tree, &point, self.radius_degrees);
        let mut results = vec![PointLocation::Outside; ras.len()];
        for id in idx {
            results[id as usize] = PointLocation::Inside;
        }
        results
        // TODO: include bounday condition?
    }
}

pub struct SphericalAnulus {
    ra_center: f64,
    dec_center: f64,
    inner_radius_deg: f64,
    outer_radius_deg: f64,
}

impl SphericalAnulus {
    pub fn locate_point(&self, ra: f64, dec: f64) -> PointLocation {
        let sep = angular_separation(self.ra_center, self.dec_center, ra, dec);
        if sep > self.outer_radius_deg || sep < self.inner_radius_deg {
            PointLocation::Outside
        } else if sep == self.outer_radius_deg || sep == self.inner_radius_deg {
            PointLocation::OnBoundary
        } else {
            PointLocation::Inside
        }
    }
    pub fn new(
        ra_center: f64,
        dec_center: f64,
        inner_radius_deg: f64,
        outer_radius_deg: f64,
    ) -> Self {
        SphericalAnulus {
            ra_center,
            dec_center,
            inner_radius_deg,
            outer_radius_deg,
        }
    }
    pub fn locate_points(&self, ras: &[f64], decs: &[f64]) -> Vec<PointLocation> {
        let tree = build_kd_tree(ras, decs);
        let point = Point {
            ra_deg: self.ra_center,
            dec_deg: self.dec_center,
        };
        let idx_inner = find_idx_within(&tree, &point, self.inner_radius_deg);
        let idx_outer = find_idx_within(&tree, &point, self.outer_radius_deg);
        let mut results = vec![PointLocation::Outside; ras.len()];
        for id in idx_outer {
            results[id as usize] = PointLocation::Inside;
        }
        for id in idx_inner {
            results[id as usize] = PointLocation::Outside;
        }
        results
    }
}

/// Represents a spherical polygon boundary
pub struct SphericalPolygon {
    ra_verticies: Vec<f64>,
    dec_verticies: Vec<f64>,
}

impl SphericalPolygon {
    pub fn new(ra_verticies: Vec<f64>, dec_verticies: Vec<f64>) -> Self {
        SphericalPolygon {
            ra_verticies,
            dec_verticies,
        }
    }
    /// Point-in-polygon using the **spherical winding number** algorithm.
    ///
    /// Returns:
    /// - Inside
    /// - Outside
    /// - OnBoundary (if exactly on an edge)
    pub fn locate_point(&self, point_ra: f64, point_dec: f64) -> PointLocation {
        // Convert P to a 3D unit vector
        let p = convert_equitorial_to_cartesian(&point_ra, &point_dec);

        // Build polygon as unit vectors
        let n = self.ra_verticies.len();
        let mut verts: Vec<[f64; 3]> = Vec::with_capacity(n);
        for (&ra, &dec) in self.ra_verticies.iter().zip(self.dec_verticies.iter()) {
            verts.push(convert_equitorial_to_cartesian(&ra, &dec));
        }

        // Compute total signed winding angle around P
        let mut total_angle = 0.0;

        for i in 0..n {
            let a = verts[i];
            let b = verts[(i + 1) % n];

            // Check if P is exactly on the edge A-B
            if point_on_gc_segment(p, a, b) {
                return PointLocation::OnBoundary;
            }

            // Project both vertices into the tangent plane of P
            let ua = project_to_tangent(a, p);
            let ub = project_to_tangent(b, p);

            // Compute signed angle from ua → ub
            let cross_term = dot(p, cross(ua, ub));
            let dot_term = dot(ua, ub);
            let angle = cross_term.atan2(dot_term);

            total_angle += angle;
        }

        // Classification:
        // If magnitude of winding is ~2π → inside
        if total_angle.abs() > std::f64::consts::PI {
            PointLocation::Inside
        } else {
            PointLocation::Outside
        }
    }
}

#[inline]
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
fn normalize(v: [f64; 3]) -> [f64; 3] {
    let n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    [v[0] / n, v[1] / n, v[2] / n]
}

/// Project vector a into the tangent plane at p
#[inline]
fn project_to_tangent(a: [f64; 3], p: [f64; 3]) -> [f64; 3] {
    let ap = dot(a, p);
    normalize([a[0] - ap * p[0], a[1] - ap * p[1], a[2] - ap * p[2]])
}

fn point_on_gc_segment(p: [f64; 3], a: [f64; 3], b: [f64; 3]) -> bool {
    const EPS: f64 = 1e-12;

    // Check if P lies on the great circle of AB
    let cross_ab = cross(a, b);
    let d = dot(cross_ab, p).abs();
    if d > EPS {
        return false;
    }

    // Check if P lies *between* a and b (not beyond endpoints)
    // by verifying angle(AP) + angle(PB) ≈ angle(AB)
    let ang_ab = dot(a, b).acos();
    let ang_ap = dot(a, p).acos();
    let ang_pb = dot(p, b).acos();

    (ang_ap + ang_pb - ang_ab).abs() < 1e-10
}

#[cfg(test)]
mod tests {
    use std::iter::zip;

    use super::*;

    #[test]
    fn test_simple_square() {
        let ras = vec![0.0, 0.0, 1.0, 1.0];
        let decs = vec![0.0, 1.0, 1.0, 0.0];

        let poly = SphericalPolygon::new(ras, decs);

        // Point inside
        assert_eq!(poly.locate_point(0.5, 0.5), PointLocation::Inside);

        // Point outside (as in your diagram)
        assert_eq!(poly.locate_point(2.0, 0.5), PointLocation::Outside);

        // Point on boundary
        assert_eq!(poly.locate_point(0.0, 0.0), PointLocation::OnBoundary);
    }

    #[test]
    fn test_square_at_ra_edge() {
        let ras = vec![359., 359., 1., 1.];
        let decs = vec![80., 82., 82., 80.];

        let poly = SphericalPolygon::new(ras, decs);

        assert_eq!(poly.locate_point(0., 81.), PointLocation::Inside);
        assert_eq!(poly.locate_point(358., 81.), PointLocation::Outside);
    }

    #[test]
    fn test_aperture() {
        let aperture = SphericalAperture::new(0., 0., 1.);
        let ras = vec![0., 0., -0.1, 1., 2.];
        let decs = vec![0., 0.5, -0.1, 0., 2.];
        let results = aperture.locate_points(&ras, &decs);
        let answers = vec![
            PointLocation::Inside,
            PointLocation::Inside,
            PointLocation::Inside,
            PointLocation::Inside,
            PointLocation::Outside,
        ];
        for (r, a) in zip(results, answers) {
            assert_eq!(r, a)
        }
        assert_eq!(aperture.locate_point(0.5, 0.5), PointLocation::Inside);
        assert_eq!(aperture.locate_point(-0.5, -0.5), PointLocation::Inside);
        assert_eq!(aperture.locate_point(-2., -0.5), PointLocation::Outside);
        assert_eq!(aperture.locate_point(0., -1.), PointLocation::OnBoundary);
    }

    #[test]
    fn test_anulus() {
        let anulus = SphericalAnulus::new(0., 0., 1., 2.);
        let ras = vec![0., 0., -0.1, 1., 2., 1.1, 2.];
        let decs = vec![0., 0.5, -0.1, 0., 2., 1.1, 0.];
        let results = anulus.locate_points(&ras, &decs);
        let answers = vec![
            PointLocation::Outside,
            PointLocation::Outside,
            PointLocation::Outside,
            PointLocation::Outside,
            PointLocation::Outside,
            PointLocation::Inside,
            PointLocation::Outside,
        ];
        for (r, a) in zip(results, answers) {
            assert_eq!(r, a)
        }
        assert_eq!(anulus.locate_point(0.5, 0.5), PointLocation::Outside);
        assert_eq!(anulus.locate_point(-0.5, -0.5), PointLocation::Outside);
        assert_eq!(anulus.locate_point(-2., 0.), PointLocation::OnBoundary);
        assert_eq!(anulus.locate_point(1.2, 0.), PointLocation::Inside);
        assert_eq!(anulus.locate_point(2.2, 0.), PointLocation::Outside);
    }
}
