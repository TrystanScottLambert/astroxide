use kiddo::ImmutableKdTree;

use crate::spherical_trig::{
    Point, angular_separation, build_kd_tree, convert_equitorial_to_cartesian, find_idx_within,
    spherical_mean,
};
fn combine_bool_vecs(list_1: Vec<bool>, list_2: Vec<bool>) -> Vec<bool> {
    list_1
        .iter()
        .zip(list_2.iter())
        .map(|(&a, &b)| a || b)
        .collect()
}
pub trait SphericalShape {
    fn is_inside(&self, ra: f64, dec: f64) -> bool;
    fn are_inside(&self, ras: &[f64], decs: &[f64]) -> Vec<bool>;
    fn are_inside_tree(
        &self,
        tree: &ImmutableKdTree<f64, 3>,
        ras: &[f64],
        decs: &[f64],
    ) -> Vec<bool>;
}

pub struct MaskingCatalog<'a> {
    ra_points: &'a [f64],
    dec_points: &'a [f64],
    kd_tree: ImmutableKdTree<f64, 3>,
    regions: &'a [Box<dyn SphericalShape>],
}
impl<'a> MaskingCatalog<'a> {
    pub fn new(
        ra_points: &'a [f64],
        dec_points: &'a [f64],
        regions: &'a [Box<dyn SphericalShape>],
    ) -> Self {
        let kd_tree = build_kd_tree(ra_points, dec_points);
        MaskingCatalog {
            ra_points,
            dec_points,
            kd_tree,
            regions,
        }
    }
    pub fn are_in_regions(&self) -> Vec<bool> {
        let mut current_results = vec![false; self.kd_tree.size()];
        for region in self.regions.iter() {
            let results = region.are_inside_tree(&self.kd_tree, self.ra_points, self.dec_points);
            current_results = combine_bool_vecs(results, current_results);
        }
        current_results
    }
}

pub struct SphericalAperture {
    ra_center: f64,
    dec_center: f64,
    radius_degrees: f64,
}
impl SphericalShape for SphericalAperture {
    fn is_inside(&self, ra: f64, dec: f64) -> bool {
        angular_separation(self.ra_center, self.dec_center, ra, dec) <= self.radius_degrees
    }

    fn are_inside_tree(
        &self,
        tree: &ImmutableKdTree<f64, 3>,
        _ras: &[f64],
        _decs: &[f64],
    ) -> Vec<bool> {
        let point = Point {
            ra_deg: self.ra_center,
            dec_deg: self.dec_center,
        };
        let idx = find_idx_within(&tree, &point, self.radius_degrees);
        let mut results = vec![false; tree.size()];
        for id in idx {
            results[id as usize] = true;
        }
        results
    }
    fn are_inside(&self, ras: &[f64], decs: &[f64]) -> Vec<bool> {
        let tree = build_kd_tree(ras, decs);
        self.are_inside_tree(&tree, ras, decs)
    }
}
impl SphericalAperture {
    pub fn new(ra_center: f64, dec_center: f64, radius_degrees: f64) -> Self {
        SphericalAperture {
            ra_center,
            dec_center,
            radius_degrees,
        }
    }
}

pub struct SphericalAnulus {
    ra_center: f64,
    dec_center: f64,
    inner_radius_deg: f64,
    outer_radius_deg: f64,
}

impl SphericalShape for SphericalAnulus {
    fn is_inside(&self, ra: f64, dec: f64) -> bool {
        let sep = angular_separation(self.ra_center, self.dec_center, ra, dec);
        sep <= self.outer_radius_deg && sep >= self.inner_radius_deg
    }

    fn are_inside_tree(
        &self,
        tree: &ImmutableKdTree<f64, 3>,
        _ras: &[f64],
        _decs: &[f64],
    ) -> Vec<bool> {
        let point = Point {
            ra_deg: self.ra_center,
            dec_deg: self.dec_center,
        };
        let idx_inner = find_idx_within(&tree, &point, self.inner_radius_deg - f64::EPSILON);
        let idx_outer = find_idx_within(&tree, &point, self.outer_radius_deg + 2. * f64::EPSILON);
        let mut results = vec![false; tree.size()];
        for id in idx_outer {
            results[id as usize] = true;
        }
        for id in idx_inner {
            results[id as usize] = false;
        }
        results
    }
    fn are_inside(&self, ras: &[f64], decs: &[f64]) -> Vec<bool> {
        let tree = build_kd_tree(ras, decs);
        self.are_inside_tree(&tree, ras, decs)
    }
}

impl SphericalAnulus {
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
}

/// Represents a spherical polygon boundary
pub struct SphericalPolygon {
    ra_verticies: Vec<f64>,
    dec_verticies: Vec<f64>,
    center: (f64, f64),
    bounding_radius: f64,
}

impl SphericalShape for SphericalPolygon {
    /// Point-in-polygon using the **spherical winding number** algorithm.
    ///
    /// Returns:
    /// - Inside
    /// - Outside
    /// - OnBoundary (if exactly on an edge)
    fn is_inside(&self, point_ra: f64, point_dec: f64) -> bool {
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
        total_angle.abs() > std::f64::consts::PI
    }

    fn are_inside_tree(
        &self,
        tree: &ImmutableKdTree<f64, 3>,
        ras: &[f64],
        decs: &[f64],
    ) -> Vec<bool> {
        let point = Point {
            ra_deg: self.center.0,
            dec_deg: self.center.1,
        };
        let idx = find_idx_within(&tree, &point, self.bounding_radius);
        let mut results = vec![false; tree.size()];
        for id in idx {
            results[id as usize] = self.is_inside(ras[id as usize], decs[id as usize]);
        }
        results
    }
    fn are_inside(&self, ra_points: &[f64], dec_points: &[f64]) -> Vec<bool> {
        let tree = build_kd_tree(ra_points, dec_points);
        self.are_inside_tree(&tree, ra_points, dec_points)
    }
}

impl SphericalPolygon {
    pub fn new(ra_verticies: Vec<f64>, dec_verticies: Vec<f64>) -> Self {
        let center = spherical_mean(&ra_verticies, &dec_verticies);
        let distances: Vec<f64> = ra_verticies
            .iter()
            .zip(dec_verticies.iter())
            .map(|(&ra, &dec)| angular_separation(center.0, center.1, ra, dec))
            .collect();
        let bounding_radius = distances.into_iter().max_by(|a, b| a.total_cmp(b)).unwrap() + 0.01;
        SphericalPolygon {
            ra_verticies,
            dec_verticies,
            center,
            bounding_radius,
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

#[cfg(test)]
mod tests {
    use std::iter::zip;

    use super::*;

    #[test]
    fn test_simple_square() {
        let ras = vec![0.0, 0.0, 1.0, 1.0];
        let decs = vec![0.0, 1.0, 1.0, 0.0];

        let poly = SphericalPolygon::new(ras, decs);

        assert!(poly.is_inside(0.5, 0.5));
        assert!(!poly.is_inside(2.0, 0.5));
    }

    #[test]
    fn test_square_at_ra_edge() {
        let ras = vec![359., 359., 1., 1.];
        let decs = vec![80., 82., 82., 80.];

        let poly = SphericalPolygon::new(ras, decs);

        assert!(poly.is_inside(0., 81.));
        assert!(!poly.is_inside(358., 81.));
    }

    #[test]
    fn test_concave_square() {
        let ras = vec![359., 359., 0., 1., 1.];
        let decs = vec![0., -2., -0.5, -2., 0.];

        let poly = SphericalPolygon::new(ras, decs);

        assert!(poly.is_inside(0., -0.3));
        assert!(poly.is_inside(359.5, -0.3));
        assert!(!poly.is_inside(0., -0.6));
    }

    #[test]
    fn test_multiple_at_pole() {
        let ra_verticies = vec![342.1537, 56.491250, 161.48667, 249.3462];
        let dec_verticies = vec![-81.250333, -74.158250, -80.6706, -78.951];
        let poly = SphericalPolygon::new(ra_verticies, dec_verticies);
        let ra_points = vec![18., 270., 133.11, 133.11];
        let dec_points = vec![-90., -90., -85.755, -60.];

        let answers = vec![true, true, true, false];
        let results = poly.are_inside(&ra_points, &dec_points);
        for (r, a) in zip(results, answers) {
            dbg!(&r);
            assert_eq!(r, a)
        }
    }

    #[test]
    fn test_aperture() {
        let aperture = SphericalAperture::new(0., 0., 1.);
        let ras = vec![0., 0., -0.1, 1., 2.];
        let decs = vec![0., 0.5, -0.1, 0., 2.];
        assert!(aperture.is_inside(0., 1.));
        assert!(aperture.is_inside(1., 0.));
        let results = aperture.are_inside(&ras, &decs);
        let answers = vec![true, true, true, true, false];
        for (r, a) in zip(results, answers) {
            assert_eq!(r, a)
        }
        assert!(aperture.is_inside(0.5, 0.5));
        assert!(aperture.is_inside(-0.5, -0.5));
        assert!(!aperture.is_inside(-2., -0.5));
    }

    #[test]
    fn test_anulus() {
        let anulus = SphericalAnulus::new(0., 0., 1., 2.);
        dbg!(anulus.is_inside(2., 0.));
        let ras = vec![0., 0., -0.1, 1., 2., 1.1, 0.];
        let decs = vec![0., 0.5, -0.1, 0., 2., 1.1, 2.];
        let results = anulus.are_inside(&ras, &decs);
        let answers = vec![false, false, false, true, false, true, true];
        for (r, a) in zip(results, answers) {
            dbg!(&r);
            assert_eq!(r, a)
        }
        assert!(!anulus.is_inside(0.5, 0.5));
        assert!(!anulus.is_inside(-0.5, -0.5));
        assert!(anulus.is_inside(1.2, 0.));
        assert!(!anulus.is_inside(2.2, 0.));
    }
    #[test]
    fn test_whole_catalog() {
        let app = SphericalAperture::new(0., 0., 1.);
        let anulus = SphericalAnulus::new(0., -90., 1., 2.);
        let poly = SphericalPolygon::new(vec![180., 180., 181., 181.], vec![0., 1., 1., 0.]);
        let regions: Vec<Box<dyn SphericalShape>> =
            vec![Box::new(app), Box::new(anulus), Box::new(poly)];
        let ra_points = vec![0., 0., 180.5, 0., 0., 179.];
        let dec_points = vec![0., -88.2, 0.5, 2., -89.5, 0.5];
        let answers = vec![true, true, true, false, false, false];
        let cat = MaskingCatalog::new(&ra_points, &dec_points, &regions);
        let results = cat.are_in_regions();
        for (r, a) in zip(results, answers) {
            dbg!(r);
            assert_eq!(r, a)
        }
    }
}
