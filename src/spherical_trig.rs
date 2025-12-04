use kiddo::ImmutableKdTree;
use kiddo::SquaredEuclidean;
use libm::{asin, atan2, hypot};

// assuming a unit sphere
pub fn convert_equitorial_to_cartesian(ra_deg: &f64, dec_deg: &f64) -> [f64; 3] {
    let ra_radians = ra_deg.to_radians();
    let dec_radians = dec_deg.to_radians();
    let x = dec_radians.cos() * ra_radians.cos();
    let y = dec_radians.cos() * ra_radians.sin();
    let z = dec_radians.sin();
    [x, y, z]
}

// assuming a unit sphere
pub fn convert_cartesian_to_equitorial(x: &f64, y: &f64, z: &f64) -> [f64; 2] {
    let radius = (x.powi(2) + y.powi(2) + z.powi(2)).sqrt();
    let ra_radian = atan2(*y, *x);
    let dec_radian = asin(z / radius);
    [ra_radian.to_degrees(), dec_radian.to_degrees()]
}

// Converts RA, Dec (degrees) and comoving distance to 3D Cartesian coordinates.
pub fn convert_equitorial_to_cartesian_scaled(
    ra_deg: f64,
    dec_deg: f64,
    distance: f64,
) -> [f64; 3] {
    let ra_rad = ra_deg.to_radians();
    let dec_rad = dec_deg.to_radians();
    let x = distance * dec_rad.cos() * ra_rad.cos();
    let y = distance * dec_rad.cos() * ra_rad.sin();
    let z = distance * dec_rad.sin();
    [x, y, z]
}

// The three dimensional euclidean distance
pub fn euclidean_distance_3d(point_1: &[f64; 3], point_2: &[f64; 3]) -> f64 {
    ((point_1[0] - point_2[0]).powi(2)
        + (point_1[1] - point_2[1]).powi(2)
        + (point_1[2] - point_2[2]).powi(2))
    .sqrt()
}

/// Chord length. Given the angular separation.
pub fn chord_distance(angular_separation_degrees: f64) -> f64 {
    assert!((0. ..=180.).contains(&angular_separation_degrees));
    2. * (angular_separation_degrees.to_radians() / 2.).sin()
}

/// Copy of the astropy angular_separation function.
///
/// Arguments must be in degrees and the separation is in degrees.
/// Using the vincenty formula which is accurate but computationally expensive.
pub fn angular_separation(lon_1: f64, lat_1: f64, lon_2: f64, lat_2: f64) -> f64 {
    let lon_1_rad = lon_1.to_radians();
    let lat_1_rad = lat_1.to_radians();
    let lon_2_rad = lon_2.to_radians();
    let lat_2_rad = lat_2.to_radians();
    let sdlon = (lon_2_rad - lon_1_rad).sin();
    let cdlon = (lon_2_rad - lon_1_rad).cos();
    let slat1 = (lat_1_rad).sin();
    let slat2 = (lat_2_rad).sin();
    let clat1 = (lat_1_rad).cos();
    let clat2 = (lat_2_rad).cos();

    let num1 = clat2 * sdlon;
    let num2 = clat1 * slat2 - slat1 * clat2 * cdlon;
    let denominator = slat1 * slat2 + clat1 * clat2 * cdlon;
    atan2(hypot(num1, num2), denominator).to_degrees()
}

/// A quick angular separation which works for small angles.
///
/// Angles must be in degrees and the separation is returned in degrees.
/// Only use for angles less than a degree.
///
pub fn angular_separation_haversine(
    lon1_deg: &f64,
    lat1_deg: &f64,
    lon2_deg: &f64,
    lat2_deg: &f64,
) -> f64 {
    // Convert to radians
    let lon1 = lon1_deg.to_radians();
    let lat1 = lat1_deg.to_radians();
    let lon2 = lon2_deg.to_radians();
    let lat2 = lat2_deg.to_radians();

    // Differences
    let dlat = lat2 - lat1;
    let dlon = lon2 - lon1;

    // Haversine formula
    let a = (dlat / 2.0).sin().powi(2) + lat1.cos() * lat2.cos() * (dlon / 2.0).sin().powi(2);
    let c = 2.0 * a.sqrt().asin();

    // Convert back to degrees
    c.to_degrees()
}

/// Very quick angular separation using the small angle formula.
///
/// Only use for very small angles. Inputs and outputs are in degrees.
pub fn angular_separation_small_angle(
    lon1_deg: &f64,
    lat1_deg: &f64,
    lon2_deg: &f64,
    lat2_deg: &f64,
) -> f64 {
    // Convert to radians
    let lon1 = lon1_deg.to_radians();
    let lat1 = lat1_deg.to_radians();
    let lon2 = lon2_deg.to_radians();
    let lat2 = lat2_deg.to_radians();

    // Differences
    let dlat = lat2 - lat1;
    let dlon = lon2 - lon1;

    // Haversine formula
    let a = (dlat / 2.0).sin().powi(2) + lat1.cos() * lat2.cos() * (dlon / 2.0).sin().powi(2);
    let c = 2.0 * a.sqrt().asin();

    // Convert back to degrees
    c.to_degrees()
}

#[derive(Debug)]
pub struct Point {
    pub ra_deg: f64,
    pub dec_deg: f64,
}

pub fn build_kd_tree(ra_array_deg: &[f64], dec_array_deg: &[f64]) -> ImmutableKdTree<f64, 3> {
    let entries: Vec<[f64; 3]> = ra_array_deg
        .iter()
        .zip(dec_array_deg)
        .map(|(ra, dec)| convert_equitorial_to_cartesian(ra, dec))
        .collect();

    ImmutableKdTree::new_from_slice(&entries)
}

pub fn find_idx_within(
    tree: &ImmutableKdTree<f64, 3>,
    point: &Point,
    angular_difference_deg: f64,
) -> Vec<u64> {
    let point_cartesian = convert_equitorial_to_cartesian(&point.ra_deg, &point.dec_deg);
    let dist = chord_distance(angular_difference_deg).powi(2); // squared chord distance.
    let neighbors = tree.within_unsorted::<SquaredEuclidean>(&point_cartesian, dist);
    neighbors.iter().map(|nn| nn.item).collect()
}

pub fn spherical_mean(ra_points: &[f64], dec_points: &[f64]) -> (f64, f64) {
    assert_eq!(ra_points.len(), dec_points.len());
    let mut sum = [0., 0., 0.];
    for (ra, dec) in ra_points.iter().zip(dec_points.iter()) {
        let v = convert_equitorial_to_cartesian(ra, dec);
        sum[0] += v[0];
        sum[1] += v[1];
        sum[2] += v[2];
    }

    let norm = (sum[0].powi(2) + sum[1].powi(2) + sum[2].powi(2)).sqrt();
    let x = sum[0] / norm;
    let y = sum[1] / norm;
    let z = sum[2] / norm;

    let center = convert_cartesian_to_equitorial(&x, &y, &z);

    let mut ra_center = center[0];
    if ra_center < 0.0 {
        ra_center += 360.;
    }
    (ra_center, center[1])
}

#[cfg(test)]
mod test {
    use super::*;
    use std::iter::zip;

    const EPSILON: f64 = 1e-6;

    #[test]
    fn test_kd_tree_building() {
        let ras = vec![0., 23., 180.];
        let decs = vec![0., -10., 20.];
        let tree = build_kd_tree(&ras, &decs);

        let test_point = convert_equitorial_to_cartesian(&ras[1], &decs[1]);
        let result = tree.within::<SquaredEuclidean>(&test_point, 1e-10);
        assert_eq!(result.len(), 1);
    }

    #[test]
    fn test_finding_within() {
        let ras = vec![120., 120.9999, 120., 180.];
        let decs = vec![0., 0., 0., -45.];
        let tree = build_kd_tree(&ras, &decs);
        let point = Point {
            ra_deg: 120.,
            dec_deg: 0.,
        };
        let mut result = find_idx_within(&tree, &point, 1.);
        assert_eq!(result.len(), 3);
        let answers = [0, 1, 2];
        result.sort();
        for (res, ans) in zip(result, answers) {
            assert_eq!(res, ans)
        }
    }

    #[test]
    fn test_finding_at_boundaries() {
        let ras = vec![0., 1., 0., 270.];
        let decs = vec![0., 0., 1., -45.];
        let tree = build_kd_tree(&ras, &decs);
        let point = Point {
            ra_deg: 0.,
            dec_deg: 0.,
        };
        let mut result = find_idx_within(&tree, &point, 1.);
        assert_eq!(result.len(), 3);
        let answers = [0, 1, 2];
        result.sort();
        for (res, ans) in zip(result, answers) {
            assert_eq!(res, ans)
        }
    }
    #[test]
    fn test_equatorial_to_cartesian_unit_sphere() {
        let ra = 0.0;
        let dec = 0.0;
        let result = convert_equitorial_to_cartesian(&ra, &dec);
        assert!((result[0] - 1.0).abs() < EPSILON);
        assert!((result[1]).abs() < EPSILON);
        assert!((result[2]).abs() < EPSILON);
    }

    #[test]
    fn test_cartesian_to_equatorial_unit_sphere() {
        let x = 1.0;
        let y = 0.0;
        let z = 0.0;
        let result = convert_cartesian_to_equitorial(&x, &y, &z);
        assert!((result[0] - 0.0).abs() < EPSILON); // RA
        assert!((result[1] - 0.0).abs() < EPSILON); // Dec
    }

    #[test]
    fn test_round_trip_conversion() {
        let ra = 123.4;
        let dec = -56.78;
        let cart = convert_equitorial_to_cartesian(&ra, &dec);
        let equi = convert_cartesian_to_equitorial(&cart[0], &cart[1], &cart[2]);

        // RA can wrap around so we normalize to [0, 360)
        let mut ra_norm = equi[0];
        if ra_norm < 0.0 {
            ra_norm += 360.0;
        }

        assert!((ra_norm - ra).abs() < EPSILON);
        assert!((equi[1] - dec).abs() < EPSILON);
    }

    #[test]
    fn test_equatorial_to_cartesian_scaled() {
        let ra = [120., 122., 124.];
        let dec = [-56., -34., -30.];
        let distance = [10., 11., 12.];
        let result_1 = convert_equitorial_to_cartesian_scaled(ra[0], dec[0], distance[0]);
        let result_2 = convert_equitorial_to_cartesian_scaled(ra[1], dec[1], distance[1]);
        let result_3 = convert_equitorial_to_cartesian_scaled(ra[2], dec[2], distance[2]);

        let answer_1 = [-2.795965, 4.842753, -8.290376];
        let answer_2 = [-4.832553, 7.733701, -6.151122];
        let answer_3 = [-5.811303, 8.615611, -6.000000];

        for (r, a) in zip(result_1, answer_1) {
            assert!((r - a).abs() < EPSILON);
        }

        for (r, a) in zip(result_2, answer_2) {
            assert!((r - a).abs() < EPSILON);
        }

        for (r, a) in zip(result_3, answer_3) {
            assert!((r - a).abs() < EPSILON);
        }
    }

    #[test]
    fn test_euclidean_distance_3d() {
        let a = [1.0, 0.0, 0.0];
        let b = [0.0, 0.0, 0.0];
        let dist = euclidean_distance_3d(&a, &b);
        assert!((dist - 1.0).abs() < EPSILON);

        let a = [1.0, 2.0, 3.0];
        let b = [4.0, 6.0, 3.0];
        let dist = euclidean_distance_3d(&a, &b);
        assert!((dist - 5.0).abs() < EPSILON);
    }

    #[test]
    fn test_chord_distance() {
        let angular_separation = [0., 20., 30., 180.];
        let results = angular_separation.iter().map(|&a| chord_distance(a));
        let answers = [0., 0.347296355, 0.517638090, 2.];
        for (res, ans) in zip(results, answers) {
            assert!((res - ans).abs() < 1e-7)
        }
    }

    #[test]
    fn test_projected_separation() {
        let ra_1s = [-20., -5., 0., 90., 359.];
        let ra_2s = [-20.5, 20., 30., 40., 50.];
        let dec_1s = [-90., -45., 0., 45., 90.];
        let dec_2s = [-45., 0., 5., 10., 20.];
        let answers = [45., 50.14429257, 30.37550653, 55.22172917, 70.];

        for i in 0..5 {
            let res = angular_separation(ra_1s[i], dec_1s[i], ra_2s[i], dec_2s[i]);
            dbg!(res);
            dbg!(answers[i]);
            assert!((res - answers[i]).abs() < 1e-5)
        }
    }
    #[test]
    fn test_separation_over_ra_edge() {
        let result = angular_separation(1., 0., 359., 0.);
        assert!((result - 2.).abs() < 1e-9);
    }

    #[test]
    fn test_small_angle_approx() {
        // testing that the small angle version is close to the accurate one for small angles.
        let ra1 = 0.1;
        let ra2 = 0.2;
        let dec1 = -20.;
        let dec2 = -20.;
        let res = angular_separation_haversine(&ra1, &dec1, &ra2, &dec2);
        dbg!(res);
        dbg!(angular_separation(ra1, dec1, ra2, dec2));
        assert!((angular_separation(ra1, dec1, ra2, dec2) - res).abs() < 1e-5);

        let res = angular_separation_small_angle(&ra1, &dec1, &ra2, &dec2);
        assert!((angular_separation(ra1, dec1, ra2, dec2) - res).abs() < 1e-5);
    }
}
