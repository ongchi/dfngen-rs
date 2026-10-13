use parry3d_f64::na::{distance, Point3, Rotation3, Unit, Vector3};
use parry3d_f64::query::details::closest_points_segment_segment_with_locations_nD;
use parry3d_f64::query::PointQuery;
use parry3d_f64::shape::Segment;

// Check if two vectors are parallel
pub fn is_parallel(v1: &Vector3<f64>, v2: &Vector3<f64>, eps: f64) -> bool {
    let dot_prod = v1.dot(v2);
    1. - eps < dot_prod && dot_prod < 1. + eps
}

/// Upward-pointing normal vector of a plane given by the trend and plunge of its pole.
///
/// Angles are in radians. Trend is clockwise from North (+y), plunge is positive downward.
/// x is East, y is North, z is up.
pub fn normal_from_trend_plunge(trend: f64, plunge: f64) -> Vector3<f64> {
    Vector3::new(
        -trend.sin() * plunge.cos(),
        -trend.cos() * plunge.cos(),
        plunge.sin(),
    )
}

/// Normal vector of a plane given by its dip and strike (right-hand rule).
///
/// Angles are in radians. Strike is clockwise from North (+y).
/// x is East, y is North, z is up.
pub fn normal_from_dip_strike(dip: f64, strike: f64) -> Vector3<f64> {
    Vector3::new(
        dip.sin() * strike.cos(),
        -dip.sin() * strike.sin(),
        dip.cos(),
    )
}

/// Rotation which turns the unit vector `from` to the unit vector `to`.
pub fn rotation_between(from: &Vector3<f64>, to: &Vector3<f64>) -> Rotation3<f64> {
    Rotation3::rotation_between(from, to).unwrap_or_else(|| {
        // Opposite vectors: half turn about an axis perpendicular to them
        let not_parallel = if from.x.abs() < 0.9 {
            Vector3::x()
        } else {
            Vector3::y()
        };
        let axis = Unit::new_normalize(from.cross(&not_parallel));
        Rotation3::from_axis_angle(&axis, std::f64::consts::PI)
    })
}

/// True if each coordinate of the two points differs by less than `eps`.
pub(crate) fn is_same_point(point1: &Point3<f64>, point2: &Point3<f64>, eps: f64) -> bool {
    (point1.x - point2.x).abs() < eps
        && (point1.y - point2.y).abs() < eps
        && (point1.z - point2.z).abs() < eps
}

/// Shortest Disatance, point to line seg
///
/// # Arguments
///
/// * `point` - Point in 3d space
/// * `line` - Line defined by two end points
///
/// # Returns
///
/// * `f64` - Shortest distance between point and line segment
pub(crate) fn point_to_line_seg(point: &Point3<f64>, line: &[Point3<f64>; 2], eps: f64) -> f64 {
    if (line[0] - line[1]).magnitude_squared() < eps {
        // Line endpoints are equal to each other
        return distance(point, &line[0]);
    }

    Segment::new(line[0], line[1]).distance_to_local_point(point, true)
}

/// Fast check for two line segments which are far from each other
///
/// The distance between the bounding boxes of two segments is not larger than
/// the distance between the segments. Thus this check finds most of the far
/// pairs without the calculation of the exact distance.
///
/// # Arguments
///
/// * `line1` - Line 1 defined by two end points
/// * `line2` - Line 2 defined by two end points
/// * `limit` - Distance limit
///
/// # Returns
///
/// true - The distance between the segments is larger than `limit`
/// false - The distance can be larger or smaller than `limit`
pub(crate) fn are_segments_farther_than(
    line1: &[Point3<f64>; 2],
    line2: &[Point3<f64>; 2],
    limit: f64,
) -> bool {
    // The margin makes the check safe against rounding of the exact distance
    const RELATIVE_MARGIN: f64 = 1e-6;

    let mut squared_box_distance = 0.;

    for axis in 0..3 {
        let (min1, max1) = (
            line1[0][axis].min(line1[1][axis]),
            line1[0][axis].max(line1[1][axis]),
        );
        let (min2, max2) = (
            line2[0][axis].min(line2[1][axis]),
            line2[0][axis].max(line2[1][axis]),
        );
        let gap = (min1 - max2).max(min2 - max1).max(0.);
        squared_box_distance += gap * gap;
    }

    let limit_with_margin = limit * (1. + RELATIVE_MARGIN);
    squared_box_distance > limit_with_margin * limit_with_margin
}

/// Closest Distance from Line Seg to Line Seg
///
/// # Arguments
///
/// * `line1` - Line 1 defined by two end points
/// * `line2` - Line 2 defined by two end points
///
/// # Returns
///
/// The minimum distance between line 1 and line 2, and the point on line 1
/// which is closest to line 2. If the lines intersect, this is the intersection point.
pub(crate) fn line_seg_to_line_seg(
    line1: &[Point3<f64>; 2],
    line2: &[Point3<f64>; 2],
) -> (f64, Point3<f64>) {
    let (location1, location2) = closest_points_segment_segment_with_locations_nD(
        (&line1[0], &line1[1]),
        (&line2[0], &line2[1]),
    );
    let point1 = Segment::new(line1[0], line1[1]).point_at(&location1);
    let point2 = Segment::new(line2[0], line2[1]).point_at(&location2);

    (distance(&point1, &point2), point1)
}

#[cfg(test)]
mod tests {
    use std::f64::consts::FRAC_PI_2;

    use super::*;

    fn assert_vector_eq(actual: Vector3<f64>, expected: [f64; 3]) {
        assert!(
            (actual - Vector3::from(expected)).magnitude() < 1e-12,
            "expected {expected:?}, found {:?}",
            actual.as_slice()
        );
    }

    /// A segment of length 2 along x, and a segment along y at the given
    /// distance from it in z
    fn crossed_segments(z_distance: f64) -> ([Point3<f64>; 2], [Point3<f64>; 2]) {
        (
            [Point3::new(-1., 0., 0.), Point3::new(1., 0., 0.)],
            [
                Point3::new(0., -1., z_distance),
                Point3::new(0., 1., z_distance),
            ],
        )
    }

    #[test]
    fn test_segments_with_a_box_distance_above_the_limit_are_far() {
        let (line1, line2) = crossed_segments(2.);

        assert!(are_segments_farther_than(&line1, &line2, 1.));
    }

    #[test]
    fn test_segments_with_a_distance_below_the_limit_are_not_far() {
        let (line1, line2) = crossed_segments(0.5);

        assert!(!are_segments_farther_than(&line1, &line2, 1.));
    }

    #[test]
    fn test_segments_with_a_distance_equal_to_the_limit_are_not_far() {
        let (line1, line2) = crossed_segments(1.);

        assert!(!are_segments_farther_than(&line1, &line2, 1.));
    }

    #[test]
    fn test_far_segments_with_touching_boxes_are_not_found_as_far() {
        // The distance is sqrt(2), but the bounding boxes have the corner (1, 1, 0) in common
        let line1 = [Point3::new(0., 0., 0.), Point3::new(1., 1., 0.)];
        let line2 = [Point3::new(1., 3., 0.), Point3::new(3., 1., 0.)];

        assert!(!are_segments_farther_than(&line1, &line2, 1.));
    }

    #[test]
    fn test_far_check_agrees_with_the_exact_distance() {
        let (line1, line2) = crossed_segments(3.);
        let (exact_distance, _) = line_seg_to_line_seg(&line1, &line2);

        assert!(are_segments_farther_than(&line1, &line2, 2.) && exact_distance > 2.);
    }

    #[test]
    fn test_vertical_plunge_gives_upward_normal() {
        assert_vector_eq(normal_from_trend_plunge(0.3, FRAC_PI_2), [0., 0., 1.]);
    }

    #[test]
    fn test_horizontal_pole_with_trend_north_gives_normal_to_south() {
        assert_vector_eq(normal_from_trend_plunge(0., 0.), [0., -1., 0.]);
    }

    #[test]
    fn test_horizontal_pole_with_trend_east_gives_normal_to_west() {
        assert_vector_eq(normal_from_trend_plunge(FRAC_PI_2, 0.), [-1., 0., 0.]);
    }

    #[test]
    fn test_zero_dip_gives_upward_normal() {
        assert_vector_eq(normal_from_dip_strike(0., 1.1), [0., 0., 1.]);
    }

    #[test]
    fn test_vertical_plane_with_strike_north_gives_normal_to_east() {
        assert_vector_eq(normal_from_dip_strike(FRAC_PI_2, 0.), [1., 0., 0.]);
    }

    #[test]
    fn test_vertical_plane_with_strike_east_gives_normal_to_south() {
        assert_vector_eq(normal_from_dip_strike(FRAC_PI_2, FRAC_PI_2), [0., -1., 0.]);
    }

    #[test]
    fn test_normals_from_angles_have_unit_length() {
        assert!((normal_from_trend_plunge(0.7, 0.4).magnitude() - 1.).abs() < 1e-12);
        assert!((normal_from_dip_strike(0.7, 0.4).magnitude() - 1.).abs() < 1e-12);
    }
}
