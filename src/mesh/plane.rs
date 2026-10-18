//! A fracture in its own plane: the outline, the lines which the mesh must
//! follow, and the target size of the triangles

use parry3d_f64::na::Point3;

use crate::fracture::poly::Poly;
use crate::geometry::PlaneFrame;
use crate::region_boundaries::{distance_to_segment, OtherSide};

/// A point in the plane of a fracture
pub type Point2 = [f64; 2];

/// A boundary line segment in the plane, with what is on its other side
pub struct BoundarySegment2 {
    pub segment: [Point2; 2],
    pub other_side: OtherSide,
}

/// A final fracture in its plane
pub struct FracturePlane {
    pub frame: PlaneFrame,
    /// Vertices of the fracture in the plane
    pub outline: Vec<Point2>,
    /// Intersections with the other final fractures
    pub intersections: Vec<[Point2; 2]>,
    pub boundary_segments: Vec<BoundarySegment2>,
    /// Triangle size near a line: `min_feature_size`
    pub min_size: f64,
    pub max_size: f64,
    pub slope: f64,
}

impl FracturePlane {
    /// # Arguments
    ///
    /// * `poly` - The fracture
    /// * `intersections` - End points of its intersections with the other final fractures
    /// * `sizes` - `min_feature_size`, the largest size and the slope
    pub fn new(poly: &Poly, intersections: &[[Point3<f64>; 2]], sizes: [f64; 3]) -> Self {
        let frame = PlaneFrame::new(&poly.vertices, &poly.normal);
        let [min_size, max_size, slope] = sizes;

        Self {
            outline: poly
                .vertices
                .iter()
                .map(|vertex| frame.to_plane(vertex))
                .collect(),
            intersections: intersections
                .iter()
                .map(|[start, end]| [frame.to_plane(start), frame.to_plane(end)])
                .collect(),
            boundary_segments: poly
                .boundary_lines
                .iter()
                .flat_map(|line| &line.segments)
                .map(|segment| BoundarySegment2 {
                    segment: [frame.to_plane(&segment.start), frame.to_plane(&segment.end)],
                    other_side: segment.other_side.clone(),
                })
                .collect(),
            frame,
            min_size,
            max_size,
            slope,
        }
    }

    /// The edges of the outline
    pub fn outline_edges(&self) -> impl Iterator<Item = [Point2; 2]> + '_ {
        (0..self.outline.len())
            .map(|i| [self.outline[i], self.outline[(i + 1) % self.outline.len()]])
    }

    /// The lines near which the triangles have the smallest size
    fn feature_segments(&self) -> impl Iterator<Item = [Point2; 2]> + '_ {
        self.intersections.iter().copied().chain(
            self.boundary_segments
                .iter()
                .map(|boundary| boundary.segment),
        )
    }

    /// Target triangle size at the point: the smallest size plus the slope
    /// times the distance to the nearest intersection or boundary line, at
    /// most the largest size
    pub fn size_at(&self, point: Point2) -> f64 {
        let distance = self
            .feature_segments()
            .map(|segment| distance_to_segment(point, segment))
            .fold(f64::INFINITY, f64::min);

        (self.min_size + self.slope * distance).min(self.max_size)
    }

    /// True if the point is in the outline (even-odd rule)
    pub fn contains(&self, [x, y]: Point2) -> bool {
        self.outline_edges()
            .filter(|[start, end]| (start[1] > y) != (end[1] > y))
            .filter(|[start, end]| {
                let crossing_x =
                    start[0] + (y - start[1]) / (end[1] - start[1]) * (end[0] - start[0]);
                x < crossing_x
            })
            .count()
            % 2
            == 1
    }

    /// Distance to the nearest line which the mesh must follow: the outline,
    /// the intersections and the boundary lines
    pub fn distance_to_lines(&self, point: Point2) -> f64 {
        self.outline_edges()
            .chain(self.feature_segments())
            .map(|segment| distance_to_segment(point, segment))
            .fold(f64::INFINITY, f64::min)
    }
}

pub fn point_on_segment([start, end]: [Point2; 2], position: f64) -> Point2 {
    [
        start[0] + (end[0] - start[0]) * position,
        start[1] + (end[1] - start[1]) * position,
    ]
}

/// Sine of the angle below which two segments are parallel
const PARALLEL_SINE: f64 = 1e-9;

/// Positions on the two segments (0 to 1) of their crossing point. `None` if
/// they are parallel or do not cross.
pub fn segment_crossing(first: [Point2; 2], second: [Point2; 2]) -> Option<(f64, f64)> {
    let first_direction = [first[1][0] - first[0][0], first[1][1] - first[0][1]];
    let second_direction = [second[1][0] - second[0][0], second[1][1] - second[0][1]];
    let denominator =
        first_direction[0] * second_direction[1] - first_direction[1] * second_direction[0];

    // Almost parallel: the crossing position has a large error. A line which
    // touches the other line is found by its ends.
    let lengths = first_direction[0].hypot(first_direction[1])
        * second_direction[0].hypot(second_direction[1]);
    if denominator.abs() <= PARALLEL_SINE * lengths {
        return None;
    }

    let offset = [second[0][0] - first[0][0], second[0][1] - first[0][1]];
    let first_position =
        (offset[0] * second_direction[1] - offset[1] * second_direction[0]) / denominator;
    let second_position =
        (offset[0] * first_direction[1] - offset[1] * first_direction[0]) / denominator;

    ((0. ..=1.).contains(&first_position) && (0. ..=1.).contains(&second_position))
        .then_some((first_position, second_position))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_crossing_of_two_diagonals_is_at_their_middle() {
        let crossing = segment_crossing([[0., 0.], [2., 2.]], [[0., 2.], [2., 0.]]);

        assert_eq!(crossing, Some((0.5, 0.5)));
    }

    #[test]
    fn test_segments_which_do_not_reach_each_other_do_not_cross() {
        let crossing = segment_crossing([[0., 0.], [1., 0.]], [[2., -1.], [2., 1.]]);

        assert_eq!(crossing, None);
    }
}
