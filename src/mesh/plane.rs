//! A fracture in its own plane: the outline, the lines which the mesh must
//! follow, and the target size of the triangles

use parry3d_f64::na::Point3;

use crate::fracture::poly::Poly;
use crate::geometry::PlaneFrame;
use crate::region_boundaries::{distance, distance_to_segment, OtherSide};

/// A point in the plane of a fracture
pub type Point2 = [f64; 2];

/// Corner angle (degrees) below which a thin fracture is a sliver
const SLIVER_ANGLE: f64 = 10.;

/// Part of `min_feature_size`: the width of a fracture limits the size to
/// no less than this (near a sharp corner the width goes to 0)
const MIN_WIDTH_SIZE_PART: f64 = 0.5;

/// Largest angle (degrees) between an outline edge and the reverse of
/// another edge for which the other edge is on the opposite side of the
/// fracture (see `FracturePlane::width_at`)
const OPPOSITE_ANGLE: f64 = 30.;

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
        let outline: Vec<Point2> = poly
            .vertices
            .iter()
            .map(|vertex| frame.to_plane(vertex))
            .collect();

        Self {
            outline,
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

    /// True if the fracture is a sliver: narrower than `min_feature_size`
    /// (2 × area / longest edge) and with a corner below `SLIVER_ANGLE`. Its
    /// triangles near that corner have small angles; no refinement can
    /// change that.
    pub fn is_sliver(&self) -> bool {
        let edges: Vec<[Point2; 2]> = self.outline_edges().collect();
        let twice_area: f64 = edges
            .iter()
            .map(|[start, end]| start[0] * end[1] - end[0] * start[1])
            .sum::<f64>()
            .abs();
        let longest = edges
            .iter()
            .map(|&[start, end]| distance(start, end))
            .fold(0., f64::max);
        let has_sharp_corner = (0..edges.len()).any(|i| {
            let [corner, next] = edges[i];
            let previous = edges[(i + edges.len() - 1) % edges.len()][0];
            let to_next = [next[0] - corner[0], next[1] - corner[1]];
            let to_previous = [previous[0] - corner[0], previous[1] - corner[1]];
            let cross = to_next[0] * to_previous[1] - to_next[1] * to_previous[0];
            let dot = to_next[0] * to_previous[0] + to_next[1] * to_previous[1];
            cross.abs().atan2(dot) < SLIVER_ANGLE.to_radians()
        });

        twice_area < self.min_size * longest && has_sharp_corner
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
    /// most the largest size, and at most the width of the fracture at the
    /// point (thus a thin fracture gets triangles of a good shape across it);
    /// the width limits the size to no less than `MIN_WIDTH_SIZE_PART` × h
    pub fn size_at(&self, point: Point2) -> f64 {
        let distance = self
            .feature_segments()
            .map(|segment| distance_to_segment(point, segment))
            .fold(f64::INFINITY, f64::min);

        (self.min_size + self.slope * distance)
            .min(self.max_size)
            .min(
                self.width_at(point)
                    .max(MIN_WIDTH_SIZE_PART * self.min_size),
            )
    }

    /// Width of the fracture at the point: the distance to the nearest
    /// outline edge which goes in the opposite direction of the outline edge
    /// nearest to the point (within `OPPOSITE_ANGLE`)
    fn width_at(&self, point: Point2) -> f64 {
        let direction = |[start, end]: [Point2; 2]| {
            let length = distance(start, end);
            [(end[0] - start[0]) / length, (end[1] - start[1]) / length]
        };
        let Some(nearest) = self.outline_edges().min_by(|a, b| {
            distance_to_segment(point, *a).total_cmp(&distance_to_segment(point, *b))
        }) else {
            return f64::INFINITY;
        };
        let nearest_direction = direction(nearest);
        let max_cosine = -OPPOSITE_ANGLE.to_radians().cos();

        self.outline_edges()
            .filter(|&edge| {
                let edge_direction = direction(edge);
                nearest_direction[0] * edge_direction[0] + nearest_direction[1] * edge_direction[1]
                    < max_cosine
            })
            .map(|edge| distance_to_segment(point, edge))
            .fold(f64::INFINITY, f64::min)
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
