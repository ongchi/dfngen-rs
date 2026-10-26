//! Lines where the boundary of the region of a fracture crosses the fracture
//!
//! The plane of a fracture cuts each tetrahedron of the cells of its region.
//! Each cross-section is clipped with the fracture polygon, and the union of
//! these convex pieces is the part of the fracture in the region. The
//! boundary of the union, without its edges on the fracture outline, gives
//! the boundary lines. The fracture itself is not cut: a mesher meshes it with
//! these lines. With FRAM, the lines must obey the minimum feature size.

use std::collections::HashMap;

use i_overlay::core::fill_rule::FillRule;
use i_overlay::core::solver::Solver;
use i_overlay::float::overlay::OverlayOptions;
use i_overlay::float::simplify::SimplifyShape;
use parry3d_f64::bounding_volume::{Aabb, BoundingVolume};
use parry3d_f64::na::{distance as distance_3d, Point3, Vector3};
use parry3d_f64::query::PointQuery;
use parry3d_f64::shape::Segment;

use crate::domain_mesh::DomainMesh;
use crate::fram::RejectReason;
use crate::geometry::{are_segments_farther_than, line_seg_to_line_seg, PlaneFrame};

/// Distance below which a vertex of the union is on the outline of the
/// fracture, as a part of the size of the fracture. The union rounds its
/// coordinates to about 1e-9 of the size.
const RELATIVE_OUTLINE_TOLERANCE: f64 = 1e-7;

/// Largest absolute cosine of the angle between an intersection and a
/// boundary line which it crosses: cos(47°), as for a triple intersection
const MAX_CROSSING_COSINE: f64 = 0.681_998_36;

/// Part of `min_feature_size`: the point at `min_feature_size` along a line
/// from its end on the fracture edge must be at least this far from the edge,
/// as for an intersection which starts on an edge (`fram::check_close_edge`)
const END_DISTANCE_PART: f64 = 0.4;

/// Options of the boundary lines (input section `[region_boundaries]`). The
/// input has them only when the option is on.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RegionBoundaryOptions {
    /// True: the region of a stochastic fracture is the region of its family.
    /// False: the whole mesh.
    pub use_family_region: bool,
}

/// What is on the other side of a boundary line
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum OtherSide {
    /// A mesh cell with these regions (N is region N of `DomainMesh::regions`).
    /// Empty: the cell has no region.
    Cell(Vec<usize>),
    /// No cell: the line is on the surface of the mesh
    MeshSurface,
}

/// What is beside a boundary segment
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SegmentSide {
    pub other_side: OtherSide,
    /// The boundary surfaces of the mesh under the segment (N is boundary N
    /// of `DomainMesh::boundary_surfaces`), sorted. Empty: no boundary surface.
    pub surface_ids: Vec<usize>,
}

/// One straight segment of a boundary line
#[derive(Clone, Debug, PartialEq)]
pub struct BoundarySegment {
    pub start: Point3<f64>,
    pub end: Point3<f64>,
    pub side: SegmentSide,
}

/// Connected segments of the region boundary on a fracture, in order. An
/// open line starts and ends on the fracture outline; a closed line is a loop
/// in the fracture (around a hole of the region).
#[derive(Clone, Debug, PartialEq)]
pub struct BoundaryLine {
    pub segments: Vec<BoundarySegment>,
    pub is_closed: bool,
}

impl BoundaryLine {
    fn length(&self) -> f64 {
        self.segments
            .iter()
            .map(|segment| distance_3d(&segment.start, &segment.end))
            .sum()
    }

    /// The vertices of the line which are not on the fracture outline: all
    /// vertices of a closed line, the inner vertices of an open line
    fn inner_vertices(&self) -> impl Iterator<Item = &Point3<f64>> {
        let skipped = usize::from(!self.is_closed);
        self.segments
            .iter()
            .skip(skipped)
            .map(|segment| &segment.start)
    }

    /// The point at the distance along the line from its start, or from its
    /// end. `None` if the line is shorter.
    fn point_along(&self, length: f64, is_from_end: bool) -> Option<Point3<f64>> {
        let mut remaining = length;
        let mut ordered: Vec<(Point3<f64>, Point3<f64>)> = self
            .segments
            .iter()
            .map(|segment| (segment.start, segment.end))
            .collect();

        if is_from_end {
            ordered = ordered
                .into_iter()
                .rev()
                .map(|(start, end)| (end, start))
                .collect();
        }

        for (start, end) in ordered {
            let segment_length = distance_3d(&start, &end);

            if remaining <= segment_length {
                return Some(start + (end - start) * (remaining / segment_length));
            }

            remaining -= segment_length;
        }

        None
    }
}

/// Finds the boundary lines of a convex fracture: the lines where the
/// boundary of the cells of one or more regions crosses it, and the lines
/// where the mesh surface crosses it in cells of other regions. The faces
/// between two cells of these regions give no line. A segment on a surface
/// cell of the mesh gets the boundary surfaces of that cell.
///
/// A line on the mesh surface in cells of other regions can end on a line of
/// the regions, where the boundary of the regions meets the mesh surface.
///
/// # Arguments
///
/// * `vertices` - Vertices of the fracture, a convex polygon
/// * `normal` - Unit normal vector of the fracture
/// * `mesh` - Domain mesh
/// * `region_ids` - Regions of the fracture (N is region N of
///   `DomainMesh::regions`). No region: all cells of the mesh.
pub fn boundary_lines(
    vertices: &[Point3<f64>],
    normal: &Vector3<f64>,
    mesh: &DomainMesh,
    region_ids: &[usize],
) -> Vec<BoundaryLine> {
    let frame = PlaneFrame::new(vertices, normal);
    let outline: Vec<[f64; 2]> = vertices
        .iter()
        .map(|vertex| frame.to_plane(vertex))
        .collect();
    let counter_clockwise_outline = counter_clockwise(outline.clone());
    let bounds = Aabb::from_points(vertices).loosened(mesh.tolerance());
    let outline_tolerance = RELATIVE_OUTLINE_TOLERANCE * bounds.extents().norm();

    // The pieces in all cells, and true for each piece in the regions
    let mut pieces: Vec<Vec<[f64; 2]>> = Vec::new();
    let mut is_piece_in_regions: Vec<bool> = Vec::new();
    let mut exact_points = ExactPoints::new(outline_tolerance);

    for (corners, is_in_regions) in mesh.plane_sections(&vertices[0], normal, &bounds, region_ids) {
        let plane_corners: Vec<[f64; 2]> = corners
            .iter()
            .map(|corner| frame.to_plane(corner))
            .collect();
        let piece = clip_convex(
            &sort_around_center(plane_corners.clone()),
            &counter_clockwise_outline,
        );

        if piece.len() < 3 || signed_area(&piece) <= 0. {
            continue;
        }

        for &point in &piece {
            // A corner of the cross-section keeps its point in space. A point of
            // the clip is on the outline of the fracture.
            let corner = plane_corners
                .iter()
                .position(|&plane_corner| plane_corner == point)
                .map(|i| corners[i]);
            exact_points.insert(point, corner);
        }

        pieces.push(piece);
        is_piece_in_regions.push(is_in_regions);
    }

    let region_pieces: Vec<Vec<[f64; 2]>> = pieces
        .iter()
        .zip(&is_piece_in_regions)
        .filter(|(_, is_in_regions)| **is_in_regions)
        .map(|(piece, _)| piece.clone())
        .collect();
    let boundary = BoundarySide {
        frame: &frame,
        outline: &outline,
        outline_tolerance,
        exact_points: &exact_points,
        mesh,
    };

    let mut lines: Vec<BoundaryLine> = union(&region_pieces)
        .iter()
        .flat_map(|contour| {
            boundary.lines_of_contour(contour, vertices, |start, end| {
                boundary.segment_side(start, end)
            })
        })
        .collect();

    if region_pieces.len() < pieces.len() {
        // The mesh surface in the cells of other regions: the edges of the
        // union of all pieces on the mesh surface, with a cell of another
        // region on their inner side
        let is_in_regions = |cell_regions: &[usize]| {
            region_ids
                .iter()
                .any(|region_id| cell_regions.contains(region_id))
        };
        lines.extend(union(&pieces).iter().flat_map(|contour| {
            boundary.lines_of_contour(contour, vertices, |start, end| {
                let side = boundary.segment_side(start, end)?;
                let inner_cell = boundary
                    .mesh
                    .regions_at(&boundary.probe(start, end, Side::Inner));
                let is_surface_of_other_regions = side.other_side == OtherSide::MeshSurface
                    && inner_cell.is_some_and(|cell_regions| !is_in_regions(cell_regions));
                is_surface_of_other_regions.then_some(side)
            })
        }));
    }

    lines
}

/// Contours of the union of the pieces. Outer contours are counter-clockwise
/// and holes clockwise: the union is on the left side of each edge.
fn union(pieces: &[Vec<[f64; 2]>]) -> Vec<Vec<[f64; 2]>> {
    let mut options = OverlayOptions::default();
    // The points where the cell on the other side of an edge can change stay
    options.preserve_output_collinear = true;

    pieces
        .simplify_shape_custom(FillRule::NonZero, options, Solver::default())
        .into_iter()
        .flatten()
        .collect()
}

/// Side of an edge of a contour of the union
#[derive(Clone, Copy)]
enum Side {
    /// In the union
    Inner,
    /// Outside of the union
    Outer,
}

/// What is on the other side of each edge of the contours of the union
struct BoundarySide<'a> {
    frame: &'a PlaneFrame,
    /// Vertices of the fracture in the plane, in the order of the fracture
    outline: &'a [[f64; 2]],
    outline_tolerance: f64,
    exact_points: &'a ExactPoints,
    mesh: &'a DomainMesh,
}

impl BoundarySide<'_> {
    /// The boundary lines of one contour of the union: the runs of edges
    /// which `line_side` accepts. Collinear edges with the same side become
    /// one segment.
    ///
    /// # Arguments
    ///
    /// * `contour` - Contour of the union
    /// * `vertices` - Vertices of the fracture
    /// * `line_side` - Side of an edge of a line, from its start and end;
    ///   `None` for an edge which is not on a line
    fn lines_of_contour(
        &self,
        contour: &[[f64; 2]],
        vertices: &[Point3<f64>],
        line_side: impl Fn([f64; 2], [f64; 2]) -> Option<SegmentSide>,
    ) -> Vec<BoundaryLine> {
        let count = contour.len();
        let other_sides: Vec<Option<SegmentSide>> = (0..count)
            .map(|i| line_side(contour[i], contour[(i + 1) % count]))
            .collect();

        // Vertex i is kept if the edges before and after it differ in kind or direction
        let kept: Vec<usize> = (0..count)
            .filter(|&i| {
                let previous = (i + count - 1) % count;
                other_sides[previous] != other_sides[i]
                    || !self.is_collinear(contour[previous], contour[i], contour[(i + 1) % count])
            })
            .collect();
        let points: Vec<Point3<f64>> = kept
            .iter()
            .map(|&i| self.to_space(contour[i], vertices))
            .collect();
        let sides: Vec<Option<SegmentSide>> =
            kept.iter().map(|&i| other_sides[i].clone()).collect();
        let segment = |i: usize, side: SegmentSide| BoundarySegment {
            start: points[i],
            end: points[(i + 1) % points.len()],
            side,
        };

        if sides.iter().all(Option::is_some) {
            return vec![BoundaryLine {
                segments: sides
                    .into_iter()
                    .enumerate()
                    .map(|(i, side)| segment(i, side.expect("all edges are boundary edges")))
                    .collect(),
                is_closed: true,
            }];
        }

        // Each line starts after an edge which is not on a line
        let mut lines = Vec::new();
        let edge_count = sides.len();
        let first_outline_edge = sides
            .iter()
            .position(Option::is_none)
            .expect("an edge is not on a line");
        let mut current: Vec<BoundarySegment> = Vec::new();

        for offset in 1..=edge_count {
            let i = (first_outline_edge + offset) % edge_count;

            match &sides[i] {
                Some(side) => current.push(segment(i, side.clone())),
                None if !current.is_empty() => lines.push(BoundaryLine {
                    segments: std::mem::take(&mut current),
                    is_closed: false,
                }),
                None => {}
            }
        }

        lines
    }

    /// What is beside the edge from `start` to `end`: the cell on its other
    /// side, and the boundary surfaces at its middle. `None` for an edge on
    /// the fracture outline.
    fn segment_side(&self, start: [f64; 2], end: [f64; 2]) -> Option<SegmentSide> {
        let is_on_outline_edge = (0..self.outline.len()).any(|i| {
            let edge = [self.outline[i], self.outline[(i + 1) % self.outline.len()]];
            distance_to_segment(start, edge) < self.outline_tolerance
                && distance_to_segment(end, edge) < self.outline_tolerance
        });

        if is_on_outline_edge {
            return None;
        }

        let other_side = match self.mesh.regions_at(&self.probe(start, end, Side::Outer)) {
            Some(region_ids) => OtherSide::Cell(region_ids.to_vec()),
            None => OtherSide::MeshSurface,
        };
        let middle = self
            .frame
            .to_space([0.5 * (start[0] + end[0]), 0.5 * (start[1] + end[1])]);
        // The union rounds the coordinates to a part of the outline tolerance
        let surface_ids = self
            .mesh
            .boundary_surfaces_at(&middle, self.outline_tolerance.max(self.mesh.tolerance()));

        Some(SegmentSide {
            other_side,
            surface_ids,
        })
    }

    /// A point a small distance from the middle of the edge from `start` to
    /// `end`, on one side. The union is on the left side of the edge.
    fn probe(&self, start: [f64; 2], end: [f64; 2], side: Side) -> Point3<f64> {
        let direction = [end[0] - start[0], end[1] - start[1]];
        let length = direction[0].hypot(direction[1]);
        let tolerance = self.mesh.tolerance();
        let offset = (1e-3 * length).clamp(100. * tolerance, 1e5 * tolerance);
        let sign = match side {
            Side::Outer => 1.,
            Side::Inner => -1.,
        };

        self.frame.to_space([
            0.5 * (start[0] + end[0]) + sign * offset * direction[1] / length,
            0.5 * (start[1] + end[1]) - sign * offset * direction[0] / length,
        ])
    }

    fn is_collinear(&self, previous: [f64; 2], point: [f64; 2], next: [f64; 2]) -> bool {
        distance_to_segment(point, [previous, next]) < self.outline_tolerance
    }

    /// Point in space of a vertex of the union. The union rounds its
    /// coordinates, thus the vertex goes back to the vertex of a piece: a
    /// vertex of the fracture, a corner of a cross-section, or a point on the
    /// outline of the fracture.
    fn to_space(&self, point: [f64; 2], vertices: &[Point3<f64>]) -> Point3<f64> {
        let count = self.outline.len();

        if let Some(i) =
            (0..count).find(|&i| distance(point, self.outline[i]) < self.outline_tolerance)
        {
            return vertices[i];
        }

        let (point, corner) = self.exact_points.nearest(point).unwrap_or((point, None));

        if let Some(corner) = corner {
            return corner;
        }

        for i in 0..count {
            let edge = [self.outline[i], self.outline[(i + 1) % count]];

            if distance_to_segment(point, edge) < self.outline_tolerance {
                let fraction = position_on_segment(point, edge);
                return vertices[i] + (vertices[(i + 1) % count] - vertices[i]) * fraction;
            }
        }

        self.frame.to_space(point)
    }
}

/// FRAM check of the boundary lines of a fracture against the fracture.
///
/// * `ShortBoundaryLine` - A line is shorter than `min_feature_size`.
/// * `BoundaryCloseToEdge` - A vertex of a line which is not on the fracture
///   outline is nearer than `min_feature_size` to the outline; or, at an end
///   of an open line, the point at `min_feature_size` along the line is
///   nearer than `END_DISTANCE_PART × min_feature_size` to the outline (the
///   line meets the outline at a small angle).
///
/// # Arguments
///
/// * `lines` - Boundary lines of the fracture
/// * `vertices` - Vertices of the fracture
/// * `min_feature_size` - Minimum feature size
/// * `is_relaxed` - Relaxed FRAM: only the check of the short lines
pub fn check_boundary_lines(
    lines: &[BoundaryLine],
    vertices: &[Point3<f64>],
    min_feature_size: f64,
    is_relaxed: bool,
) -> Result<(), RejectReason> {
    let distance_to_outline = |point: &Point3<f64>| {
        (0..vertices.len())
            .map(|i| {
                Segment::new(vertices[i], vertices[(i + 1) % vertices.len()])
                    .distance_to_local_point(point, true)
            })
            .fold(f64::MAX, f64::min)
    };

    for line in lines {
        if line.length() < min_feature_size {
            return Err(RejectReason::ShortBoundaryLine);
        }

        if is_relaxed {
            continue;
        }

        if line
            .inner_vertices()
            .any(|vertex| distance_to_outline(vertex) < min_feature_size)
        {
            return Err(RejectReason::BoundaryCloseToEdge);
        }

        if !line.is_closed {
            for is_from_end in [false, true] {
                let point = line
                    .point_along(min_feature_size, is_from_end)
                    .expect("the line is not shorter than min_feature_size");

                if distance_to_outline(&point) < END_DISTANCE_PART * min_feature_size {
                    return Err(RejectReason::BoundaryCloseToEdge);
                }
            }
        }
    }

    Ok(())
}

/// True if an intersection on a fracture is too near to a boundary line of
/// the fracture (FRAM, `IntersectionCloseToBoundary`): nearer than
/// `min_feature_size` without a crossing; or it crosses the line nearer than
/// `min_feature_size` to an end of the intersection, or at an angle outside
/// of 47° to 133°.
///
/// # Arguments
///
/// * `lines` - Boundary lines of the fracture
/// * `intersection` - The two end points of the intersection
/// * `min_feature_size` - Minimum feature size
/// * `eps` - Epsilon value for floating point comparisons
pub fn is_intersection_near_boundary(
    lines: &[BoundaryLine],
    intersection: &[Point3<f64>; 2],
    min_feature_size: f64,
    eps: f64,
) -> bool {
    lines.iter().flat_map(|line| &line.segments).any(|segment| {
        let boundary = [segment.start, segment.end];

        if are_segments_farther_than(intersection, &boundary, min_feature_size) {
            return false;
        }

        let (distance, point) = line_seg_to_line_seg(intersection, &boundary);

        if distance > eps {
            return distance < min_feature_size;
        }

        let is_near_an_end = intersection
            .iter()
            .any(|end| distance_3d(end, &point) < min_feature_size);
        let cosine = (intersection[1] - intersection[0])
            .normalize()
            .dot(&(boundary[1] - boundary[0]).normalize());

        is_near_an_end || cosine.abs() > MAX_CROSSING_COSINE
    })
}

/// A vertex of a piece in the plane, and its exact point in space, if any
type ExactPoint = ([f64; 2], Option<Point3<f64>>);

/// The vertices of the pieces before the union, found by their position in
/// the plane. Each vertex can have its exact point in space.
struct ExactPoints {
    tolerance: f64,
    /// Points by the grid cell of size `tolerance` which has them
    cells: HashMap<(i64, i64), Vec<ExactPoint>>,
}

impl ExactPoints {
    fn new(tolerance: f64) -> Self {
        Self {
            tolerance,
            cells: HashMap::new(),
        }
    }

    fn cell_of(&self, [x, y]: [f64; 2]) -> (i64, i64) {
        (
            (x / self.tolerance).floor() as i64,
            (y / self.tolerance).floor() as i64,
        )
    }

    fn insert(&mut self, point: [f64; 2], space_point: Option<Point3<f64>>) {
        let cell = self.cell_of(point);
        self.cells
            .entry(cell)
            .or_default()
            .push((point, space_point));
    }

    /// The nearest vertex nearer than the tolerance, with its point in space, if any
    fn nearest(&self, point: [f64; 2]) -> Option<ExactPoint> {
        let (x, y) = self.cell_of(point);

        (x - 1..=x + 1)
            .flat_map(|cell_x| (y - 1..=y + 1).map(move |cell_y| (cell_x, cell_y)))
            .filter_map(|cell| self.cells.get(&cell))
            .flatten()
            .filter(|(exact, _)| distance(*exact, point) < self.tolerance)
            .min_by(|a, b| distance(a.0, point).total_cmp(&distance(b.0, point)))
            .copied()
    }
}

/// Twice the signed area of a polygon. Positive for counter-clockwise order.
fn signed_area(polygon: &[[f64; 2]]) -> f64 {
    (0..polygon.len())
        .map(|i| {
            let [x, y] = polygon[i];
            let [next_x, next_y] = polygon[(i + 1) % polygon.len()];
            x * next_y - next_x * y
        })
        .sum()
}

fn counter_clockwise(mut polygon: Vec<[f64; 2]>) -> Vec<[f64; 2]> {
    if signed_area(&polygon) < 0. {
        polygon.reverse();
    }
    polygon
}

/// The corners of a convex polygon in counter-clockwise order
fn sort_around_center(mut corners: Vec<[f64; 2]>) -> Vec<[f64; 2]> {
    let count = corners.len() as f64;
    let center_x = corners.iter().map(|corner| corner[0]).sum::<f64>() / count;
    let center_y = corners.iter().map(|corner| corner[1]).sum::<f64>() / count;
    let angle = |corner: &[f64; 2]| (corner[1] - center_y).atan2(corner[0] - center_x);

    corners.sort_by(|a, b| angle(a).total_cmp(&angle(b)));
    corners
}

/// Part of a convex polygon in a convex counter-clockwise polygon
/// (Sutherland-Hodgman)
fn clip_convex(polygon: &[[f64; 2]], clip: &[[f64; 2]]) -> Vec<[f64; 2]> {
    let mut result = polygon.to_vec();

    for i in 0..clip.len() {
        let [start, end] = [clip[i], clip[(i + 1) % clip.len()]];
        // Positive on the left side of the edge: inside
        let side = |point: [f64; 2]| {
            (end[0] - start[0]) * (point[1] - start[1])
                - (end[1] - start[1]) * (point[0] - start[0])
        };
        let input = std::mem::take(&mut result);

        for j in 0..input.len() {
            let current = input[j];
            let next = input[(j + 1) % input.len()];
            let (current_side, next_side) = (side(current), side(next));

            if current_side >= 0. {
                result.push(current);
            }

            if (current_side >= 0.) != (next_side >= 0.) {
                let fraction = current_side / (current_side - next_side);
                result.push([
                    current[0] + (next[0] - current[0]) * fraction,
                    current[1] + (next[1] - current[1]) * fraction,
                ]);
            }
        }
    }

    result
}

pub(crate) fn distance(a: [f64; 2], b: [f64; 2]) -> f64 {
    (a[0] - b[0]).hypot(a[1] - b[1])
}

/// Position of the nearest point of the segment, from 0 at its start to 1 at its end
pub(crate) fn position_on_segment(point: [f64; 2], [start, end]: [[f64; 2]; 2]) -> f64 {
    let direction = [end[0] - start[0], end[1] - start[1]];
    let length_squared = direction[0] * direction[0] + direction[1] * direction[1];

    if length_squared == 0. {
        return 0.;
    }

    (((point[0] - start[0]) * direction[0] + (point[1] - start[1]) * direction[1]) / length_squared)
        .clamp(0., 1.)
}

pub(crate) fn distance_to_segment(point: [f64; 2], segment: [[f64; 2]; 2]) -> f64 {
    let [start, end] = segment;
    let fraction = position_on_segment(point, segment);
    let nearest = [
        start[0] + (end[0] - start[0]) * fraction,
        start[1] + (end[1] - start[1]) * fraction,
    ];

    distance(point, nearest)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain_mesh::test_meshes::{box_tetrahedra, l_shaped_mesh, region_names};
    use crate::domain_mesh::MeshSurfaceTriangle;

    const MIN_FEATURE_SIZE: f64 = 0.1;
    const EPS: f64 = 1e-9;

    /// A square with half side `half` in the plane z = 0, centered at (x, y)
    fn square(x: f64, y: f64, half: f64) -> Vec<Point3<f64>> {
        [(-1., -1.), (1., -1.), (1., 1.), (-1., 1.)]
            .iter()
            .map(|(dx, dy)| Point3::new(x + dx * half, y + dy * half, 0.))
            .collect()
    }

    fn lines(
        vertices: &[Point3<f64>],
        mesh: &DomainMesh,
        region_ids: &[usize],
    ) -> Vec<BoundaryLine> {
        boundary_lines(vertices, &Vector3::z(), mesh, region_ids)
    }

    /// The other side of each segment of all lines, and the length of each line
    fn sides(lines: &[BoundaryLine]) -> Vec<&OtherSide> {
        lines
            .iter()
            .flat_map(|line| &line.segments)
            .map(|segment| &segment.side.other_side)
            .collect()
    }

    /// The boundary surfaces of each segment of all lines
    fn surface_ids(lines: &[BoundaryLine]) -> Vec<&[usize]> {
        lines
            .iter()
            .flat_map(|line| &line.segments)
            .map(|segment| segment.side.surface_ids.as_slice())
            .collect()
    }

    /// A box of 10 x 10 x 10 with its center at the origin. Its face at
    /// x = -5 has the surface triangles of the diagonal of its tetrahedra
    /// (from (-5, -5, -5) to (-5, 5, 5)): the triangle where z > y is on the
    /// boundary surface `west` (boundary 1) if `is_west_marked`, the
    /// triangle where z < y on no boundary surface.
    fn box_with_a_marked_triangle(is_west_marked: bool) -> DomainMesh {
        let tetrahedra = box_tetrahedra([-5.; 3], [5.; 3], &[]);
        let triangles = [
            MeshSurfaceTriangle {
                vertices: [
                    Point3::new(-5., -5., -5.),
                    Point3::new(-5., 5., 5.),
                    Point3::new(-5., -5., 5.),
                ],
                boundaries: if is_west_marked { vec![0] } else { Vec::new() },
            },
            MeshSurfaceTriangle {
                vertices: [
                    Point3::new(-5., -5., -5.),
                    Point3::new(-5., 5., -5.),
                    Point3::new(-5., 5., 5.),
                ],
                boundaries: Vec::new(),
            },
        ];

        DomainMesh::new(&tetrahedra, &[])
            .unwrap()
            .with_boundary_surfaces(&triangles, &region_names(&["west"]))
    }

    fn lengths(lines: &[BoundaryLine]) -> Vec<f64> {
        lines.iter().map(BoundaryLine::length).collect()
    }

    /// 3 x 3 boxes of 2 x 2 x 10 in x and y, without the center box
    fn mesh_with_a_hole() -> DomainMesh {
        let mut tetrahedra = Vec::new();
        for i in 0..3 {
            for j in 0..3 {
                if (i, j) != (1, 1) {
                    let (x, y) = (-3. + 2. * i as f64, -3. + 2. * j as f64);
                    tetrahedra.extend(box_tetrahedra([x, y, -5.], [x + 2., y + 2., 5.], &[]));
                }
            }
        }
        DomainMesh::new(&tetrahedra, &[]).unwrap()
    }

    /// An open line from (x0, y0) to (x1, y1) in the plane z = 0
    fn straight_line(x0: f64, y0: f64, x1: f64, y1: f64) -> BoundaryLine {
        BoundaryLine {
            segments: vec![BoundarySegment {
                start: Point3::new(x0, y0, 0.),
                end: Point3::new(x1, y1, 0.),
                side: SegmentSide {
                    other_side: OtherSide::MeshSurface,
                    surface_ids: Vec::new(),
                },
            }],
            is_closed: false,
        }
    }

    fn check(lines: &[BoundaryLine]) -> Result<(), RejectReason> {
        check_boundary_lines(lines, &square(0., 0., 1.), MIN_FEATURE_SIZE, false)
    }

    fn is_near(lines: &[BoundaryLine], from: [f64; 2], to: [f64; 2]) -> bool {
        let intersection = [
            Point3::new(from[0], from[1], 0.),
            Point3::new(to[0], to[1], 0.),
        ];
        is_intersection_near_boundary(lines, &intersection, MIN_FEATURE_SIZE, EPS)
    }

    #[test]
    fn test_fracture_in_the_cells_has_no_line() {
        assert!(lines(&square(-2.5, -2.5, 1.), &l_shaped_mesh(), &[]).is_empty());
    }

    #[test]
    fn test_fracture_outside_the_cells_has_no_line() {
        assert!(lines(&square(3., 3., 1.), &l_shaped_mesh(), &[]).is_empty());
    }

    #[test]
    fn test_fracture_across_the_corner_of_the_l_shape_has_one_line_on_the_mesh_surface() {
        // The mesh has no cell where x and y are positive: the line goes from
        // (2, 0) to (0, 0) to (0, 2)
        let lines = lines(&square(0., 0., 2.), &l_shaped_mesh(), &[]);

        assert_eq!(lengths(&lines), [4.]);
        assert_eq!(sides(&lines), [&OtherSide::MeshSurface; 2]);
    }

    #[test]
    fn test_line_of_one_region_has_the_other_region_beside_it() {
        // Region 1 (south) has the cells where y < 0; region 2 (north) is at y > 0, x < 0.
        // The line is at y = 0: north on the other side for x < 0, no cell for x > 0.
        let lines = lines(&square(0., 0., 2.), &l_shaped_mesh(), &[1]);

        let mut sides = sides(&lines[..1]);
        sides.sort_by_key(|side| matches!(side, OtherSide::MeshSurface));
        assert_eq!(sides, [&OtherSide::Cell(vec![2]), &OtherSide::MeshSurface]);
    }

    #[test]
    fn test_mesh_surface_in_cells_of_another_region_gives_a_line() {
        // The fracture of region 1 (south) crosses the mesh surface at x = 0,
        // y > 0 in the cells of region 2 (north). The line ends on the line of
        // region 1 at (0, 0).
        let lines = lines(&square(0., 0., 2.), &l_shaped_mesh(), &[1]);

        assert_eq!(lines.len(), 2);
        assert_eq!(sides(&lines[1..]), [&OtherSide::MeshSurface]);
        assert!((lines[1].length() - 2.).abs() < 1e-9);
    }

    #[test]
    fn test_fracture_outside_its_region_gets_its_mesh_surface_line() {
        // x from 3.5 to 5.5 and y < 0: only in cells of region 1 (south), and
        // across the mesh surface at x = 5. The region of the fracture is 2 (north).
        let lines = lines(&square(4.5, -2.5, 1.), &l_shaped_mesh(), &[2]);

        assert_eq!(sides(&lines), [&OtherSide::MeshSurface]);
    }

    #[test]
    fn test_two_regions_together_have_no_line_between_them() {
        let lines = lines(&square(0., 0., 2.), &l_shaped_mesh(), &[1, 2]);

        assert_eq!(sides(&lines), [&OtherSide::MeshSurface; 2]);
    }

    #[test]
    fn test_segment_on_a_marked_surface_cell_has_its_boundary_surface() {
        // The line at x = -5 from y = -1 to y = 0 is in the triangle where z > y
        let lines = lines(
            &square(-5., -0.5, 0.5),
            &box_with_a_marked_triangle(true),
            &[],
        );

        assert_eq!(surface_ids(&lines), [&[1][..]]);
    }

    #[test]
    fn test_segment_on_an_unmarked_surface_cell_has_no_boundary_surface() {
        let lines = lines(
            &square(-5., -0.5, 0.5),
            &box_with_a_marked_triangle(false),
            &[],
        );

        assert_eq!(surface_ids(&lines), [&[] as &[usize]]);
    }

    #[test]
    fn test_line_across_two_surface_cells_has_a_segment_for_each() {
        // The line at x = -5 from y = -1 to y = 1 crosses the diagonal at y = 0
        let lines = lines(&square(-5., 0., 1.), &box_with_a_marked_triangle(true), &[]);
        let mut ids = surface_ids(&lines);
        ids.sort();

        assert_eq!(ids, [&[] as &[usize], &[1]]);
    }

    #[test]
    fn test_hole_of_the_mesh_gives_a_closed_line() {
        let lines = lines(&square(0., 0., 2.), &mesh_with_a_hole(), &[]);

        assert_eq!(lines.len(), 1);
        assert!(lines[0].is_closed);
        assert!((lines[0].length() - 8.).abs() < 1e-9);
    }

    #[test]
    fn test_line_vertices_are_on_the_cell_faces() {
        let lines = lines(&square(0., 0., 2.), &l_shaped_mesh(), &[]);

        let corner = lines[0].segments[0].end;
        assert_eq!(corner, Point3::new(0., 0., 0.));
    }

    #[test]
    fn test_short_line_is_rejected() {
        let line = straight_line(1., 0.95, 0.95, 1.);

        assert_eq!(check(&[line]), Err(RejectReason::ShortBoundaryLine));
    }

    #[test]
    fn test_short_line_is_rejected_with_relaxed_fram() {
        let line = straight_line(1., 0.95, 0.95, 1.);

        let result = check_boundary_lines(&[line], &square(0., 0., 1.), MIN_FEATURE_SIZE, true);

        assert_eq!(result, Err(RejectReason::ShortBoundaryLine));
    }

    #[test]
    fn test_line_across_the_fracture_is_accepted() {
        assert_eq!(check(&[straight_line(0., -1., 0., 1.)]), Ok(()));
    }

    #[test]
    fn test_line_vertex_near_the_edge_is_rejected() {
        // From the edge y = -1 to (0.5, 0.95), 0.05 from the edge y = 1, and back
        let line = BoundaryLine {
            segments: vec![
                straight_line(0., -1., 0.5, 0.95).segments[0].clone(),
                straight_line(0.5, 0.95, 1., -0.5).segments[0].clone(),
            ],
            is_closed: false,
        };

        assert_eq!(check(&[line]), Err(RejectReason::BoundaryCloseToEdge));
    }

    #[test]
    fn test_line_at_a_small_angle_to_the_edge_is_rejected() {
        // From the edge y = -1 at less than 1° to it
        let line = straight_line(-0.5, -1., 1., -0.98);

        assert_eq!(check(&[line]), Err(RejectReason::BoundaryCloseToEdge));
    }

    #[test]
    fn test_line_at_a_large_angle_to_the_edge_is_accepted() {
        // From the edge x = -1 at 80° to it, to the edge x = 1
        let line = straight_line(-1., -0.8, 1., -0.8 + 2. * 10_f64.to_radians().tan());

        assert_eq!(check(&[line]), Ok(()));
    }

    #[test]
    fn test_intersection_near_a_line_is_rejected() {
        let line = straight_line(0., -1., 0., 1.);

        assert!(is_near(&[line], [0.05, -0.5], [0.05, 0.5]));
    }

    #[test]
    fn test_intersection_far_from_a_line_is_accepted() {
        let line = straight_line(0., -1., 0., 1.);

        assert!(!is_near(&[line], [0.5, -0.5], [0.5, 0.5]));
    }

    #[test]
    fn test_crossing_at_a_right_angle_far_from_the_ends_is_accepted() {
        let line = straight_line(0., -1., 0., 1.);

        assert!(!is_near(&[line], [-0.5, 0.], [0.5, 0.]));
    }

    #[test]
    fn test_crossing_near_an_end_of_the_intersection_is_rejected() {
        let line = straight_line(0., -1., 0., 1.);

        assert!(is_near(&[line], [-0.05, 0.], [0.5, 0.]));
    }

    #[test]
    fn test_crossing_at_a_small_angle_is_rejected() {
        let line = straight_line(0., -1., 0., 1.);
        let direction = [30_f64.to_radians().sin(), 30_f64.to_radians().cos()];

        assert!(is_near(
            &[line],
            [-0.5 * direction[0], -0.5 * direction[1]],
            [0.5 * direction[0], 0.5 * direction[1]]
        ));
    }
}
