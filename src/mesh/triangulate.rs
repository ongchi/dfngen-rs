//! Triangle mesh of one fracture in its plane

use std::collections::{HashMap, HashSet};

use spade::handles::{
    DirectedEdgeHandle, FaceHandle, FixedFaceHandle, FixedVertexHandle, InnerTag,
};
use spade::{
    mitigate_underflow, CdtEdge, ConstrainedDelaunayTriangulation, HasPosition, InsertionError,
    Point2 as SpadePoint, PositionInTriangulation, Triangulation,
};

use super::plane::{point_on_segment, FracturePlane, Point2};
use super::shared_points::{Join, SharedNetwork};
use super::{JOIN_DISTANCE_PART, MERGE_DISTANCE_PART};
use crate::region_boundaries::{distance, distance_to_segment, position_on_segment, OtherSide};

/// A triangle gets a new point if its circumradius is above this part of the
/// target size. An equilateral triangle with sides of the target size has
/// 1/√3 ≈ 0.58.
const REFINE_RADIUS_PART: f64 = 0.65;

/// Smallest angle of a triangle (degrees) which the refinement accepts,
/// except for an angle of the input (between two constraint edges)
const MIN_ANGLE: f64 = 20.;

/// Part of `min_feature_size`: a new point is at least this far from the
/// other points, thus the refinement ends at small features of the input
const MIN_SPACING_PART: f64 = 0.7;

/// Largest angle (degrees) between a joinable shared line and an outline
/// edge for which the line is along the edge (see `NearLines::of_edge`)
const ALONG_ANGLE: f64 = 30.;

/// Limit of the refinement passes; a normal mesh needs far fewer
const MAX_REFINEMENT_PASSES: usize = 200;

/// A point of the mesh of one fracture
#[derive(Clone, Copy, Debug)]
pub struct LocalVertex {
    pub position: Point2,
    /// Id in `SharedPoints`, for a point which other fractures use too
    pub shared_id: Option<usize>,
}

/// The mesh of one fracture
pub struct LocalMesh {
    pub vertices: Vec<LocalVertex>,
    /// Indices into `vertices`, counter-clockwise around the fracture normal
    pub triangles: Vec<[usize; 3]>,
    /// The vertices on the boundary lines, with what is on the other side
    pub boundary_sides: Vec<(usize, OtherSide)>,
    /// Number of constraint edges which crossed another one and were not added
    pub skipped_constraints: usize,
    /// Edges of shared lines (pairs of shared ids) which the refinement
    /// needs split; the network splits them in all fractures of the line
    pub shared_splits: Vec<[usize; 2]>,
}

/// The vertices of the mesh, found by their position in the plane
struct LocalVertices {
    vertices: Vec<LocalVertex>,
    cell_size: f64,
    cells: HashMap<(i64, i64), Vec<usize>>,
}

impl LocalVertices {
    fn new(cell_size: f64) -> Self {
        Self {
            vertices: Vec::new(),
            cell_size,
            cells: HashMap::new(),
        }
    }

    fn cell_of(&self, [x, y]: Point2) -> (i64, i64) {
        (
            (x / self.cell_size).floor() as i64,
            (y / self.cell_size).floor() as i64,
        )
    }

    fn add(&mut self, vertex: LocalVertex) -> usize {
        let index = self.vertices.len();
        self.cells
            .entry(self.cell_of(vertex.position))
            .or_default()
            .push(index);
        self.vertices.push(vertex);
        index
    }

    /// The nearest vertex nearer than the distance
    fn nearest_within(&self, point: Point2, max_distance: f64) -> Option<usize> {
        let (x, y) = self.cell_of(point);
        let reach = (max_distance / self.cell_size).ceil() as i64;

        (x - reach..=x + reach)
            .flat_map(|cell_x| (y - reach..=y + reach).map(move |cell_y| (cell_x, cell_y)))
            .filter_map(|cell| self.cells.get(&cell))
            .flatten()
            .copied()
            .filter(|&index| distance(self.vertices[index].position, point) < max_distance)
            .min_by(|&a, &b| {
                distance(self.vertices[a].position, point)
                    .total_cmp(&distance(self.vertices[b].position, point))
            })
    }

    /// The vertex nearer than the merge distance, or a new vertex
    fn find_or_add(&mut self, point: Point2, merge_distance: f64) -> usize {
        self.nearest_within(point, merge_distance)
            .unwrap_or_else(|| {
                self.add(LocalVertex {
                    position: point,
                    shared_id: None,
                })
            })
    }
}

#[derive(Clone, Copy)]
struct CdtVertex {
    position: SpadePoint<f64>,
    index: usize,
}

impl HasPosition for CdtVertex {
    type Scalar = f64;

    fn position(&self) -> SpadePoint<f64> {
        self.position
    }
}

type Cdt = ConstrainedDelaunayTriangulation<CdtVertex>;

/// The point for spade. A coordinate which is not 0 but nearer to 0 than
/// spade accepts becomes 0.
fn spade_point([x, y]: Point2) -> SpadePoint<f64> {
    mitigate_underflow(SpadePoint::new(x, y))
}

/// Meshes one fracture.
///
/// # Arguments
///
/// * `plane` - The fracture in its plane
/// * `fracture` - Position of the fracture in the final fractures
/// * `shared` - The shared points and lines of the network
pub fn mesh_fracture(
    plane: &FracturePlane,
    fracture: usize,
    shared: &SharedNetwork,
) -> Result<LocalMesh, InsertionError> {
    let merge_distance = MERGE_DISTANCE_PART * plane.min_size;
    let mut vertices = LocalVertices::new(plane.min_size);

    // The shared points first: the other points merge with them. Two shared
    // points can be at one position in this plane, but apart in space (a
    // user-defined polygon is not exactly planar): they are one vertex here.
    let mut local_of_shared: HashMap<usize, usize> = HashMap::new();
    for &id in &shared.points_of_fracture[fracture] {
        let position = plane.frame.to_plane(&shared.shared_points.points[id]);
        let index = vertices
            .nearest_within(position, merge_distance)
            .unwrap_or_else(|| {
                vertices.add(LocalVertex {
                    position,
                    shared_id: Some(id),
                })
            });
        local_of_shared.insert(id, index);
    }
    let mut shared_vertices: Vec<usize> = local_of_shared.values().copied().collect();
    shared_vertices.sort_unstable();
    shared_vertices.dedup();

    let mut polylines: Vec<Vec<usize>> = shared
        .lines
        .iter()
        .filter(|line| line.fractures.contains(&fracture))
        .map(|line| {
            line.point_ids
                .iter()
                .map(|id| local_of_shared[id])
                .collect()
        })
        .collect();

    let shared_line_count = polylines.len();
    let joinable_polylines: Vec<&Vec<usize>> = shared
        .lines
        .iter()
        .filter(|line| line.fractures.contains(&fracture))
        .zip(&polylines)
        .filter(|(line, _)| line.join == Join::Joinable)
        .map(|(_, polyline)| polyline)
        .collect();

    let boundary_ends: Vec<Point2> = plane
        .boundary_segments
        .iter()
        .flat_map(|boundary| boundary.segment)
        .collect();

    // Shared points on no shared line: corners of pieces, and points where
    // two fractures only touch
    let line_vertices: HashSet<usize> = polylines.iter().flatten().copied().collect();
    let contact_vertices: Vec<usize> = shared_vertices
        .iter()
        .copied()
        .filter(|index| !line_vertices.contains(index))
        .collect();

    let mut outline_polylines = Vec::new();
    for edge in plane.outline_edges() {
        if is_shared_edge(edge, &polylines, &vertices, merge_distance) {
            continue;
        }
        let near = NearLines::of_edge(
            edge,
            &joinable_polylines,
            &contact_vertices,
            &vertices,
            plane,
        );
        outline_polylines.push(divide_segment(
            edge,
            &boundary_ends,
            &shared_vertices,
            &mut vertices,
            plane,
            &near,
        ));
    }
    polylines.extend(outline_polylines);

    let mut edges = ConstraintEdges {
        shared: polylines[..shared_line_count]
            .iter()
            .flat_map(|polyline| polyline.windows(2))
            .map(|pair| edge_key(pair[0], pair[1]))
            .collect(),
        boundary: HashMap::new(),
    };

    let mut boundary_sides = Vec::new();
    for boundary in &plane.boundary_segments {
        let polyline = divide_segment(
            boundary.segment,
            &[],
            &shared_vertices,
            &mut vertices,
            plane,
            &NearLines::default(),
        );
        boundary_sides.extend(
            polyline
                .iter()
                .map(|&index| (index, boundary.other_side.clone())),
        );
        for pair in polyline.windows(2) {
            edges
                .boundary
                .insert(edge_key(pair[0], pair[1]), boundary.other_side.clone());
        }
        polylines.push(polyline);
    }

    let mut cdt = Cdt::new();
    let mut handles = Vec::with_capacity(vertices.vertices.len());
    for (index, vertex) in vertices.vertices.iter().enumerate() {
        handles.push(cdt.insert(CdtVertex {
            position: spade_point(vertex.position),
            index,
        })?);
    }

    let mut skipped_constraints = 0;
    for polyline in &polylines {
        for pair in polyline.windows(2) {
            let (from, to) = (handles[pair[0]], handles[pair[1]]);
            if from != to
                && !cdt.exists_constraint(from, to)
                && cdt.try_add_constraint(from, to).is_empty()
            {
                skipped_constraints += 1;
            }
        }
    }

    let shared_splits = refine(
        &mut cdt,
        &mut vertices,
        plane,
        &mut edges,
        &mut boundary_sides,
        &mut skipped_constraints,
    )?;
    let shared_splits = shared_splits
        .into_iter()
        .filter_map(|(first, second)| {
            Some([
                vertices.vertices[first].shared_id?,
                vertices.vertices[second].shared_id?,
            ])
        })
        .collect();

    // A skipped constraint edge can be a gap in the outline, and a path from
    // the hull can then go into the fracture. Then the centroid test.
    let outside = outside_faces(&cdt);
    let triangles: Vec<[usize; 3]> = cdt
        .inner_faces()
        .filter(|face| {
            if skipped_constraints == 0 {
                !outside.contains(&face.fix())
            } else {
                is_centroid_inside(face.positions(), plane, merge_distance)
            }
        })
        .map(|face| face.vertices().map(|vertex| vertex.data().index))
        .collect();

    Ok(LocalMesh {
        vertices: vertices.vertices,
        triangles,
        boundary_sides,
        skipped_constraints,
        shared_splits,
    })
}

/// The triangles outside the outline: those which a path from the convex
/// hull reaches without a crossing of a constraint edge. A test of the
/// centroid is not sufficient: points on a straight outline edge are not
/// exactly on one line, and the triangles between them and the hull have
/// almost no area.
fn outside_faces(cdt: &Cdt) -> HashSet<FixedFaceHandle<InnerTag>> {
    let mut outside = HashSet::new();
    let mut stack: Vec<FixedFaceHandle<InnerTag>> = cdt
        .convex_hull()
        .filter(|edge| !edge.is_constraint_edge())
        .filter_map(|edge| edge.rev().face().as_inner())
        .map(|face| face.fix())
        .collect();

    while let Some(face) = stack.pop() {
        if !outside.insert(face) {
            continue;
        }

        stack.extend(
            cdt.face(face)
                .adjacent_edges()
                .into_iter()
                .filter(|edge| !edge.is_constraint_edge())
                .filter_map(|edge| edge.rev().face().as_inner())
                .map(|next| next.fix())
                .filter(|next| !outside.contains(next)),
        );
    }

    outside
}

/// True if the centroid of the triangle is in the outline, and not on it
/// (the triangles between the outline and the hull have their centroid on it)
fn is_centroid_inside(
    [a, b, c]: [SpadePoint<f64>; 3],
    plane: &FracturePlane,
    merge_distance: f64,
) -> bool {
    let centroid = [(a.x + b.x + c.x) / 3., (a.y + b.y + c.y) / 3.];

    plane.contains(centroid)
        && plane
            .outline_edges()
            .all(|edge| distance_to_segment(centroid, edge) > merge_distance)
}

/// True if a shared line (a common edge of two pieces) is on the edge
fn is_shared_edge(
    [start, end]: [Point2; 2],
    polylines: &[Vec<usize>],
    vertices: &LocalVertices,
    merge_distance: f64,
) -> bool {
    let is_near = |index: usize, point: Point2| {
        distance(vertices.vertices[index].position, point) < merge_distance
    };

    polylines.iter().any(|polyline| {
        let (first, last) = (polyline[0], polyline[polyline.len() - 1]);
        (is_near(first, start) && is_near(last, end))
            || (is_near(first, end) && is_near(last, start))
    })
}

/// Joinable shared lines near an outline edge (tiny features of
/// user-defined fractures which FRAM did not check)
#[derive(Default)]
struct NearLines {
    /// Vertices which the edge must go through
    vertices: Vec<usize>,
    /// Pairs of these vertices which are consecutive on a shared line: the
    /// edge follows the line there and gets no other points
    covered: HashSet<(usize, usize)>,
}

impl NearLines {
    /// A line with both ends within the join distance of the edge, longer
    /// than that distance and within `ALONG_ANGLE` of the edge direction is
    /// along the edge: the edge uses all its points. Of another line, only
    /// an end within the join distance is on the edge, as is a shared point
    /// on no line (a point contact) within the join distance.
    fn of_edge(
        edge: [Point2; 2],
        joinable_polylines: &[&Vec<usize>],
        contact_vertices: &[usize],
        vertices: &LocalVertices,
        plane: &FracturePlane,
    ) -> Self {
        let join_distance = JOIN_DISTANCE_PART * plane.min_size;
        let is_near = |index: usize| {
            distance_to_segment(vertices.vertices[index].position, edge) < join_distance
        };
        let is_along = |first: usize, last: usize| {
            let [start, end] = [first, last].map(|index| vertices.vertices[index].position);
            let line = [end[0] - start[0], end[1] - start[1]];
            let edge_direction = [edge[1][0] - edge[0][0], edge[1][1] - edge[0][1]];
            let length = distance(start, end);
            let cosine = (line[0] * edge_direction[0] + line[1] * edge_direction[1])
                / (length * distance(edge[0], edge[1]));

            is_near(first)
                && is_near(last)
                && length > join_distance
                && cosine.abs() > ALONG_ANGLE.to_radians().cos()
        };
        let mut near = Self::default();

        for polyline in joinable_polylines {
            let (first, last) = (polyline[0], polyline[polyline.len() - 1]);
            if is_along(first, last) {
                near.vertices.extend(polyline.iter());
                near.covered
                    .extend(polyline.windows(2).map(|pair| edge_key(pair[0], pair[1])));
            } else {
                near.vertices
                    .extend([first, last].into_iter().filter(|&end| is_near(end)));
            }
        }

        near.vertices.extend(
            contact_vertices
                .iter()
                .copied()
                .filter(|&index| is_near(index)),
        );

        near
    }
}

/// The vertices along a segment: its ends, the points on it (shared
/// vertices and the extra points), the vertices of the near lines, and points
/// between them so that no part is longer than the target size at its middle
fn divide_segment(
    segment: [Point2; 2],
    extra_points: &[Point2],
    shared_vertices: &[usize],
    vertices: &mut LocalVertices,
    plane: &FracturePlane,
    near: &NearLines,
) -> Vec<usize> {
    let merge_distance = MERGE_DISTANCE_PART * plane.min_size;

    let mut on_segment: Vec<(f64, usize)> = shared_vertices
        .iter()
        .map(|&index| (vertices.vertices[index].position, index))
        .filter(|(position, _)| distance_to_segment(*position, segment) < merge_distance)
        .map(|(position, index)| (position_on_segment(position, segment), index))
        .collect();
    on_segment.extend(near.vertices.iter().map(|&index| {
        (
            position_on_segment(vertices.vertices[index].position, segment),
            index,
        )
    }));

    for point in segment.iter().chain(extra_points) {
        if distance_to_segment(*point, segment) < merge_distance {
            let index = vertices.find_or_add(*point, merge_distance);
            on_segment.push((position_on_segment(*point, segment), index));
        }
    }

    on_segment.sort_by(|a, b| a.0.total_cmp(&b.0));
    on_segment.dedup_by_key(|(_, index)| *index);

    let mut polyline = vec![on_segment[0].1];
    for pair in on_segment.windows(2) {
        if near.covered.contains(&edge_key(pair[0].1, pair[1].1)) {
            polyline.push(pair[1].1);
            continue;
        }
        let start = vertices.vertices[pair[0].1].position;
        let end = vertices.vertices[pair[1].1].position;
        for point in divide_2d(start, end, &|point| plane.size_at(point)) {
            polyline.push(vertices.add(LocalVertex {
                position: point,
                shared_id: None,
            }));
        }
        polyline.push(pair[1].1);
    }

    polyline
}

/// Points between the two points (not the two points), so that no part is
/// longer than the size at its middle (by halving)
fn divide_2d(start: Point2, end: Point2, size: &impl Fn(Point2) -> f64) -> Vec<Point2> {
    let middle = point_on_segment([start, end], 0.5);

    if distance(start, end) <= size(middle) {
        return Vec::new();
    }

    let mut points = divide_2d(start, middle, size);
    points.push(middle);
    points.extend(divide_2d(middle, end, size));
    points
}

/// The constraint edges of a fracture which the refinement must know
struct ConstraintEdges {
    /// Edges of the shared lines (local vertex indices, smaller first): only
    /// a split in all fractures of the line keeps the mesh connected
    shared: HashSet<(usize, usize)>,
    /// Edges of the boundary lines, with what is on their other side
    boundary: HashMap<(usize, usize), OtherSide>,
}

fn edge_key(first: usize, second: usize) -> (usize, usize) {
    (first.min(second), first.max(second))
}

/// What the refinement found to do for one bad triangle
enum Refinement {
    /// Insert the circumcenter
    Insert(Point2),
    /// Split a constraint edge at its middle (the circumcenter encroaches it)
    Split([FixedVertexHandle; 2]),
}

/// Adds points until no triangle in the outline is larger than the target
/// size or has an angle below `MIN_ANGLE` (Ruppert type): the circumcenter
/// of a bad triangle, or, if it encroaches a constraint edge (it is in the
/// diametral circle of the edge, or the edge is between it and the
/// triangle), the middle of that edge. An edge of a shared line is not split
/// here: its split is returned, and the network splits it in all fractures.
///
/// Exceptions, thus the refinement ends: an angle between two constraint
/// edges is an angle of the input and is kept; no new point is nearer than
/// `MIN_SPACING_PART` × h to another point; no edge shorter than two times
/// that is split.
///
/// Returns the shared edges to split (local vertex indices).
fn refine(
    cdt: &mut Cdt,
    vertices: &mut LocalVertices,
    plane: &FracturePlane,
    edges: &mut ConstraintEdges,
    boundary_sides: &mut Vec<(usize, OtherSide)>,
    skipped_constraints: &mut usize,
) -> Result<HashSet<(usize, usize)>, InsertionError> {
    let min_spacing = MIN_SPACING_PART * plane.min_size;
    let merge_distance = MERGE_DISTANCE_PART * plane.min_size;
    let mut shared_splits = HashSet::new();

    for _ in 0..MAX_REFINEMENT_PASSES {
        let outside = outside_faces(cdt);
        let refinements: Vec<Refinement> = cdt
            .inner_faces()
            .filter(|face| {
                if *skipped_constraints == 0 {
                    !outside.contains(&face.fix())
                } else {
                    is_centroid_inside(face.positions(), plane, merge_distance)
                }
            })
            .filter(|face| is_bad(face, plane))
            .filter_map(|face| {
                let [a, b, c] = face.positions();
                let centroid = SpadePoint::new((a.x + b.x + c.x) / 3., (a.y + b.y + c.y) / 3.);
                let circumcenter = face.circumcenter();
                // A face with its points on one line has no circumcenter; spade
                // cannot locate such a point
                if !(circumcenter.x.is_finite() && circumcenter.y.is_finite()) {
                    return None;
                }
                Some(match encroached_edge(cdt, centroid, circumcenter) {
                    Some(edge) => Refinement::Split(edge),
                    None => Refinement::Insert([circumcenter.x, circumcenter.y]),
                })
            })
            .collect();

        let mut changes = 0;
        for refinement in refinements {
            match refinement {
                Refinement::Insert(point) => {
                    if !plane.contains(point)
                        || vertices.nearest_within(point, min_spacing).is_some()
                    {
                        continue;
                    }
                    let index = vertices.add(LocalVertex {
                        position: point,
                        shared_id: None,
                    });
                    cdt.insert(CdtVertex {
                        position: spade_point(point),
                        index,
                    })?;
                    changes += 1;
                }
                Refinement::Split([from, to]) => {
                    // The edge can be gone after an earlier change of this pass
                    if !cdt.exists_constraint(from, to) {
                        continue;
                    }
                    let (first, second) =
                        (cdt.vertex(from).data().index, cdt.vertex(to).data().index);
                    let length = distance(
                        vertices.vertices[first].position,
                        vertices.vertices[second].position,
                    );
                    if length < 2. * min_spacing {
                        continue;
                    }
                    if edges.shared.contains(&edge_key(first, second)) {
                        shared_splits.insert(edge_key(first, second));
                        continue;
                    }
                    split_edge(
                        cdt,
                        vertices,
                        [from, to],
                        edges,
                        boundary_sides,
                        skipped_constraints,
                    )?;
                    changes += 1;
                }
            }
        }

        if changes == 0 {
            break;
        }
    }

    Ok(shared_splits)
}

/// True if the triangle is larger than the target size, or has an angle
/// below `MIN_ANGLE` which is not between two constraint edges
fn is_bad(
    face: &FaceHandle<InnerTag, CdtVertex, (), CdtEdge<()>, ()>,
    plane: &FracturePlane,
) -> bool {
    let positions = face.positions().map(|point| [point.x, point.y]);
    let [a, b, c] = positions;
    let centroid = [(a[0] + b[0] + c[0]) / 3., (a[1] + b[1] + c[1]) / 3.];
    let circumcenter = face.circumcenter();
    let radius = distance([circumcenter.x, circumcenter.y], a);

    if radius > REFINE_RADIUS_PART * plane.size_at(centroid) {
        return true;
    }

    // Edge i goes from vertex i to vertex i + 1; the angle at vertex i is
    // between edge i and edge i + 2
    let edges = face.adjacent_edges();
    (0..3).any(|i| {
        let vertex = positions[i];
        let next = positions[(i + 1) % 3];
        let previous = positions[(i + 2) % 3];
        let angle = angle_between(
            [next[0] - vertex[0], next[1] - vertex[1]],
            [previous[0] - vertex[0], previous[1] - vertex[1]],
        );
        let is_input_angle =
            edges[i].is_constraint_edge() && edges[(i + 2) % 3].is_constraint_edge();
        angle < MIN_ANGLE.to_radians() && !is_input_angle
    })
}

fn angle_between(first: Point2, second: Point2) -> f64 {
    let cross = first[0] * second[1] - first[1] * second[0];
    let dot = first[0] * second[0] + first[1] * second[1];
    cross.abs().atan2(dot)
}

/// The constraint edge which the point encroaches: a constraint edge between
/// the triangle (its centroid) and the point, or a constraint edge near the
/// point which has the point in its diametral circle
fn encroached_edge(
    cdt: &Cdt,
    centroid: SpadePoint<f64>,
    point: SpadePoint<f64>,
) -> Option<[FixedVertexHandle; 2]> {
    let vertices_of = |edge: DirectedEdgeHandle<CdtVertex, (), CdtEdge<()>, ()>| {
        let [from, to] = edge.vertices();
        [from.fix(), to.fix()]
    };

    if let Some(edge) = cdt
        .get_conflicting_edges_between_points(centroid, point)
        .next()
    {
        return Some(vertices_of(edge));
    }

    let is_in_diametral_circle = |edge: &DirectedEdgeHandle<CdtVertex, (), CdtEdge<()>, ()>| {
        let [from, to] = edge.positions();
        let middle = [(from.x + to.x) / 2., (from.y + to.y) / 2.];
        distance([point.x, point.y], middle) < distance([from.x, from.y], [to.x, to.y]) / 2.
    };

    let faces: Vec<FixedFaceHandle<InnerTag>> = match cdt.locate(point) {
        PositionInTriangulation::OnFace(face) => vec![face],
        PositionInTriangulation::OnEdge(edge) => {
            let edge = cdt.directed_edge(edge);
            if edge.is_constraint_edge() {
                return Some(vertices_of(edge));
            }
            [edge.face(), edge.rev().face()]
                .into_iter()
                .filter_map(|face| face.as_inner())
                .map(|face| face.fix())
                .collect()
        }
        _ => return None,
    };

    // The faces at the point and their neighbors
    let mut near_faces = faces.clone();
    for &face in &faces {
        near_faces.extend(
            cdt.face(face)
                .adjacent_edges()
                .into_iter()
                .filter_map(|edge| edge.rev().face().as_inner())
                .map(|neighbor| neighbor.fix()),
        );
    }

    near_faces
        .into_iter()
        .flat_map(|face| cdt.face(face).adjacent_edges())
        .find(|edge| edge.is_constraint_edge() && is_in_diametral_circle(edge))
        .map(vertices_of)
}

/// Splits a constraint edge of this fracture at its middle. The vertex on a
/// boundary line gets the other side of the line.
fn split_edge(
    cdt: &mut Cdt,
    vertices: &mut LocalVertices,
    [from, to]: [FixedVertexHandle; 2],
    edges: &mut ConstraintEdges,
    boundary_sides: &mut Vec<(usize, OtherSide)>,
    skipped_constraints: &mut usize,
) -> Result<(), InsertionError> {
    let (first, second) = (cdt.vertex(from).data().index, cdt.vertex(to).data().index);
    let middle = point_on_segment(
        [
            vertices.vertices[first].position,
            vertices.vertices[second].position,
        ],
        0.5,
    );

    if let Some(edge) = cdt.get_edge_from_neighbors(from, to) {
        let edge = edge.fix().as_undirected();
        cdt.remove_constraint_edge(edge);
    }

    let index = vertices.add(LocalVertex {
        position: middle,
        shared_id: None,
    });
    let middle_handle = cdt.insert(CdtVertex {
        position: spade_point(middle),
        index,
    })?;

    for (start, end) in [(from, middle_handle), (middle_handle, to)] {
        if start != end
            && !cdt.exists_constraint(start, end)
            && cdt.try_add_constraint(start, end).is_empty()
        {
            *skipped_constraints += 1;
        }
    }

    if let Some(other_side) = edges.boundary.remove(&edge_key(first, second)) {
        edges
            .boundary
            .insert(edge_key(first, index), other_side.clone());
        edges
            .boundary
            .insert(edge_key(index, second), other_side.clone());
        boundary_sides.push((index, other_side));
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use parry3d_f64::na::{Point3, Vector3};

    use super::*;
    use crate::fracture::poly::Poly;
    use crate::mesh::shared_points::{build_shared_network, LineOnFractures};

    /// Sizes: `min_feature_size`, largest size, slope
    const SIZES: [f64; 3] = [0.5, 2., 0.5];

    /// Two squares of 10 × 10 with the center at the origin: one in the plane
    /// z = 0 and one in the plane y = 0. They intersect along the x axis.
    fn crossing_squares() -> (Vec<FracturePlane>, SharedNetwork) {
        let square = |vertices: [[f64; 3]; 4], normal: Vector3<f64>| Poly {
            vertices: vertices.map(Point3::from).to_vec(),
            normal,
            ..Default::default()
        };
        let horizontal = square(
            [[-5., -5., 0.], [5., -5., 0.], [5., 5., 0.], [-5., 5., 0.]],
            Vector3::z(),
        );
        let vertical = square(
            [[-5., 0., -5.], [5., 0., -5.], [5., 0., 5.], [-5., 0., 5.]],
            -Vector3::y(),
        );
        let intersection = [Point3::new(-5., 0., 0.), Point3::new(5., 0., 0.)];

        let planes = vec![
            FracturePlane::new(&horizontal, &[intersection], SIZES),
            FracturePlane::new(&vertical, &[intersection], SIZES),
        ];
        let lines = [LineOnFractures {
            fractures: [0, 1],
            ends: intersection,
            join: Join::None,
        }];
        let network = build_shared_network(&planes, &[0, 1], &lines, &[], &[]);

        (planes, network)
    }

    fn triangle_area(mesh: &LocalMesh, triangle: &[usize; 3]) -> f64 {
        let [a, b, c] = triangle.map(|index| mesh.vertices[index].position);
        ((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])) / 2.
    }

    #[test]
    fn test_triangles_cover_the_square_counter_clockwise() {
        let (planes, network) = crossing_squares();

        let mesh = mesh_fracture(&planes[0], 0, &network).unwrap();

        let areas: Vec<f64> = mesh
            .triangles
            .iter()
            .map(|triangle| triangle_area(&mesh, triangle))
            .collect();
        assert!(areas.iter().all(|&area| area > 0.));
        assert!((areas.iter().sum::<f64>() - 100.).abs() < 1e-9);
    }

    #[test]
    fn test_both_fractures_use_each_point_of_the_intersection() {
        let (planes, network) = crossing_squares();

        for (fracture, plane) in planes.iter().enumerate() {
            let mesh = mesh_fracture(plane, fracture, &network).unwrap();
            let used_ids: Vec<usize> = mesh
                .triangles
                .iter()
                .flatten()
                .filter_map(|&index| mesh.vertices[index].shared_id)
                .collect();

            for id in &network.lines[0].point_ids {
                assert!(used_ids.contains(id), "fracture {fracture}: point {id}");
            }
        }
    }

    #[test]
    fn test_no_angle_is_below_the_minimum_angle() {
        let (planes, network) = crossing_squares();

        for (fracture, plane) in planes.iter().enumerate() {
            let mesh = mesh_fracture(plane, fracture, &network).unwrap();

            for triangle in &mesh.triangles {
                let corners = triangle.map(|index| mesh.vertices[index].position);
                let smallest = (0..3)
                    .map(|i| {
                        let [vertex, next, previous] =
                            [corners[i], corners[(i + 1) % 3], corners[(i + 2) % 3]];
                        angle_between(
                            [next[0] - vertex[0], next[1] - vertex[1]],
                            [previous[0] - vertex[0], previous[1] - vertex[1]],
                        )
                    })
                    .fold(f64::INFINITY, f64::min);
                assert!(smallest.to_degrees() >= MIN_ANGLE - 1e-9, "{smallest}");
            }
        }
    }

    #[test]
    fn test_no_triangle_is_larger_than_the_largest_size() {
        let (planes, network) = crossing_squares();

        let mesh = mesh_fracture(&planes[0], 0, &network).unwrap();

        for triangle in &mesh.triangles {
            let [a, b, c] = triangle.map(|index| mesh.vertices[index].position);
            let product = distance(a, b) * distance(b, c) * distance(c, a);
            let circumradius = product / (4. * triangle_area(&mesh, triangle));
            assert!(circumradius <= REFINE_RADIUS_PART * SIZES[1] + 1e-9);
        }
    }
}
