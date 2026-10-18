//! Triangle mesh of one fracture in its plane

use std::collections::{HashMap, HashSet};

use spade::handles::{FixedFaceHandle, InnerTag};
use spade::{
    mitigate_underflow, ConstrainedDelaunayTriangulation, HasPosition, InsertionError,
    Point2 as SpadePoint, Triangulation,
};

use super::plane::{point_on_segment, FracturePlane, Point2};
use super::shared_points::SharedNetwork;
use super::MERGE_DISTANCE_PART;
use crate::region_boundaries::{distance, distance_to_segment, position_on_segment, OtherSide};

/// A triangle gets a new point if its circumradius is above this part of the
/// target size. An equilateral triangle with sides of the target size has
/// 1/√3 ≈ 0.58.
const REFINE_RADIUS_PART: f64 = 0.65;

/// Part of the target size: a new point must be at least this far from the
/// lines and from the other points
const NEW_POINT_DISTANCE_PART: f64 = 0.4;

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

    let boundary_ends: Vec<Point2> = plane
        .boundary_segments
        .iter()
        .flat_map(|boundary| boundary.segment)
        .collect();

    for edge in plane.outline_edges() {
        if is_shared_edge(edge, &polylines, &vertices, merge_distance) {
            continue;
        }
        let polyline = divide_segment(edge, &boundary_ends, &shared_vertices, &mut vertices, plane);
        polylines.push(polyline);
    }

    let mut boundary_sides = Vec::new();
    for boundary in &plane.boundary_segments {
        let polyline = divide_segment(
            boundary.segment,
            &[],
            &shared_vertices,
            &mut vertices,
            plane,
        );
        boundary_sides.extend(
            polyline
                .iter()
                .map(|&index| (index, boundary.other_side.clone())),
        );
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

    refine(&mut cdt, &mut vertices, plane)?;

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

/// The vertices along a segment: its ends, the points on it (shared
/// vertices and the extra points), and points between them so that no part
/// is longer than the target size at its middle
fn divide_segment(
    segment: [Point2; 2],
    extra_points: &[Point2],
    shared_vertices: &[usize],
    vertices: &mut LocalVertices,
    plane: &FracturePlane,
) -> Vec<usize> {
    let merge_distance = MERGE_DISTANCE_PART * plane.min_size;

    let mut on_segment: Vec<(f64, usize)> = shared_vertices
        .iter()
        .map(|&index| (vertices.vertices[index].position, index))
        .filter(|(position, _)| distance_to_segment(*position, segment) < merge_distance)
        .map(|(position, index)| (position_on_segment(position, segment), index))
        .collect();

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

/// Adds points until no triangle in the outline is larger than the target
/// size: the circumcenter of a large triangle, or its centroid if the
/// circumcenter is outside or near a line. A point near a line or near
/// another point is not added, thus the lines are never divided again.
fn refine(
    cdt: &mut Cdt,
    vertices: &mut LocalVertices,
    plane: &FracturePlane,
) -> Result<(), InsertionError> {
    let is_valid = |point: Point2, size: f64| {
        plane.contains(point) && plane.distance_to_lines(point) >= NEW_POINT_DISTANCE_PART * size
    };

    for _ in 0..MAX_REFINEMENT_PASSES {
        let candidates: Vec<Point2> = cdt
            .inner_faces()
            .filter_map(|face| {
                let [a, b, c] = face.positions();
                let centroid = [(a.x + b.x + c.x) / 3., (a.y + b.y + c.y) / 3.];

                if !plane.contains(centroid) {
                    return None;
                }

                let circumcenter = face.circumcenter();
                let circumcenter = [circumcenter.x, circumcenter.y];
                let radius = distance(circumcenter, [a.x, a.y]);

                if radius <= REFINE_RADIUS_PART * plane.size_at(centroid) {
                    return None;
                }

                [circumcenter, centroid]
                    .into_iter()
                    .find(|&point| is_valid(point, plane.size_at(point)))
            })
            .collect();

        let mut added = 0;
        for point in candidates {
            let min_distance = NEW_POINT_DISTANCE_PART * plane.size_at(point);
            if vertices.nearest_within(point, min_distance).is_some() {
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
            added += 1;
        }

        if added == 0 {
            return Ok(());
        }
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
