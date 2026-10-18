//! Triangle mesh of the final fractures, the mesh file of fehm-rs
//!
//! Fractures which intersect share the points of their intersection line.
//! The triangle size is `min_feature_size` near the intersections and the
//! boundary lines of the regions, and grows with the distance from them.

mod plane;
mod shared_points;
mod triangulate;

use std::collections::HashMap;

use parry3d_f64::na::{center, distance, Point3};
use rand_mt::Mt64;
use tracing::{info, warn};

use crate::dfn::Dfn;
use crate::distribution::ApertureDistribution;
use crate::error::DfngenError;
use crate::fracture::poly::{FractureOrigin, Poly};
use crate::io::input::Input;
use crate::region_boundaries::OtherSide;
use plane::FracturePlane;
use shared_points::{build_shared_network, LineOnFractures, PointOnFractures};
use triangulate::mesh_fracture;

/// Part of `min_feature_size`: points nearer than this are one point. It
/// removes only the rounding differences of a point which is calculated
/// two times (in the planes of two fractures); it does not join features
/// which are near to each other (user-defined pieces can be thinner than
/// `min_feature_size`).
const MERGE_DISTANCE_PART: f64 = 1e-6;

/// Options of the mesh (input section `[mesh]`). The input has them only
/// when the option is on.
#[derive(Debug, Clone, PartialEq)]
pub struct MeshOptions {
    /// Largest triangle size, as a multiple of `min_feature_size`
    pub max_size_factor: f64,

    /// Growth of the triangle size with the distance from an intersection or
    /// a boundary line
    pub slope: f64,

    /// Aperture of each stochastic family, in the order of the families
    pub family_apertures: Vec<ApertureDistribution>,

    /// Aperture of the user-defined fractures. `None` if there are none.
    pub user_fracture_aperture: Option<ApertureDistribution>,
}

/// Triangle mesh of the final fractures
pub struct FractureMesh {
    /// Points in the coordinates of the generation
    pub points: Vec<Point3<f64>>,

    /// Indices into `points`, counter-clockwise around the normal of the fracture
    pub triangles: Vec<[usize; 3]>,

    /// Position of the fracture of each triangle in the final fractures
    pub triangle_fractures: Vec<usize>,

    /// Aperture of each final fracture (m). The pieces of a user-defined
    /// fracture have one aperture.
    pub apertures: Vec<f64>,

    /// For each point: what is on the other side of the boundary lines
    /// through it. Empty for a point which is on no boundary line.
    pub point_sides: Vec<Vec<OtherSide>>,
}

/// Meshes the final fractures.
///
/// # Arguments
///
/// * `input` - Input options
/// * `options` - Options of the mesh
/// * `dfn` - The generated network
/// * `final_fractures` - Sorted indices of the final fractures in `dfn.accepted_polys`
pub fn mesh_network(
    input: &Input,
    options: &MeshOptions,
    dfn: &Dfn,
    final_fractures: &[usize],
) -> Result<FractureMesh, DfngenError> {
    let polys: Vec<&Poly> = final_fractures
        .iter()
        .map(|&i| &dfn.accepted_polys[i])
        .collect();
    let groups = fracture_groups(&polys);
    let apertures = sample_apertures(&polys, &groups, options, input.seed)?;

    let mut final_position = vec![None; dfn.accepted_polys.len()];
    for (position, &index) in final_fractures.iter().enumerate() {
        final_position[index] = Some(position);
    }

    // An intersection of almost no length is a point where the two
    // fractures touch: one shared point keeps them connected
    let (intersections, touching): (Vec<LineOnFractures>, Vec<LineOnFractures>) = dfn
        .intersections
        .iter()
        .filter_map(|intersection| {
            Some(LineOnFractures {
                fractures: [
                    final_position[intersection.fracture1]?,
                    final_position[intersection.fracture2]?,
                ],
                ends: [intersection.p1, intersection.p2],
            })
        })
        .partition(|line| distance(&line.ends[0], &line.ends[1]) > input.eps);

    let sizes = [
        input.min_feature_size,
        options.max_size_factor * input.min_feature_size,
        options.slope,
    ];
    let planes: Vec<FracturePlane> = polys
        .iter()
        .enumerate()
        .map(|(position, poly)| {
            let ends: Vec<[Point3<f64>; 2]> = intersections
                .iter()
                .filter(|line| line.fractures.contains(&position))
                .map(|line| line.ends)
                .collect();
            FracturePlane::new(poly, &ends, sizes)
        })
        .collect();

    let (piece_edges, mut shared_points) = piece_contacts(&polys, &groups);
    shared_points.extend(touching.iter().map(|line| PointOnFractures {
        fractures: line.fractures.to_vec(),
        point: center(&line.ends[0], &line.ends[1]),
    }));
    let shared = build_shared_network(
        &planes,
        &groups,
        &intersections,
        &piece_edges,
        &shared_points,
    );

    let mut mesh = FractureMesh {
        points: shared.shared_points.points.clone(),
        triangles: Vec::new(),
        triangle_fractures: Vec::new(),
        apertures,
        point_sides: Vec::new(),
    };
    let mut point_sides: HashMap<usize, Vec<OtherSide>> = HashMap::new();
    let mut skipped_constraints = 0;

    for (position, plane) in planes.iter().enumerate() {
        let local = mesh_fracture(plane, position, &shared)
            .map_err(|error| DfngenError::Mesh(format!("fracture {}: {error}", position + 1)))?;
        skipped_constraints += local.skipped_constraints;

        let global_ids: Vec<usize> = local
            .vertices
            .iter()
            .map(|vertex| {
                vertex.shared_id.unwrap_or_else(|| {
                    mesh.points.push(plane.frame.to_space(vertex.position));
                    mesh.points.len() - 1
                })
            })
            .collect();

        for (index, other_side) in local.boundary_sides {
            let sides = point_sides.entry(global_ids[index]).or_default();
            if !sides.contains(&other_side) {
                sides.push(other_side);
            }
        }

        for triangle in local.triangles {
            mesh.triangles.push(triangle.map(|index| global_ids[index]));
            mesh.triangle_fractures.push(position);
        }
    }

    if skipped_constraints > 0 {
        warn!(
            "Mesh: {} edges of lines which cross other lines were not kept; the mesh does not follow these lines there",
            skipped_constraints
        );
    }

    mesh.point_sides = (0..mesh.points.len())
        .map(|id| point_sides.remove(&id).unwrap_or_default())
        .collect();

    let omitted_fractures = remove_parts_without_boundary(&mut mesh);
    if omitted_fractures > 0 {
        info!(
            "Mesh: {} fractures are not in the mesh: their clusters have no point on a boundary (fehm-rs needs one)",
            omitted_fractures
        );
    }

    remove_unused_points(&mut mesh);

    info!(
        "Mesh: {} points, {} triangles",
        mesh.points.len(),
        mesh.triangles.len()
    );

    Ok(mesh)
}

/// Group of each fracture: the position of the first piece of its
/// user-defined fracture, or its own position
fn fracture_groups(polys: &[&Poly]) -> Vec<usize> {
    let mut first_piece: HashMap<u64, usize> = HashMap::new();

    polys
        .iter()
        .enumerate()
        .map(|(position, poly)| match poly.user_fracture_id {
            Some(id) => *first_piece.entry(id).or_insert(position),
            None => position,
        })
        .collect()
}

/// The aperture of each final fracture. The apertures have their own
/// generator from the seed, thus the network does not change.
fn sample_apertures(
    polys: &[&Poly],
    groups: &[usize],
    options: &MeshOptions,
    seed: u64,
) -> Result<Vec<f64>, DfngenError> {
    let mut generator = Mt64::new(seed);
    let mut apertures: Vec<f64> = Vec::with_capacity(polys.len());

    for (position, poly) in polys.iter().enumerate() {
        if groups[position] != position {
            apertures.push(apertures[groups[position]]);
            continue;
        }

        let distribution = match poly.origin {
            FractureOrigin::Family(index) => options.family_apertures.get(index),
            FractureOrigin::UserDefined => options.user_fracture_aperture.as_ref(),
        }
        .ok_or_else(|| {
            DfngenError::Mesh(format!(
                "fracture {}: no aperture distribution",
                position + 1
            ))
        })?;

        apertures.push(distribution.sample(poly.x_radius, &mut generator));
    }

    Ok(apertures)
}

/// The common edges of two pieces of one user-defined fracture (same two
/// vertices), and the vertices of the pieces of fractures with two or more pieces
fn piece_contacts(
    polys: &[&Poly],
    groups: &[usize],
) -> (Vec<LineOnFractures>, Vec<PointOnFractures>) {
    let mut piece_count: HashMap<usize, usize> = HashMap::new();
    for &group in groups {
        *piece_count.entry(group).or_default() += 1;
    }

    let key_of = |point: &Point3<f64>| [point.x, point.y, point.z].map(f64::to_bits);
    let mut edges: HashMap<([u64; 3], [u64; 3]), usize> = HashMap::new();
    let mut piece_edges = Vec::new();
    let mut piece_vertices = Vec::new();

    for (position, poly) in polys.iter().enumerate() {
        if piece_count[&groups[position]] < 2 {
            continue;
        }

        let vertices = &poly.vertices;
        for (i, vertex) in vertices.iter().enumerate() {
            piece_vertices.push(PointOnFractures {
                fractures: vec![position],
                point: *vertex,
            });

            let next = &vertices[(i + 1) % vertices.len()];
            let key = if key_of(vertex) < key_of(next) {
                (key_of(vertex), key_of(next))
            } else {
                (key_of(next), key_of(vertex))
            };

            match edges.get(&key) {
                Some(&other) if groups[other] == groups[position] && other != position => {
                    piece_edges.push(LineOnFractures {
                        fractures: [other, position],
                        ends: [*vertex, *next],
                    });
                }
                _ => {
                    edges.insert(key, position);
                }
            }
        }
    }

    (piece_edges, piece_vertices)
}

/// True if the point is on a boundary of the mesh file: on the mesh surface,
/// or beside a cell of a region (a cell of no region has no array)
fn is_on_boundary(sides: &[OtherSide]) -> bool {
    sides.iter().any(|side| match side {
        OtherSide::MeshSurface => true,
        OtherSide::Cell(region_ids) => !region_ids.is_empty(),
    })
}

/// Removes the triangles of each connected part of the mesh which has no
/// point on a boundary. Returns the number of fractures which are then not
/// in the mesh.
fn remove_parts_without_boundary(mesh: &mut FractureMesh) -> usize {
    let mut parent: Vec<usize> = (0..mesh.points.len()).collect();
    fn root(parent: &mut [usize], mut point: usize) -> usize {
        while parent[point] != point {
            parent[point] = parent[parent[point]];
            point = parent[point];
        }
        point
    }

    for triangle in &mesh.triangles {
        for pair in [[triangle[0], triangle[1]], [triangle[1], triangle[2]]] {
            let (first, second) = (root(&mut parent, pair[0]), root(&mut parent, pair[1]));
            parent[first] = second;
        }
    }

    let mut has_boundary = vec![false; mesh.points.len()];
    for (point, sides) in mesh.point_sides.iter().enumerate() {
        if is_on_boundary(sides) {
            has_boundary[root(&mut parent, point)] = true;
        }
    }

    let fracture_count = |fractures: &[usize]| {
        let mut fractures = fractures.to_vec();
        fractures.sort_unstable();
        fractures.dedup();
        fractures.len()
    };
    let count_before = fracture_count(&mesh.triangle_fractures);

    let (triangles, fractures): (Vec<[usize; 3]>, Vec<usize>) = mesh
        .triangles
        .iter()
        .zip(&mesh.triangle_fractures)
        .filter(|(triangle, _)| has_boundary[root(&mut parent, triangle[0])])
        .map(|(triangle, fracture)| (*triangle, *fracture))
        .unzip();
    mesh.triangles = triangles;
    mesh.triangle_fractures = fractures;

    count_before - fracture_count(&mesh.triangle_fractures)
}

/// Removes the points which are in no triangle, and numbers the others again.
fn remove_unused_points(mesh: &mut FractureMesh) {
    let mut new_ids = vec![None; mesh.points.len()];
    for &id in mesh.triangles.iter().flatten() {
        new_ids[id] = Some(0);
    }

    for (count, new_id) in new_ids.iter_mut().flatten().enumerate() {
        *new_id = count;
    }

    let points = std::mem::take(&mut mesh.points);
    let point_sides = std::mem::take(&mut mesh.point_sides);
    for ((point, sides), new_id) in points.into_iter().zip(point_sides).zip(&new_ids) {
        if new_id.is_some() {
            mesh.points.push(point);
            mesh.point_sides.push(sides);
        }
    }

    for id in mesh.triangles.iter_mut().flatten() {
        *id = new_ids[*id].expect("each point of a triangle is used");
    }
}
