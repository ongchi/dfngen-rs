use std::{
    fs::File,
    io::{self, BufWriter, Write},
};

use parry3d_f64::na::{center, distance, Point3, Vector3};
use tracing::info;

use crate::fracture::poly::{DomainBoundaryMode, Poly};
use crate::io::input::Input;
use crate::io::vtk::{
    write_fractures_vtk, write_rejected_fractures_vtk, write_removed_fractures_vtk,
};
use crate::structures::{IntersectionPoints, RejectedFracture, Stats};

fn create_file(path: String) -> io::Result<BufWriter<File>> {
    Ok(BufWriter::new(File::create(path)?))
}

/// Moves a point from the coordinates of the generation to those of the output
///
/// With a domain mesh, the generation has the center of the mesh at the origin,
/// and the output has the coordinates of the mesh file.
///
/// # Arguments
///
/// * `point` - Point in the coordinates of the generation
/// * `mesh_center` - Center of the domain mesh, if the domain has a mesh
pub fn to_output_coordinates(
    point: &Point3<f64>,
    mesh_center: Option<&Vector3<f64>>,
) -> Point3<f64> {
    match mesh_center {
        Some(mesh_center) => point + mesh_center,
        None => *point,
    }
}

/// Polygons which remain after isolated fracture removal, in output order.
fn final_polys<'a>(
    final_fractures: &'a [usize],
    accepted_poly: &'a [Poly],
) -> impl Iterator<Item = &'a Poly> {
    final_fractures.iter().map(|&i| &accepted_poly[i])
}

/// Writes all output
///
/// # Arguments
///
/// * `input` - Input structure, contains all user input
/// * `output` - Path to output folder
/// * `accepted_poly` - Vector of all accepted polygons
/// * `rejected_poly` - Vector of all rejected polygons
/// * `int_pts` - Vector of all intersections from accepted polygons
/// * `triple_points` - Vector of all accepted triple intersection points
/// * `pstats` - Stats structure, running program statistics
/// * `final_fractures` - sorted indices into the Poly array of accepted polgons
///   which remain after isolated fractures (polys) were removed
#[allow(clippy::too_many_arguments)]
pub fn write_output(
    input: &Input,
    output: &str,
    accepted_poly: &[Poly],
    rejected_poly: &[RejectedFracture],
    int_pts: &mut [IntersectionPoints],
    triple_points: &[Point3<f64>],
    pstats: &Stats,
    final_fractures: &[usize],
) -> io::Result<()> {
    info!("{}", output);
    std::fs::create_dir_all(output)?;

    let mesh_center = input.domain_mesh.as_ref().map(|mesh| mesh.center());
    let mesh_center = mesh_center.as_ref();

    let mut vtk_file = create_file(format!("{}/fractures.vtp", output))?;
    let boundary_face_counts: Vec<usize> = accepted_poly
        .iter()
        .map(|poly| find_boundary_face_lines_for_mode(poly, input).len())
        .collect();
    write_fractures_vtk(
        final_fractures,
        accepted_poly,
        &boundary_face_counts,
        mesh_center,
        &mut vtk_file,
    )?;
    vtk_file.flush()?;

    let removed_fractures: Vec<usize> = (0..accepted_poly.len())
        .filter(|i| final_fractures.binary_search(i).is_err())
        .collect();

    if !removed_fractures.is_empty() {
        let mut removed_file = create_file(format!("{}/fractures_removed.vtp", output))?;
        write_removed_fractures_vtk(
            &removed_fractures,
            accepted_poly,
            &boundary_face_counts,
            mesh_center,
            &mut removed_file,
        )?;
        removed_file.flush()?;
    }

    if input.output_rejected_fractures {
        let mut rejected_file = create_file(format!("{}/fractures_rejected.vtp", output))?;
        write_rejected_fractures_vtk(rejected_poly, mesh_center, &mut rejected_file)?;
        rejected_file.flush()?;
    }

    adjust_int_fract_ids(final_fractures, accepted_poly, int_pts);
    write_intersection_list(
        input,
        final_fractures,
        accepted_poly,
        int_pts,
        mesh_center,
        output,
    )?;
    write_rejection_stats(pstats, output)?;
    write_user_rejected_fracture_information(pstats, output)?;
    write_connectivity(final_fractures, accepted_poly, int_pts, output)?;
    write_rejects_per_attempt(pstats, output)?;

    if input.triple_intersections {
        info!("Writing Triple Intersection Points File");
        write_triple_pts(
            triple_points,
            final_fractures,
            accepted_poly,
            int_pts,
            mesh_center,
            output,
        )?;
    }

    Ok(())
}

/// Adjust the intersectins fracture numbers. The finalFracture list is the indexes of the final polygons.
/// If finalFractures = {4, 2, 8}, the fracture ID's must be 1, 2, 3 respectively. Easiest way
/// is to adjust the fracture ID's in the corresponding intersections first, before writing output files
/// Use the negative value to keep to not loose track of what fracture id's are what (prevent aliasing).
///
/// EXAMPLE: For triple intersections, each intersection lists three fractures which intersect.
/// Say fractures 1, 5, and 6 intersect and once adjusted the the fracture IDs become 1, 4, and 5.
/// We access the intersection structure through the polygons, so for this particular
/// intersection structure, it will be accessed three times (once for each polygon). If we adjust ID 5 to be id 4 during
/// the second access, then during the thrid access, ID 4 may be aliasing another ID and be adjusted again.
/// To prevent this, we use the negative value to prevent aliasing.
///
/// # Arguments
///
/// * `final_fractures` - Vector array of indices of fractures left after isolated fracture removal
/// * `all_polys` - Vector array of all accepted polygons
/// * `int_pts` - Vector array of all intersections
fn adjust_int_fract_ids(
    final_fractures: &[usize],
    all_polys: &[Poly],
    int_pts: &mut [IntersectionPoints],
) {
    for (i, &poly_idx) in final_fractures.iter().enumerate() {
        let new_id = -(i as isize + 1);

        for &int_idx in &all_polys[poly_idx].intersection_index {
            let intersection = &mut int_pts[int_idx];

            if intersection.fract1 == poly_idx as isize {
                intersection.fract1 = new_id;
            } else if intersection.fract2 == poly_idx as isize {
                intersection.fract2 = new_id;
            }
        }
    }
}

/// Returns the output ids `(own, other)` of the two fractures of an intersection,
/// where `own` is the fracture with output id `fracture_id`.
///
/// NOTE: At this stage in the program, the fracture ID's have been changed to
/// negative. See adjust_int_fract_ids() for more details.
fn intersecting_fracture_ids(
    intersection: &IntersectionPoints,
    fracture_id: usize,
) -> (isize, isize) {
    if -intersection.fract1 == fracture_id as isize {
        (-intersection.fract1, -intersection.fract2)
    } else {
        (-intersection.fract2, -intersection.fract1)
    }
}

/// Writes triple intersection points to file (triple_Points.dat)
///
/// # Arguments
///
/// * `triple_points` - Vector array of all triple intersection points
/// * `final_fractures` - Vector array of indices of fractures left after isolated fracture removal
/// * `accepted_poly` - Vector array of all accetped fractures
/// * `int_pts` - Vector array of all intersections
/// * `mesh_center` - Center of the domain mesh, if the domain has a mesh
/// * `output` - Path to output folder
fn write_triple_pts(
    triple_points: &[Point3<f64>],
    final_fractures: &[usize],
    accepted_poly: &[Poly],
    int_pts: &[IntersectionPoints],
    mesh_center: Option<&Vector3<f64>>,
    output: &str,
) -> io::Result<()> {
    let mut file = create_file(format!("{}/triple_points.dat", output))?;

    // Each triple point is referenced by several intersections, write it only once
    let mut triple_pts_list: Vec<usize> = final_polys(final_fractures, accepted_poly)
        .flat_map(|poly| &poly.intersection_index)
        .flat_map(|&int_idx| &int_pts[int_idx].triple_points_idx)
        .copied()
        .collect();
    triple_pts_list.sort();
    triple_pts_list.dedup();

    for idx in triple_pts_list {
        let point = to_output_coordinates(&triple_points[idx], mesh_center);
        writeln!(file, "{:.17} {:.17} {:.17}", point.x, point.y, point.z)?;
    }

    file.flush()
}

/// Write rejections.dat, rejection statistics
///
/// # Arguments
///
/// * `pstats` - Stats structure of program statistics
/// * `output` - Path to output folder
fn write_rejection_stats(pstats: &Stats, output: &str) -> io::Result<()> {
    info!("Writing Rejection Statistics File (rejections.dat)");
    let mut file = create_file(format!("{}/rejections.dat", output))?;
    let reasons = &pstats.rejection_reasons;

    writeln!(file, "Short Intersection: {}", reasons.short_intersection)?;
    writeln!(file, "Close to Node: {}", reasons.close_to_node)?;
    writeln!(file, "Close to Edge: {}", reasons.close_to_edge)?;
    writeln!(
        file,
        "Vertex Close to Edge: {}",
        reasons.close_point_to_edge
    )?;
    writeln!(file, "Outside of Domain: {}", reasons.outside)?;
    writeln!(file, "Triple Intersection: {}", reasons.triple)?;
    writeln!(
        file,
        "Intersections Too Close: {}",
        reasons.inter_close_to_inter
    )?;
    file.flush()
}

/// Write userFractureRejections.dat, the user defined fractures which were rejected
///
/// # Arguments
///
/// * `pstats` - Stats structure of program statistics
/// * `output` - Path to output folder
fn write_user_rejected_fracture_information(pstats: &Stats, output: &str) -> io::Result<()> {
    if pstats.rejected_user_fracture.is_empty() {
        return Ok(());
    }

    info!("Writing User Fracture Rejection File (userFractureRejections.dat)");
    let mut file = create_file(format!("{}/userFractureRejections.dat", output))?;
    writeln!(file, "Fracture id")?;

    for fracture_id in &pstats.rejected_user_fracture {
        writeln!(file, "{}", fracture_id)?;
    }

    file.flush()
}

/// Writes fracture connectivity edge graph
///
/// # Arguments
///
/// * `final_fractures` - Vector array of indices of fractures left after isolated fracture removal
/// * `accepted_poly` - Vector array of all accetped fractures
/// * `int_pts` - Vector array of all intersections
/// * `output` - Path to output folder
fn write_connectivity(
    final_fractures: &[usize],
    accepted_poly: &[Poly],
    int_pts: &[IntersectionPoints],
    output: &str,
) -> io::Result<()> {
    info!("Writing Connectivity Data (connectivity.dat)");
    let mut file = create_file(format!("{}/connectivity.dat", output))?;

    // Format for connectivity: Line number = poly number, integers on line are intersecting fractures
    //  eg:
    // 2 3
    // 1 3
    // 1 2
    // Fracture 1 intersects frac 2 and 3, fract 2 intersects 1 and 3 etc

    for (i, poly) in final_polys(final_fractures, accepted_poly).enumerate() {
        for &int_idx in &poly.intersection_index {
            let (_, other) = intersecting_fracture_ids(&int_pts[int_idx], i + 1);
            write!(file, "{} ", other)?;
        }

        writeln!(file)?;
    }

    file.flush()
}

/// Writes rejectsPerAttempt.dat
///
/// Outputs a file that contains a list of integers of the number
/// of attempts per fracture.
/// For example, if the 5th number in the list is 100, it means that
/// it took 100 fracture insertion attemps for before the 5th fracture
/// was accepted.
///
/// # Arguments
///
/// * `pstats` - Stats structure of program statistics
/// * `output` - Path to output folder
fn write_rejects_per_attempt(pstats: &Stats, output: &str) -> io::Result<()> {
    let mut file = create_file(format!("{}/rejectsPerAttempt.dat", output))?;
    info!("Writing Rotation Data File (rejectsPerAttempt.dat)");

    for rejects in &pstats.rejects_per_attempt {
        writeln!(file, "{}", rejects)?;
    }

    file.flush()
}

/// Domain boundary faces as `(axis, id of the + face, id of the - face)`.
/// The ids are the negated FEHM face numbers:
/// top / Z+ / 1, bottom / Z- / 2, left / X- / 3, front / Y+ / 4, right / X+ / 5, back / Y- / 6
const BOUNDARY_FACE_IDS: [(usize, isize, isize); 3] = [(0, -5, -3), (1, -4, -6), (2, -1, -2)];

/// Finds the lines where a fracture meets the domain boundary faces, if the
/// fractures are truncated
///
/// A whole fracture has no edge on a domain face, thus it has no such line.
fn find_boundary_face_lines_for_mode(poly: &Poly, input: &Input) -> Vec<(isize, usize, usize)> {
    match input.domain_boundary_mode {
        DomainBoundaryMode::Truncate => {
            find_boundary_face_lines(poly, &input.domain_size, input.eps)
        }
        DomainBoundaryMode::KeepWhole => Vec::new(),
    }
}

/// Finds the lines where a fracture meets the domain boundary faces
///
/// There is at most one line for each face: between the first vertex found
/// on the face and the first other vertex on the same face.
///
/// # Arguments
///
/// * `poly` - Fracture
/// * `domain_size` - Domain size in x, y, z
/// * `eps` - Epsilon value for floating point comparisons
///
/// # Returns
///
/// `(face id, index of the first vertex, index of the second vertex)` for each line
pub fn find_boundary_face_lines(
    poly: &Poly,
    domain_size: &Vector3<f64>,
    eps: f64,
) -> Vec<(isize, usize, usize)> {
    let half_domain = domain_size * 0.5;
    let number_of_nodes = poly.vertices.len();
    let mut is_face_found = [[false; 2]; 3];
    let mut lines = Vec::new();

    for k in 0..number_of_nodes {
        for (axis, plus_face_id, minus_face_id) in BOUNDARY_FACE_IDS {
            let is_on_plus_face =
                |vertex: usize| poly.vertices[vertex][axis] >= half_domain[axis] - eps;
            let is_on_minus_face =
                |vertex: usize| poly.vertices[vertex][axis] <= -half_domain[axis] + eps;

            let (side, face_id) = if is_on_plus_face(k) {
                (0, plus_face_id)
            } else if is_on_minus_face(k) {
                (1, minus_face_id)
            } else {
                continue;
            };

            if is_face_found[axis][side] {
                continue;
            }

            let other_vertex = (0..number_of_nodes).find(|&kk| {
                kk != k
                    && if side == 0 {
                        is_on_plus_face(kk)
                    } else {
                        is_on_minus_face(kk)
                    }
            });

            if let Some(kk) = other_vertex {
                lines.push((face_id, k, kk));
                is_face_found[axis][side] = true;
            }
        }
    }

    lines
}

/// Writes intersection_list.dat, the intersections between the final fractures
/// and between a final fracture and a domain boundary face
///
/// # Arguments
///
/// * `input` - Input structure, contains all user input
/// * `final_fractures` - Vector array of indices of fractures left after isolated fracture removal
/// * `accepted_poly` - Vector array of all accetped fractures (before isolated fracture removal)
/// * `int_pts` - Vector array of all intersections
/// * `mesh_center` - Center of the domain mesh, if the domain has a mesh
/// * `output` - Path to output folder
fn write_intersection_list(
    input: &Input,
    final_fractures: &[usize],
    accepted_poly: &[Poly],
    int_pts: &[IntersectionPoints],
    mesh_center: Option<&Vector3<f64>>,
    output: &str,
) -> io::Result<()> {
    info!("Writing Intersection List (intersection_list.dat)");
    let mut int_file = create_file(format!("{}/intersection_list.dat", output))?;
    writeln!(int_file, "f1 f2 x y z length")?;

    for (i, poly) in final_polys(final_fractures, accepted_poly).enumerate() {
        let fracture_id = i + 1;

        // Write each intersection once, from the fracture with the smaller id
        for &int_idx in &poly.intersection_index {
            let intersection = &int_pts[int_idx];
            let (fract1, fract2) = intersecting_fracture_ids(intersection, fracture_id);

            if fract1 < fract2 {
                write_mid_point(
                    &mut int_file,
                    fract1,
                    fract2,
                    &intersection.p1,
                    &intersection.p2,
                    mesh_center,
                )?;
            }
        }

        for (face_id, k, kk) in find_boundary_face_lines_for_mode(poly, input) {
            write_mid_point(
                &mut int_file,
                fracture_id as isize,
                face_id,
                &poly.vertices[k],
                &poly.vertices[kk],
                mesh_center,
            )?;
        }
    }

    int_file.flush()
}

/// Writes mid point and length of the line between two points into file fp
///
/// # Arguments
///
/// * `fp` - File to write information into
/// * `fract1` - Fracture 1
/// * `fract2` - Fracture 2
/// * `point1` - First endpoint
/// * `point2` - Second endpoint
/// * `mesh_center` - Center of the domain mesh, if the domain has a mesh
fn write_mid_point(
    fp: &mut impl Write,
    fract1: isize,
    fract2: isize,
    point1: &Point3<f64>,
    point2: &Point3<f64>,
    mesh_center: Option<&Vector3<f64>>,
) -> io::Result<()> {
    let mid_point = to_output_coordinates(&center(point1, point2), mesh_center);

    writeln!(
        fp,
        "{} {} {:.10} {:.10} {:.10} {:.10}",
        fract1,
        fract2,
        mid_point.x,
        mid_point.y,
        mid_point.z,
        distance(point1, point2)
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    fn write_mid_point_of_a_line_on_the_x_axis(mesh_center: Option<&Vector3<f64>>) -> String {
        let mut line = Vec::new();
        write_mid_point(
            &mut line,
            1,
            2,
            &Point3::new(1., 0., 0.),
            &Point3::new(3., 0., 0.),
            mesh_center,
        )
        .unwrap();
        String::from_utf8(line).unwrap()
    }

    #[test]
    fn test_mid_point_without_mesh_has_the_generation_coordinates() {
        assert_eq!(
            write_mid_point_of_a_line_on_the_x_axis(None),
            "1 2 2.0000000000 0.0000000000 0.0000000000 2.0000000000\n"
        );
    }

    #[test]
    fn test_mid_point_with_mesh_has_the_mesh_coordinates() {
        assert_eq!(
            write_mid_point_of_a_line_on_the_x_axis(Some(&Vector3::new(10., 20., 30.))),
            "1 2 12.0000000000 20.0000000000 30.0000000000 2.0000000000\n"
        );
    }

    #[test]
    fn test_fracture_inside_the_domain_meets_no_boundary_face() {
        let rectangle = Poly::new(4, std::f64::consts::SQRT_2, 2.);

        let lines = find_boundary_face_lines(&rectangle, &Vector3::new(10., 10., 10.), EPS);

        assert!(lines.is_empty());
    }

    #[test]
    fn test_fracture_with_two_vertices_on_a_face_meets_that_face() {
        // The rectangle has x = 1 at vertices 0 and 3, and x = -1 at vertices 1 and 2
        let rectangle = Poly::new(4, std::f64::consts::SQRT_2, 2.);

        let lines = find_boundary_face_lines(&rectangle, &Vector3::new(2., 10., 10.), EPS);

        assert_eq!(lines, [(-5, 0, 3), (-3, 1, 2)]);
    }

    #[test]
    fn test_whole_fracture_across_a_face_gives_no_boundary_face_line() {
        // The rectangle has x = 1 at vertices 0 and 3, and x = -1 at vertices 1 and 2
        let rectangle = Poly::new(4, std::f64::consts::SQRT_2, 2.);
        let input = Input {
            domain_size: Vector3::new(1., 10., 10.),
            eps: EPS,
            domain_boundary_mode: DomainBoundaryMode::KeepWhole,
            ..Default::default()
        };

        let lines = find_boundary_face_lines_for_mode(&rectangle, &input);

        assert!(lines.is_empty());
    }
}
