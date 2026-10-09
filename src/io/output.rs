use std::{
    fs::File,
    io::{self, BufWriter, Seek, Write},
};

use parry3d_f64::na::{center, distance, Point3, Vector3};
use tracing::info;

use crate::distribution::generating_points::discretize_line_of_intersection;
use crate::fracture::fracture_family::FractureFamilyCollection;
use crate::io::input::Input;
use crate::io::vtk::write_fractures_vtk;
use crate::math_functions::sorted_index;
use crate::structures::{IntersectionPoints, Stats};
use crate::{cg::poly_and_intersection_rotation_to_xy, fracture::poly::Poly};

const RADII_FORMAT_HEADER: &str =
    "(-2 = userPolygon, -1 = userRectangle, 0 = userEllipse, > 0 is family in order of famProb)";

fn create_file(path: String) -> io::Result<BufWriter<File>> {
    Ok(BufWriter::new(File::create(path)?))
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
/// * `output_folder` - Path to output folder
/// * `accepted_poly` - Vector of all accepted polygons
/// * `int_pts` - Vector of all intersections from accepted polygons
/// * `triple_points` - Vector of all accepted triple intersection points
/// * `pstats` - Stats structure, running program statistics
/// * `final_fractures` - sorted indices into the Poly array of accepted polgons
///   which remain after isolated fractures (polys) were removed
/// * `frac_families` - Family structure array of all stocastic families defined by user input
#[allow(clippy::too_many_arguments)]
pub fn write_output(
    input: &Input,
    output_folder: &str,
    accepted_poly: &mut [Poly],
    int_pts: &mut [IntersectionPoints],
    triple_points: &[Point3<f64>],
    pstats: &mut Stats,
    final_fractures: &[usize],
    frac_families: &FractureFamilyCollection,
) -> io::Result<()> {
    let output = format!("{}/dfnGen_output", output_folder);
    info!("{}", output);
    std::fs::create_dir_all(&output)?;

    let intersection_folder = format!("{}/intersections", output_folder);
    std::fs::create_dir_all(&intersection_folder)?;

    let poly_folder = format!("{}/polys", output_folder);
    std::fs::create_dir_all(&poly_folder)?;

    let radii_folder = format!("{}/radii/", output);
    if input.outputAcceptedRadiiPerFamily || input.outputFinalRadiiPerFamily {
        std::fs::create_dir_all(&radii_folder)?;
    }

    // Must come before write_intersection_files(), which rotates the polygons to the x-y plane
    let mut vtk_file = create_file(format!("{}/fractures.vtp", output))?;
    write_fractures_vtk(final_fractures, accepted_poly, &mut vtk_file)?;
    vtk_file.flush()?;

    adjust_int_fract_ids(final_fractures, accepted_poly, int_pts);
    write_graph_data(input, final_fractures, accepted_poly, int_pts, &output)?;
    write_polys(final_fractures, accepted_poly, &output)?;
    // Rotates polys to x-y plane, must come before write_polys_inp()
    write_intersection_files(
        input,
        final_fractures,
        accepted_poly,
        int_pts,
        triple_points,
        &intersection_folder,
        pstats,
    )?;
    write_polys_inp(final_fractures, accepted_poly, &poly_folder)?;
    write_params_file(input, final_fractures, pstats, &output)?;
    write_radii_file(final_fractures, accepted_poly, &output)?;
    write_rejection_stats(pstats, &output)?;
    write_user_rejected_fracture_information(pstats, &output)?;
    write_frac_fams(input, frac_families, &output)?;
    write_fracture_translations(final_fractures, accepted_poly, &output)?;
    write_connectivity(final_fractures, accepted_poly, int_pts, &output)?;
    write_rotation_data(
        input,
        accepted_poly,
        final_fractures,
        frac_families,
        &output,
    )?;
    write_normal_vectors(accepted_poly, final_fractures, &output)?;
    write_rejects_per_attempt(pstats, &output)?;
    write_final_poly_radii(final_fractures, accepted_poly, &output)?;
    write_final_poly_area(final_fractures, accepted_poly, &output)?;
    write_boundary_files(final_fractures, accepted_poly, &output)?;

    let user_files = &input.ext_fracture_files;
    // Fractures are marked -2 for user rects, -1 for user ellipses, -3 for user polygons
    let mut family_nums: Vec<isize> = (0..frac_families.families.len() as isize).collect();
    if user_files.user_rect_file.is_some() {
        family_nums.push(-2);
    }
    if user_files.user_ell_file.is_some() {
        family_nums.push(-1);
    }
    if user_files.user_poly_by_coord_file.is_some() {
        family_nums.push(-3);
    }

    if input.outputAcceptedRadiiPerFamily {
        info!("Writing Accepted Radii Files Per Family");

        for &family_num in &family_nums {
            write_family_radii(
                format!(
                    "{}/radii_AllAccepted_Fam_{}.dat",
                    radii_folder,
                    family_num + 1
                ),
                "Fracture Radii List Before Isolated Fracture and Cluster Removal",
                family_num,
                accepted_poly.iter(),
            )?;
        }
    }

    if input.outputFinalRadiiPerFamily {
        info!("Writing Final Radii Files Per Family");

        for family_num in family_nums {
            write_family_radii(
                format!("{}/radii_Final_Fam_{}.dat", radii_folder, family_num + 1),
                "Fracture Radii List After Isolated Fracture and Cluster Removal",
                family_num,
                final_polys(final_fractures, accepted_poly),
            )?;
        }
    }

    if input.tripleIntersections {
        info!("Writing Triple Intersection Points File");
        write_triple_pts(
            triple_points,
            final_fractures,
            accepted_poly,
            int_pts,
            &output,
        )?;
    }

    Ok(())
}

/// Helper function for writing discretized intersections
///
/// Function writes n points to file
///
/// # Arguments
///
/// * `output` - Output file stream object (file we are writing to)
/// * `points` - Vector of Point3<f64> array (discretized points)
/// * `start` - Index to point to start output to file.
///   When discretizing points, the program can discretize from an end
///   point to a triple intersection point, and then from the triple intersection point to
///   another end point. This means you have two arrays of points. The last point in the first
///   array will be the same as the first point in the second array. In this situation, you
///   would set start = 1 while writing the second array to avoid duplicate points in the output
/// * `count` - Counter of poitns written. Used to rember at which node a new intersection starts
fn write_points(
    output: &mut impl Write,
    points: &[Point3<f64>],
    start: usize,
    count: &mut usize,
) -> io::Result<()> {
    for point in points.iter().skip(start) {
        writeln!(
            output,
            "{:.12} {:.12} {:.12} {:.12}",
            count, point.x, point.y, point.z
        )?;
        *count += 1;
    }

    Ok(())
}

/// Helper function for writing discretized intersection points
/// Writes header and line connections after points have been written
///
/// # Arguments
///
/// * `fract_int_file` - Output file stream object (intersection file)
/// * `fract1` - Number, or index, of fracture whos intersection is being written
/// * `num_points` - Number of intersection points on fracture
/// * `num_intersections` - Number of intersections on fracture
/// * `int_start` - Vector array of node numbers which start an intersection.
///   Used to generate "line" connections in intersection inp files
/// * `intersecting_fractures` - Vector array of fracture id's (indices) who intersect fract1
fn finish_writing_int_file(
    fract_int_file: &mut BufWriter<File>,
    fract1: usize,
    num_points: usize,
    num_intersections: usize,
    int_start: &[usize],
    intersecting_fractures: &[usize],
) -> io::Result<()> {
    let mut count = 1;
    let mut idx = 0;

    // Lines
    for i in 0..num_points - num_intersections {
        if int_start[idx] == count + 1 {
            count += 1;
            idx += 1;
        }

        writeln!(
            fract_int_file,
            "{} {} line {} {}",
            i + 1,
            fract1,
            count,
            count + 1
        )?;
        count += 1;
    }

    writeln!(fract_int_file, "2 1 1")?;
    writeln!(fract_int_file, "a_b, integer")?;
    writeln!(fract_int_file, "b_a, integer")?;
    idx = 0;

    for i in 0..num_points {
        if int_start[idx] == i + 1 {
            idx += 1;
        }

        writeln!(
            fract_int_file,
            "{} {} {}",
            i + 1,
            fract1,
            intersecting_fractures[idx]
        )?;
    }

    // Move to header postion
    fract_int_file.rewind()?;
    write!(
        fract_int_file,
        "{} {} 2 0 0",
        num_points,
        num_points - num_intersections
    )?;
    fract_int_file.flush()
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

/// Writes intersection inp files to output folder
///
/// Rotates intersections, and triple intersection points to x-y plane
/// Also rotates polygons to x-y plane to save on computation during writePolysInp()
/// Rotating polygons here increases performance as we do not need to recalculate rotation matricies
///
/// # Arguments
///
/// * `input` - Input structure, contains all user input
/// * `final_fractures` - Vector array of indices of fractures left after isolated fracture removal
/// * `accepted_poly` - Vector array of all accetped fractures (before isolated fracture removal)
/// * `int_pts` - Vector array of all intersections
/// * `triple_points` - Vector array of all triple intersection points
/// * `intersection_folder` - Path to intersections folder
/// * `pstats` - Stats structure, running program statistics
fn write_intersection_files(
    input: &Input,
    final_fractures: &[usize],
    accepted_poly: &mut [Poly],
    int_pts: &[IntersectionPoints],
    triple_points: &[Point3<f64>],
    intersection_folder: &str,
    pstats: &mut Stats,
) -> io::Result<()> {
    let h = input.h;
    let visualization_mode = input.visualizationMode;

    // Keeps track of current un-rotated points we are working with
    let mut temp_point1;
    let mut temp_point2;
    info!("Writing Intersection Files");

    // Go through finalFractures. Rotate poly, intersections, and triple intersection points
    // to XY plane. Discretize and write to file
    for (i, &poly_idx) in final_fractures.iter().enumerate() {
        // Starting positions for each intersection. Lets us know how to make the line connections
        let mut int_start = Vec::new();
        let mut count = 1;
        // Counter for intersection header  number of points
        let mut num_int_pts = 0;
        // Used in writing which nodes belong to which two fractures in finishWritingOutput()
        let mut intersecting_fractures = Vec::new();
        let mut fract_int_file = create_file(format!(
            "{}/intersections_{}.inp",
            intersection_folder,
            i + 1
        ))?;
        // Buffer first line to leave space to write header
        writeln!(fract_int_file, "{:63}", "")?;
        // Go through each final fracture's intersections and write to output
        let size = accepted_poly[poly_idx].intersection_index.len();

        for j in 0..size {
            // tempTripPts holds rotated triple points for an intersection. Triple pts must be rotated 3 different
            // ways so we cannot change the original data
            let mut temp_trip_pts = Vec::new();
            // Used to measure current length before rotation (used to calculate number of points
            // to discretize) based on original intersection. This fixes any precision errors we
            // when calculating length after rotation, both rotations for the same intersection will
            // always have the same step size and same number of discretized points.
            let mut cur_length;
            let intersection = &int_pts[accepted_poly[poly_idx].intersection_index[j]];
            // Similarly to above, the intersection must be rotated two different ways,
            // one for each intersecting poly. We can't change the original data so we must use temp data
            let temp_intersection = poly_and_intersection_rotation_to_xy(
                intersection,
                &mut accepted_poly[poly_idx],
                triple_points,
                &mut temp_trip_pts,
                input.eps,
            );
            // poly and intersection now rotated
            let triple_pts_size = temp_trip_pts.len();

            let (_, fract2) = intersecting_fracture_ids(intersection, i + 1);
            intersecting_fractures.push(fract2 as usize);

            // If triple points exist on intersection, discretize from endpoint to closest triple point,
            // from triple to next triple point, and finally to other end point
            if triple_pts_size != 0 {
                // Keep track of number of triple points which will be in the
                // DFN (this is after isolated fracture removal)
                // NOTE: This will need to be divided by six to get correct value.
                // Division by six was determined through testing.
                pstats.triple_node_count += triple_pts_size;
                // Order the triple points by distances to know to discretize from point to next closest point
                let mut pt1 = temp_intersection.p1;
                temp_point1 = intersection.p1;

                // Distances from the first end point to the triple points
                let distances: Vec<f64> = temp_trip_pts
                    .iter()
                    .map(|triple_point| distance(&pt1, triple_point))
                    .collect();

                // Order the indices of the distances array shortest to largest distance
                // this lets us know which point to discritize to next
                let s = sorted_index(&distances);
                // Discretize from end point1 to first triple pt
                // pt1 already = enpoint1
                let mut pt2 = temp_trip_pts[s[0]];
                temp_point2 = triple_points[intersection.triple_points_idx[s[0]]];
                cur_length = distance(&temp_point1, &temp_point2);
                let mut points =
                    discretize_line_of_intersection(h, visualization_mode, &pt1, &pt2, cur_length);
                // Write points to file
                num_int_pts += points.len();
                write_points(&mut fract_int_file, &points, 0, &mut count)?;

                // If one trip pt, set up points to discretize from only triple pt to other end point
                if triple_pts_size == 1 {
                    pt1 = pt2;
                    pt2 = temp_intersection.p2;
                    temp_point1 = temp_point2;
                    temp_point2 = intersection.p2;
                } else {
                    // More than 1 triple point
                    for jj in 0..triple_pts_size - 1 {
                        pt1 = temp_trip_pts[s[jj]];
                        pt2 = temp_trip_pts[s[jj + 1]];
                        temp_point1 = triple_points[intersection.triple_points_idx[s[jj]]];
                        temp_point2 = triple_points[intersection.triple_points_idx[s[jj + 1]]];
                        cur_length = distance(&temp_point1, &temp_point2);
                        points = discretize_line_of_intersection(
                            h,
                            visualization_mode,
                            &pt1,
                            &pt2,
                            cur_length,
                        );
                        // Write points for first fracture to file, save second set of points to temp
                        num_int_pts += points.len() - 1;
                        write_points(&mut fract_int_file, &points, 1, &mut count)?;
                    }

                    // Set up points to go from last triple point to last endpoint
                    pt1 = pt2;
                    pt2 = temp_intersection.p2;
                    temp_point1 = temp_point2;
                    temp_point2 = intersection.p2;
                }

                cur_length = distance(&temp_point1, &temp_point2);
                points =
                    discretize_line_of_intersection(h, visualization_mode, &pt1, &pt2, cur_length);
                num_int_pts += points.len() - 1;
                write_points(&mut fract_int_file, &points, 1, &mut count)?;
            } else {
                // No triple intersection points on intersection line
                cur_length = distance(&intersection.p1, &intersection.p2);
                let points = discretize_line_of_intersection(
                    h,
                    visualization_mode,
                    &temp_intersection.p1,
                    &temp_intersection.p2,
                    cur_length,
                );
                num_int_pts += points.len();
                write_points(&mut fract_int_file, &points, 0, &mut count)?;
            }

            int_start.push(count);
        }

        // Done with fracture and intersections
        pstats.intersection_node_count += num_int_pts;
        // Write line connectivity and header
        finish_writing_int_file(
            &mut fract_int_file,
            i + 1,
            num_int_pts,
            size,
            &int_start,
            &intersecting_fractures,
        )?;
    }

    //Divide by 6 to remove the duplicate counts
    pstats.triple_node_count /= 6;

    Ok(())
}

/// Writes polygons.dat, the vertices of all final polygons
///
/// # Arguments
///
/// * `final_fractures` - Vector array of indices of fractures left after isolated fracture removal
/// * `accepted_poly` - Vector array of all accetped fractures
/// * `output` - Path to output folder
fn write_polys(final_fractures: &[usize], accepted_poly: &[Poly], output: &str) -> io::Result<()> {
    info!("Writing Polygon Files");
    let mut poly_output = create_file(format!("{}/polygons.dat", output))?;
    writeln!(poly_output, "nPolygons: {}", final_fractures.len())?;

    for poly in final_polys(final_fractures, accepted_poly) {
        write!(poly_output, "{} ", poly.vertices.len())?;

        for vertex in &poly.vertices {
            write!(
                poly_output,
                "{{{:.12}, {:.12}, {:.12}}} ",
                vertex.x, vertex.y, vertex.z
            )?;
        }

        writeln!(poly_output)?;
    }

    info!("Writing Polygon Files Complete");
    poly_output.flush()
}

/// Parses and writes all poly_x.inp files containing polygon (fracture) vertice and connectivity data
///
/// # Arguments
///
/// * `final_fractures` - Vector array of indices of fractures left after isolated fracture removal
/// * `accepted_poly` - Vector array of all accetped fractures
/// * `output` - Path to output folder
fn write_polys_inp(
    final_fractures: &[usize],
    accepted_poly: &[Poly],
    output: &str,
) -> io::Result<()> {
    info!("Writing poly inp files");

    for (j, poly) in final_polys(final_fractures, accepted_poly).enumerate() {
        let mut poly_output = create_file(format!("{}/poly_{}.inp", output, j + 1))?;
        let number_of_nodes = poly.vertices.len();

        // Write header
        writeln!(
            poly_output,
            "{} {} 0 0 0",
            number_of_nodes,
            number_of_nodes - 1
        )?;

        // Write vertices
        for (i, vertex) in poly.vertices.iter().enumerate() {
            writeln!(
                poly_output,
                "{:.12} {:.12} {:.12} {:.12}",
                i + 1,
                vertex.x,
                vertex.y,
                vertex.z
            )?;
        }

        // Write line connectivity
        for i in 1..number_of_nodes {
            writeln!(poly_output, "{} {} line {} {}", i, j + 1, i, i + 1)?;
        }

        poly_output.flush()?;
    }

    Ok(())
}

/// Writes params.txt
///
/// # Arguments
///
/// * `input` - Input structure, contains all user input
/// * `final_fractures` - Vector array of indices of fractures left after isolated fracture removal
/// * `pstats` - Stats structure, running program statistics
/// * `output` - Path to output folder
fn write_params_file(
    input: &Input,
    final_fractures: &[usize],
    pstats: &Stats,
    output: &str,
) -> io::Result<()> {
    let params_output_file = format!("{}/../params.txt", output);
    let mut params = create_file(params_output_file.clone())?;
    info!("Writing {}", &params_output_file);

    writeln!(params, "{}", final_fractures.len())?;
    writeln!(params, "{}", input.h)?;
    writeln!(params, "{}", input.visualizationMode as u8)?; // Production mode
    writeln!(
        params,
        "{}",
        pstats.intersection_node_count / 2 - pstats.triple_node_count
    )?;
    writeln!(params, "{}", input.domainSize[0])?;
    writeln!(params, "{}", input.domainSize[1])?;
    writeln!(params, "{}", input.domainSize[2])?;
    params.flush()
}

/// Line ending for a row of an all-accepted-fractures file.
/// Fractures removed due to isolation are marked with `R`.
///
/// `final_fractures` must be sorted.
fn removed_marker(final_fractures: &[usize], poly_index: usize) -> &'static str {
    if final_fractures.binary_search(&poly_index).is_ok() {
        ""
    } else {
        " R"
    }
}

/// Writes radii.dat (Radii Data)
///
/// # Arguments
///
/// * `final_fractures` - Vector array of indices of fractures left after isolated fracture removal
/// * `accepted_poly` - Vector array of all accetped fractures
/// * `output` - Path to output folder
fn write_radii_file(
    final_fractures: &[usize],
    accepted_poly: &[Poly],
    output: &str,
) -> io::Result<()> {
    info!("Writing Radii File (radii.dat)");
    let mut radii = create_file(format!("{}/radii.dat", output))?;
    writeln!(
        radii,
        "Format: xRadius yRadius Family# Removed {}",
        RADII_FORMAT_HEADER
    )?;

    for (i, poly) in accepted_poly.iter().enumerate() {
        writeln!(
            radii,
            "{:.8} {:.8} {:.8}{}",
            poly.xradius,
            poly.yradius,
            poly.family_num + 1,
            removed_marker(final_fractures, i)
        )?;
    }

    radii.flush()
}

/// Wrtes translations file (translations.dat)
///
/// # Arguments
///
/// * `final_fractures` - Vector array of indices of fractures left after isolated fracture removal
/// * `accepted_poly` - Vector array of all accetped fractures
/// * `output` - Path to output folder
fn write_fracture_translations(
    final_fractures: &[usize],
    accepted_poly: &[Poly],
    output: &str,
) -> io::Result<()> {
    info!("Writing Fracture Translations File (translations.dat)");
    let mut file = create_file(format!("{}/translations.dat", output))?;
    writeln!(
        file,
        "Format: x y z  (R = removed from domain due to fracture isolation)"
    )?;

    for (i, poly) in accepted_poly.iter().enumerate() {
        writeln!(
            file,
            "{:.10} {:.10} {:.10}{}",
            poly.translation[0],
            poly.translation[1],
            poly.translation[2],
            removed_marker(final_fractures, i)
        )?;
    }

    file.flush()
}

/// Deprecated Function
/// Writes final radii file (after isoloated fractures have been removed)
///
/// # Arguments
///
/// * `final_fractures` - Vector array of indices of fractures left after isolated fracture removal
/// * `accepted_poly` - Vector array of all accetped fractures
/// * `output` - Path to output folder
fn write_final_poly_radii(
    final_fractures: &[usize],
    accepted_poly: &[Poly],
    output: &str,
) -> io::Result<()> {
    let mut radii_final = create_file(format!("{}/radii_Final.dat", output))?;
    writeln!(
        radii_final,
        "Fracture Radii List After Isolated Fracture and Cluster Removal"
    )?;
    writeln!(
        radii_final,
        "Format: xRadius yRadius Family# {}",
        RADII_FORMAT_HEADER
    )?;

    for poly in final_polys(final_fractures, accepted_poly) {
        writeln!(
            radii_final,
            "{} {} {}",
            poly.xradius,
            poly.yradius,
            poly.family_num + 1
        )?;
    }

    radii_final.flush()
}

/// Deprecated Function
/// Writes final surface area file (after isoloated fractures have been removed)
///
/// # Arguments
///
/// * `final_fractures` - Vector array of indices of fractures left after isolated fracture removal
/// * `accepted_poly` - Vector array of all accetped fractures
/// * `output` - Path to output folder
fn write_final_poly_area(
    final_fractures: &[usize],
    accepted_poly: &[Poly],
    output: &str,
) -> io::Result<()> {
    let mut area_final = create_file(format!("{}/surface_area_Final.dat", output))?;
    writeln!(
        area_final,
        "Fracture Surface Area After Isolated Fracture and Cluster Removal"
    )?;

    for poly in final_polys(final_fractures, accepted_poly) {
        writeln!(area_final, "{}", poly.area)?;
    }

    area_final.flush()
}

/// Writes a radii file for the fractures of one family
///
/// # Arguments
///
/// * `path` - Path of the file to write
/// * `title` - First line of the file, without the family number
/// * `family_num` - Family number for which radii file will be written for
///   -2 - User Rectangles, -1 - User Ellipses, Family# >= 0 - Family in order of 'FamProb' in input file
/// * `polys` - Fractures to select from
fn write_family_radii<'a>(
    path: String,
    title: &str,
    family_num: isize,
    polys: impl Iterator<Item = &'a Poly>,
) -> io::Result<()> {
    let mut file = create_file(path)?;
    writeln!(file, "{} (Family {})", title, family_num + 1)?;
    writeln!(
        file,
        "Format: xRadius yRadius Distrubution# {}",
        RADII_FORMAT_HEADER
    )?;

    for poly in polys.filter(|poly| poly.family_num == family_num) {
        writeln!(
            file,
            "{} {} {}",
            poly.xradius,
            poly.yradius,
            poly.family_num + 1
        )?;
    }

    file.flush()
}

/// Writes triple intersection points to file (triple_Points.dat)
///
/// # Arguments
///
/// * `triple_points` - Vector array of all triple intersection points
/// * `final_fractures` - Vector array of indices of fractures left after isolated fracture removal
/// * `accepted_poly` - Vector array of all accetped fractures
/// * `int_pts` - Vector array of all intersections
/// * `output` - Path to output folder
fn write_triple_pts(
    triple_points: &[Point3<f64>],
    final_fractures: &[usize],
    accepted_poly: &[Poly],
    int_pts: &[IntersectionPoints],
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
        let point = &triple_points[idx];
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
    writeln!(file, "Fracture id,User Fracture Type")?;

    for rejected in &pstats.rejected_user_fracture {
        writeln!(file, "{},{}", rejected.id, rejected.user_fracture_type)?;
    }

    file.flush()
}

/// Writes families.dat, Shape families definition file
///
/// # Arguments
///
/// * `input` - Input structure, contains all user input
/// * `frac_families` - Vector array of all fracture shape families
/// * `output` - Path to output folder
fn write_frac_fams(
    input: &Input,
    frac_families: &FractureFamilyCollection,
    output: &str,
) -> io::Result<()> {
    let rad_to_deg = 180. / std::f64::consts::PI;
    info!("Writing Family Definitions File (families.dat)");
    let mut file = create_file(format!("{}/families.dat", output))?;
    let user_files = &input.ext_fracture_files;

    //TODO: add stub code in families.dat for userDefined fractures, IF there are user defined fractures

    if user_files.user_ell_file.is_some() {
        write!(file, "UserDefined Ellipse Family: 0\n\n")?;
    }

    if user_files.user_rect_file.is_some() {
        write!(file, "UserDefined Rectangle Family: -1\n\n")?;
    }

    if user_files.user_poly_by_coord_file.is_some() {
        write!(file, "UserDefined Polygon Family: -2\n\n")?;
    }

    for (i, shape) in frac_families.families.iter().enumerate() {
        //name(rect or ell) and number of family
        writeln!(file, "{} Family: {}", shape.shape, shape.id)?;
        writeln!(file, "Global Family: {}", i + 1)?;
        writeln!(
            file,
            "Number of Vertices: {}",
            shape.shape.number_of_nodes()
        )?;
        writeln!(file, "Aspect Ratio: {}", shape.aspect_ratio)?;

        if input.stopCondition == 1 {
            writeln!(
                file,
                "P32 (Fracture Intensity) Target: {}",
                shape.p32_target
            )?;
        }

        // beta distribution, rotation around normal vector
        if !shape.beta_distribution {
            writeln!(
                file,
                "Beta Distribution (Rotation Around Normal Vector): Uniform"
            )?;
        } else {
            writeln!(
                file,
                "Beta (Rotation Around Normal Vector)-rad: {}",
                shape.beta
            )?;
            writeln!(
                file,
                "Beta (Rotation Around Normal Vector)-deg: {}",
                shape.beta * rad_to_deg
            )?;
        }

        write!(file, "{}", shape.orientation)?;

        // Print layer family belongs to
        if shape.layer_id == 0 {
            writeln!(file, "Layer: Entire domain")?;
        } else {
            let idx = (shape.layer_id - 1) * 2;
            writeln!(file, "Layer Number: {}", shape.layer_id)?;
            writeln!(
                file,
                "Layer: {{{}, {}}}",
                input.layers[idx],
                input.layers[idx + 1]
            )?;
        }

        // Print region family belongs to
        if shape.region_id == 0 {
            writeln!(file, "Region: Entire domain")?;
        } else {
            let idx = (shape.region_id - 1) * 6;
            let region = &input.regions[idx..idx + 6];
            writeln!(file, "Region Number: {}", shape.region_id)?;
            writeln!(
                file,
                "Region: {{{}, {}, {}, {}, {}, {}}}",
                region[0], region[1], region[2], region[3], region[4], region[5]
            )?;
        }

        // Print distribution data
        write!(file, "{}", shape.radius)?;

        write!(
            file,
            "Family Insertion Probability: {}\n\n",
            frac_families.original_probabilities[i]
        )?;
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

/// Writes poly_info.dat
/// Writes fracture rotation data. Also includes shape families each fracture belongs to.
///
/// # Arguments
///
/// * `input` - Input structure, contains all user input
/// * `accepted_poly` - Array off all polygons in DFN before isolated fracture removal
/// * `final_fractures` - Array of indices to polys in 'accepted_poly' which are left after isolated fracture removal
/// * `frac_families` - Array of all fracture shape families
/// * `output` - Path to output folder
fn write_rotation_data(
    input: &Input,
    accepted_poly: &[Poly],
    final_fractures: &[usize],
    frac_families: &FractureFamilyCollection,
    output: &str,
) -> io::Result<()> {
    let mut file = create_file(format!("{}/../poly_info.dat", output))?;
    info!("Writing Rotation Data File (poly_info.dat)");

    let domain_size = &input.domainSize;
    let max_domain_size = 10. * domain_size[0].max(domain_size[1]).max(domain_size[2]);
    let total_families = frac_families.families.len() as isize;
    let e3 = Vector3::new(0., 0., 1.);

    for (i, poly) in final_polys(final_fractures, accepted_poly).enumerate() {
        // poly's normal is already normalized at this point
        let normal = poly.normal;
        // Rotation angle in degrees
        let theta = normal.dot(&e3).acos() * (180.0 / std::f64::consts::PI);
        // Rotation into xy plane
        let mut v = e3.cross(&normal);

        if !(v[0].abs() < input.eps && v[1].abs() < input.eps && v[2].abs() < input.eps) {
            //if not zero vector
            v = v.normalize();
        }

        let x0 = 1.1 * (-max_domain_size * v[0]);
        let y0 = 1.1 * (-max_domain_size * v[1]);
        let z0 = 1.1 * (-max_domain_size * v[2]);
        let x1 = 1.1 * max_domain_size * v[0];
        let y1 = 1.1 * max_domain_size * v[1];
        let z1 = 1.1 * max_domain_size * v[2];
        // The last number is the shape family number. Throughout this program,
        // -1 and -2 are used to denote user defined ell. and rectangles.
        // -1 and -2 are not good numbers to use as material IDs in Lagrit.
        // If the family number is -1 or -2 we change these numbers to be
        // number of stochastic families + 1 and number of stochastic families + 2
        let fam_num = match poly.family_num {
            -1 => total_families + 1,
            -2 => total_families + 2,
            family_num => family_num + 1,
        };

        // Format: fracture#, family#, theta, x0, y0, z0, x1, y1, z1
        writeln!(
            file,
            "{} {} {:.15} {:.15} {:.15} {:.15} {:.15} {:.15} {:.15}",
            i + 1,
            fam_num,
            theta,
            x0,
            y0,
            z0,
            x1,
            y1,
            z1
        )?;
    }

    file.flush()
}

/// Writes normal_vectors.dat
///
/// # Arguments
///
/// * `accepted_poly` - Array off all polygons in DFN before isolated fracture removal
/// * `final_fractures` - Array of indices to polys in 'accepted_poly' which are left after isolated fracture removal
/// * `output` - Path to output folder
fn write_normal_vectors(
    accepted_poly: &[Poly],
    final_fractures: &[usize],
    output: &str,
) -> io::Result<()> {
    let mut file = create_file(format!("{}/normal_vectors.dat", output))?;
    info!("Writing Normal Vectors into File (normal_vectors.dat)");

    for poly in final_polys(final_fractures, accepted_poly) {
        // poly's normal is already normalized at this point
        // Format: nx, ny, nz
        writeln!(
            file,
            " {:.15} {:.15} {:.15}",
            poly.normal[0], poly.normal[1], poly.normal[2]
        )?;
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

/// Writes graph data files to intersections_list.dat and fracture_info.dat
///
/// # Arguments
///
/// * `input` - Input structure, contains all user input
/// * `final_fractures` - Vector array of indices of fractures left after isolated fracture removal
/// * `accepted_poly` - Vector array of all accetped fractures (before isolated fracture removal)
/// * `int_pts` - Vector array of all intersections
/// * `output` - Path to output folder
fn write_graph_data(
    input: &Input,
    final_fractures: &[usize],
    accepted_poly: &[Poly],
    int_pts: &[IntersectionPoints],
    output: &str,
) -> io::Result<()> {
    let eps = input.eps;
    let half_domain = [
        input.domainSize[0] * 0.5,
        input.domainSize[1] * 0.5,
        input.domainSize[2] * 0.5,
    ];

    info!("Writing Graph Data Files");
    let mut int_file = create_file(format!("{}/intersection_list.dat", output))?;
    writeln!(int_file, "f1 f2 x y z length")?;
    let mut fract_file = create_file(format!("{}/fracture_info.dat", output))?;
    writeln!(fract_file, "num_connections perm aperture")?;

    for (i, poly) in final_polys(final_fractures, accepted_poly).enumerate() {
        let fracture_id = i + 1;
        let mut num_conn = 0;

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
                )?;
                num_conn += 1;
            }
        }

        // Find intersections with domain boundaries.
        // A line of intersection is written only once per face: between the first
        // vertex found on the face and the first other vertex on the same face.
        let number_of_nodes = poly.vertices.len();
        let mut is_face_written = [[false; 2]; 3];

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

                if is_face_written[axis][side] {
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
                    write_mid_point(
                        &mut int_file,
                        fracture_id as isize,
                        face_id,
                        &poly.vertices[k],
                        &poly.vertices[kk],
                    )?;
                    is_face_written[axis][side] = true;
                    num_conn += 1;
                }
            }
        }

        writeln!(fract_file, "{} 0 0", num_conn)?;
    }

    int_file.flush()?;
    fract_file.flush()
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
fn write_mid_point(
    fp: &mut impl Write,
    fract1: isize,
    fract2: isize,
    point1: &Point3<f64>,
    point2: &Point3<f64>,
) -> io::Result<()> {
    let mid_point = center(point1, point2);

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

/// Writes fracture numbers into ASCII files corresponding to which boundary they touch
///
/// # Arguments
///
/// * `final_fractures` - Vector array of indices of fractures left after isolated fracture removal
/// * `accepted_poly` - Vector array of all accepted fractures (before isolated fracture removal)
fn write_boundary_files(
    final_fractures: &[usize],
    accepted_poly: &[Poly],
    output: &str,
) -> io::Result<()> {
    info!("Writing Boundary Files");

    // File name for each index of Poly::faces
    let boundary_files = ["right", "left", "front", "back", "top", "bottom"];

    for (face, name) in boundary_files.iter().enumerate() {
        let mut file = create_file(format!("{}/{}.dat", output, name))?;

        for (i, poly) in final_polys(final_fractures, accepted_poly).enumerate() {
            if poly.faces[face] {
                writeln!(file, "{}", i + 1)?;
            }
        }

        file.flush()?;
    }

    Ok(())
}
