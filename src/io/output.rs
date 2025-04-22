use std::{
    fs::File,
    io::{Seek, Write},
};

use parry3d_f64::na::{distance, Point3, Vector3};
use tracing::info;

use super::input::Input;
use crate::{
    computational_geometry::poly_and_intersection_rotation_to_xy,
    distribution::generating_points::discretize_line_of_intersection,
    fracture::insert_shape::get_family_number,
    math_functions::sorted_index,
    structures::{FractureFamilyOption, IntersectionPoints, Poly, Shape, Stats},
};

/// Writes all output for DFNGen
///
/// # Arguments
///
/// * `input` - Input structure, contains all user input
/// * `output_folder` - Path to output folder
/// * `accepted_poly` - Vector of all accepted polygons
/// * `int_pts` - Vector of all intersections from accepted polygons
/// * `triple_points` - Vector of all accepted triple intersection points
/// * `pstats` - Stats structure, running program statistics
/// * `final_fractures` - indices into the Poly array of accepted polgons
///     which remain after isolated fractures (polys) were removed
/// * `frac_families` - Family structure array of all stocastic families defined by user input
#[allow(clippy::too_many_arguments)]
pub fn write_output(
    input: &Input,
    output_folder: &str,
    accepted_poly: &mut [Poly],
    int_pts: &mut [IntersectionPoints],
    triple_points: &mut [Point3<f64>],
    pstats: &mut Stats,
    final_fractures: &[usize],
    frac_families: &FractureFamilyOption,
) {
    let output = format!("{}/dfnGen_output", output_folder);
    info!("{}", output);
    let _ = std::fs::create_dir_all(&output);

    let intersection_folder = format!("{}/intersections", output_folder);
    let _ = std::fs::create_dir_all(&intersection_folder);

    let poly_folder = format!("{}/polys", output_folder);
    let _ = std::fs::create_dir_all(&poly_folder);

    let radii_folder = format!("{}/radii/", output);

    // Adjust Fracture numbering
    adjust_int_fract_ids(final_fractures, accepted_poly, int_pts);
    // Write out graph information
    write_graph_data(
        input.eps,
        &input.domainSize,
        final_fractures,
        accepted_poly,
        int_pts,
        &output,
    );
    // Write polygon.dat file
    write_polys(final_fractures, accepted_poly, &output);
    // Write intersection files (must be first file written, rotates polys to x-y plane)
    write_intersection_files(
        input.h,
        input.eps,
        input.keepIsolatedFractures,
        input.visualizationMode,
        final_fractures,
        accepted_poly,
        int_pts,
        triple_points,
        &intersection_folder,
        pstats,
    );
    // Write polys.inp
    write_polys_inp(final_fractures, accepted_poly, &poly_folder);
    // Write params.txt
    write_params_file(
        input.h,
        input.visualizationMode,
        &input.domainSize,
        final_fractures,
        pstats,
        &output,
    );
    // Write aperture file
    // writeApertureFile(finalFractures, acceptedPoly, output);
    // Write permability file
    // writePermFile(finalFractures, acceptedPoly, output);
    // Write radii file
    write_radii_file(final_fractures, accepted_poly, &output);
    // Write rejection stats file
    write_rejection_stats(pstats, &output);
    // write out userRejetedFracture information
    write_user_rejected_fracture_information(pstats, &output);
    // Write families to output Files
    write_frac_fams(
        input.userEllipsesOnOff,
        input.userRectanglesOnOff,
        input.userPolygonByCoord,
        input.stopCondition,
        input.nFamEll,
        &input.layers,
        &input.regions,
        frac_families,
        &output,
    );
    // Write fracture translations file
    write_fracture_translations(final_fractures, accepted_poly, &output);
    // Write fracture connectivity (edge graph) file
    write_connectivity(final_fractures, accepted_poly, int_pts, &output);
    // Write rotation data
    write_rotation_data(
        input.eps,
        &input.domainSize,
        accepted_poly,
        final_fractures,
        frac_families,
        &output,
    );
    // Write normal vectors
    write_normal_vectors(accepted_poly, final_fractures, &output);
    // Write rejects per fracture insertion attempt data
    write_rejects_per_attempt(pstats, &output);
    // Write all accepted radii
    write_final_poly_radii(final_fractures, accepted_poly, &output);
    // Write all accepted Surface Area
    write_final_poly_area(final_fractures, accepted_poly, &output);
    // Write out which fractures touch which boundaries
    write_boundary_files(final_fractures, accepted_poly, &output);
    if input.outputAcceptedRadiiPerFamily {
        info!("Writing Accepted Radii Files Per Family");
        // Creates radii files per family, before isolated fracture removal.
        let size = frac_families.families.len();

        for i in 0..size {
            write_all_accepted_radii_of_family(i as isize, accepted_poly, &radii_folder);
        }

        if input.userRectanglesOnOff {
            // Fractures are marked -2 for user rects
            write_all_accepted_radii_of_family(-2, accepted_poly, &radii_folder);
        }

        if input.userEllipsesOnOff {
            // Fractures are marked -1 for user ellipses
            write_all_accepted_radii_of_family(-1, accepted_poly, &radii_folder);
        }

        if input.userPolygonByCoord {
            // Fractures are marked -3 for user user polygons
            write_all_accepted_radii_of_family(-3, accepted_poly, &radii_folder);
        }
    }

    if input.outputFinalRadiiPerFamily {
        info!("Writing Final Radii Files Per Family");
        let size = frac_families.families.len();

        for i in 0..size {
            write_final_radii_of_family(final_fractures, i as isize, accepted_poly, &radii_folder);
        }

        if input.userRectanglesOnOff {
            write_final_radii_of_family(final_fractures, -1, accepted_poly, &radii_folder);
        }

        if input.userEllipsesOnOff {
            write_final_radii_of_family(final_fractures, -2, accepted_poly, &radii_folder);
        }

        if input.userPolygonByCoord {
            write_final_radii_of_family(final_fractures, -3, accepted_poly, &radii_folder);
        }
    }

    // If triple intersections are on, write triple intersection points file
    if input.tripleIntersections {
        info!("Writing Triple Intersection Points File");
        write_triple_pts(
            triple_points,
            final_fractures,
            accepted_poly,
            int_pts,
            &output,
        );
    }
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
///     When discretizing points, the program can discretize from an end
///     point to a triple intersection point, and then from the triple intersection point to
///     another end point. This means you have two arrays of points. The last point in the first
///     array will be the same as the first point in the second array. In this situation, you
///     would set start = 1 while writing the second array to avoid duplicate points in the output
/// * `count` - Counter of poitns written. Used to rember at which node a new intersection starts
fn write_points(output: &mut File, points: &[Point3<f64>], start: usize, count: &mut usize) {
    let n = points.len();

    for point in points.iter().take(n).skip(start) {
        output
            .write_all(
                format!(
                    "{:.12} {:.12} {:.12} {:.12}\n",
                    count, point.x, point.y, point.z
                )
                .as_bytes(),
            )
            .unwrap();
        *count += 1;
    }
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
///     Used to generate "line" connections in intersection inp files
/// * `intersecting_fractures` - Vector array of fracture id's (indices) who intersect fract1
fn finish_writing_int_file(
    fract_int_file: &mut File,
    fract1: usize,
    num_points: usize,
    num_intersections: usize,
    int_start: &[usize],
    intersecting_fractures: &[usize],
) {
    let mut count = 1;
    let mut idx = 0;

    // Lines
    for i in 0..num_points - num_intersections {
        if int_start[idx] == count + 1 {
            count += 1;
            idx += 1;
        }

        fract_int_file
            .write_all(format!("{} {} line {} {}\n", i + 1, fract1, count, count + 1).as_bytes())
            .unwrap();
        count += 1;
    }

    fract_int_file.write_all("2 1 1\n".as_bytes()).unwrap();
    fract_int_file
        .write_all("a_b, integer\n".as_bytes())
        .unwrap();
    fract_int_file
        .write_all("b_a, integer\n".as_bytes())
        .unwrap();
    idx = 0;

    for i in 0..num_points {
        if int_start[idx] == i + 1 {
            idx += 1;
        }

        fract_int_file
            .write_all(format!("{} {} {}\n", i + 1, fract1, intersecting_fractures[idx]).as_bytes())
            .unwrap();
    }

    // Move to header postion
    fract_int_file.rewind().unwrap();
    fract_int_file
        .write_all(format!("{} {} 2 0 0", num_points, num_points - num_intersections).as_bytes())
        .unwrap();
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
    //go through all final fractures
    for i in 0..final_fractures.len() {
        //go through each final fractures intersections
        for j in 0..all_polys[final_fractures[i]].intersection_index.len() {
            //change matching fracture numbers to new number (order of finalFractures list) for output
            let int_idx = all_polys[final_fractures[i]].intersection_index[j];

            if int_pts[int_idx].fract1 == final_fractures[i] as isize {
                int_pts[int_idx].fract1 = -(i as isize + 1);
            } else if int_pts[int_idx].fract2 == final_fractures[i] as isize {
                int_pts[int_idx].fract2 = -(i as isize + 1);
            }
        }
    }
}

// Writes intersection inp files to output folder
//
// Rotates intersections, and triple intersection points to x-y plane
// Also rotates polygons to x-y plane to save on computation during writePolysInp()
// Rotating polygons here increases performance as we do not need to recalculate rotation matricies
//
// # Arguments
//
/// * `h` - Minimum feature size
/// * `eps` - Epsilon value for floating point comparison
/// * `keep_isolated_fractures` - Flag to keep isolated fractures
/// * `visualization_mode` - If false, creates a fine mesh, according to h parameter.
///     If true, produce only first round of triangulations. In this case no modeling of flow and transport is possible.
/// * `final_fractures` - Vector array of indices of fractures left after isolated fracture removal
/// * `accepted_poly` - Vector array of all accetped fractures (before isolated fracture removal)
/// * `int_pts` - Vector array of all intersections
/// * `triple_points` - Vector array of all triple intersection points
/// * `intersection_folder` - Path to intersections folder
/// * `pstats` - Stats structure, running program statistics
#[allow(clippy::too_many_arguments)]
fn write_intersection_files(
    h: f64,
    eps: f64,
    keep_isolated_fractures: bool,
    visualization_mode: bool,
    final_fractures: &[usize],
    accepted_poly: &mut [Poly],
    int_pts: &[IntersectionPoints],
    triple_points: &mut [Point3<f64>],
    intersection_folder: &str,
    pstats: &mut Stats,
) {
    // Keeps track of current un-rotated points we are working with
    let mut temp_point1;
    let mut temp_point2;
    info!("Writing Intersection Files");
    // let mut fractIntFile = File::open();

    // Go through finalFractures. Rotate poly, intersections, and triple intersection points
    // to XY plane. Discretize and write to file
    for i in 0..final_fractures.len() {
        // Starting positions for each intersection. Lets us know how to make the line connections
        let mut int_start = Vec::new();
        let mut count = 1;
        // Counter for intersection header  number of points
        let mut num_int_pts = 0;
        // Used in writing which nodes belong to which two fractures in finishWritingOutput()
        let mut intersecting_fractures = Vec::new();
        // Make new intersection file in intersections folder
        let file = format!("{}/intersections_{}.inp", intersection_folder, i + 1);
        let mut fract_int_file = File::create(file).unwrap();
        // Buffer first line to leave space to write header
        fract_int_file
            .write_all(
                "                                                               \n".as_bytes(),
            )
            .unwrap();
        // Go through each final fracture's intersections and write to output
        let size = accepted_poly[final_fractures[i]].intersection_index.len();

        if size > 0 || !keep_isolated_fractures {
            for j in 0..size {
                // tempTripPts holds rotated triple points for an intersection. Triple pts must be rotated 3 different
                // ways so we cannot change the original data
                let mut temp_trip_pts = Vec::new();
                // Used to measure current length before rotation (used to calculate number of points
                // to discretize) based on original intersection. This fixes any precision errors we
                // when calculating length after rotation, both rotations for the same intersection will
                // always have the same step size and same number of discretized points.
                let mut cur_length;
                let poly_int_idx = accepted_poly[final_fractures[i]].intersection_index[j];
                // Similarly to above, the intersection must be rotated two different ways,
                // one for each intersecting poly. We can't change the original data so we must use temp data
                let temp_intersection = poly_and_intersection_rotation_to_xy(
                    &int_pts[poly_int_idx],
                    &mut accepted_poly[final_fractures[i]],
                    triple_points,
                    &mut temp_trip_pts,
                    eps,
                );
                // poly and intersection now rotated
                let triple_pts_size = temp_trip_pts.len();
                // fracture 1 is i
                // fracture 2 is the other intersecting fracture
                let fract2;

                if -int_pts[poly_int_idx].fract1 == i as isize + 1 {
                    fract2 = -int_pts[poly_int_idx].fract2;
                    intersecting_fractures.push(fract2 as usize);
                } else {
                    fract2 = -int_pts[poly_int_idx].fract1;
                    intersecting_fractures.push(fract2 as usize);
                }

                // If triple points exist on intersection, discretize from endpoint to closest triple point,
                // from triple to next triple point, and finally to other end point
                if triple_pts_size != 0 {
                    // Keep track of number of triple points which will be in the
                    // DFN (this is after isolated fracture removal)
                    // NOTE: This will need to be divided by six to get correct value.
                    // Division by six was determined through testing.
                    pstats.triple_node_count += triple_pts_size;
                    // Order the triple points by distances to know to discretize from point to next closest point
                    let mut distances = Vec::with_capacity(triple_pts_size);
                    let mut pt1 = temp_intersection.p1;
                    temp_point1 = int_pts[poly_int_idx].p1;

                    // Create array of distances first end point to triple points
                    for k in 0..triple_pts_size {
                        //loop through triple points on  intersection i
                        distances[k] = distance(&pt1, &temp_trip_pts[k]); //create array of distances
                    }

                    // Order the indices of the distances array shortest to largest distance
                    // this lets us know which point to discritize to next
                    let s = sorted_index(&distances);
                    // Discretize from end point1 to first triple pt
                    // pt1 already = enpoint1
                    let mut pt2 = temp_trip_pts[s[0]];
                    temp_point2 = triple_points[int_pts[poly_int_idx].triple_points_idx[s[0]]];
                    cur_length = distance(&temp_point1, &temp_point2);
                    let mut points = discretize_line_of_intersection(
                        h,
                        visualization_mode,
                        &pt1,
                        &pt2,
                        cur_length,
                    );
                    // Write points to file
                    num_int_pts += points.len();
                    write_points(&mut fract_int_file, &points, 0, &mut count);

                    // If one trip pt, set up points to discretize from only triple pt to other end point
                    if triple_pts_size == 1 {
                        pt1 = pt2;
                        pt2 = temp_intersection.p2;
                        temp_point1 = temp_point2;
                        temp_point2 = int_pts[poly_int_idx].p2;
                    } else {
                        // More than 1 triple point
                        for jj in 0..triple_pts_size - 1 {
                            pt1[0] = temp_trip_pts[s[jj]].x;
                            pt1[1] = temp_trip_pts[s[jj]].y;
                            pt1[2] = temp_trip_pts[s[jj]].z;
                            pt2[0] = temp_trip_pts[s[jj + 1]].x;
                            pt2[1] = temp_trip_pts[s[jj + 1]].y;
                            pt2[2] = temp_trip_pts[s[jj + 1]].z;
                            temp_point1 =
                                triple_points[int_pts[poly_int_idx].triple_points_idx[s[jj]]];
                            temp_point2 =
                                triple_points[int_pts[poly_int_idx].triple_points_idx[s[jj + 1]]];
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
                            write_points(&mut fract_int_file, &points, 1, &mut count);
                        }

                        // Set up points to go from last triple point to last endpoint
                        pt1 = pt2;
                        pt2 = temp_intersection.p2;
                        temp_point1 = temp_point2;
                        temp_point2 = int_pts[poly_int_idx].p2;
                    }

                    cur_length = distance(&temp_point1, &temp_point2);
                    points = discretize_line_of_intersection(
                        h,
                        visualization_mode,
                        &pt1,
                        &pt2,
                        cur_length,
                    );
                    num_int_pts += points.len() - 1;
                    write_points(&mut fract_int_file, &points, 1, &mut count);
                } else {
                    // No triple intersection points on intersection line
                    let pt1 = temp_intersection.p1;
                    let pt2 = temp_intersection.p2;
                    temp_point1 = int_pts[poly_int_idx].p1;
                    temp_point2 = int_pts[poly_int_idx].p2;
                    cur_length = distance(&temp_point1, &temp_point2);
                    let points = discretize_line_of_intersection(
                        h,
                        visualization_mode,
                        &pt1,
                        &pt2,
                        cur_length,
                    );
                    num_int_pts += points.len();
                    write_points(&mut fract_int_file, &points, 0, &mut count);
                }

                int_start.push(count);
            }
        } else {
            // let mut tempTripPts = Vec::new();
            // let temp_intersection = polyAndIntersection_RotationToXY(
            //     &intPts[0],
            //     &mut acceptedPoly[finalFractures[i]],
            //     triplePoints,
            //     &mut tempTripPts,
            // );
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
        );
        intersecting_fractures.clear();
        int_start.clear();
    }

    //Divide by 6 to remove the duplicate counts
    pstats.triple_node_count /= 6;
}

/// Parses and writes all poly_x.inp files containing polygon (fracture) vertice and connectivity data
///
/// # Arguments
///
/// * `final_fractures` - Vector array of indices of fractures left after isolated fracture removal
/// * `accepted_poly` - Vector array of all accetped fractures
/// * `output` - Path to output folder
fn write_polys(final_fractures: &[usize], accepted_poly: &[Poly], output: &str) {
    info!("Writing Polygon Files");
    let poly_count = final_fractures.len();
    let poly_output_file = format!("{}/polygons.dat", output);
    let mut poly_output = File::create(poly_output_file).unwrap();
    poly_output
        .write_all(format!("nPolygons: {}\n", poly_count).as_bytes())
        .unwrap();

    for j in 0..poly_count {
        // Write vertices
        let number_of_nodes = accepted_poly[final_fractures[j]].number_of_nodes;
        //polyOutput << acceptedPoly[finalFractures[j]].familyNum << " ";
        poly_output
            .write_all(format!("{} ", number_of_nodes).as_bytes())
            .unwrap();

        for i in 0..number_of_nodes {
            let idx = i as usize * 3;
            poly_output
                .write_all(
                    format!(
                        "{{{:.12}, {:.12}, {:.12}}} ",
                        accepted_poly[final_fractures[j]].vertices[idx],
                        accepted_poly[final_fractures[j]].vertices[idx + 1],
                        accepted_poly[final_fractures[j]].vertices[idx + 2],
                    )
                    .as_bytes(),
                )
                .unwrap();
        }

        poly_output.write_all("\n".as_bytes()).unwrap();
    }

    info!("Writing Polygon Files Complete");
}

/// Parses and writes all poly_x.inp files containing polygon (fracture) vertice and connectivity data
///
/// # Arguments
///
/// * `final_fractures` - Vector array of indices of fractures left after isolated fracture removal
/// * `accepted_poly` - Vector array of all accetped fractures
/// * `output` - Path to output folder
fn write_polys_inp(final_fractures: &[usize], accepted_poly: &[Poly], output: &str) {
    info!("Writing poly inp files");
    let poly_count = final_fractures.len();

    for j in 0..poly_count {
        let poly_output_file = format!("{}/poly_{}.inp", output, j + 1);
        let mut poly_output = File::create(poly_output_file).unwrap();
        // Write header
        poly_output
            .write_all(
                format!(
                    "{} {} 0 0 0\n",
                    accepted_poly[final_fractures[j]].number_of_nodes,
                    accepted_poly[final_fractures[j]].number_of_nodes - 1,
                )
                .as_bytes(),
            )
            .unwrap();
        // Write vertices
        let number_of_nodes = accepted_poly[final_fractures[j]].number_of_nodes;

        for i in 0..number_of_nodes {
            let idx = i as usize * 3;
            poly_output
                .write_all(
                    format!(
                        "{:.12} {:.12} {:.12} {:.12}\n",
                        i + 1,
                        accepted_poly[final_fractures[j]].vertices[idx],
                        accepted_poly[final_fractures[j]].vertices[idx + 1],
                        accepted_poly[final_fractures[j]].vertices[idx + 2]
                    )
                    .as_bytes(),
                )
                .unwrap();
        }

        // Write line connectivity
        for i in 1..number_of_nodes {
            poly_output
                .write_all(format!("{} {} line {} {}\n", i, j + 1, i, i + 1).as_bytes())
                .unwrap();
        }
    }
}

// Writes params.txt
//
// # Arguments
//
/// * `h` - Minimum feature size
/// * `visualization_mode` - If false, creates a fine mesh, according to h parameter.
///     If true, produce only first round of triangulations. In this case no modeling of flow and transport is possible.
/// * `domain_size` - Vector3<f64> of domain size
/// * `final_fractures` - Vector array of indices of fractures left after isolated fracture removal
/// * `pstats` - Stats structure, running program statistics
/// * `output` - Path to output folder
fn write_params_file(
    h: f64,
    visualization_mode: bool,
    domain_size: &Vector3<f64>,
    final_fractures: &[usize],
    pstats: &Stats,
    output: &str,
) {
    let params_output_file = format!("{}/../params.txt", output);
    let mut params = File::create(&params_output_file).unwrap();
    info!("Writing {}", &params_output_file);
    params
        .write_all(format!("{}\n", final_fractures.len()).as_bytes())
        .unwrap();
    params.write_all(format!("{}\n", h).as_bytes()).unwrap();
    params
        .write_all(format!("{}\n", visualization_mode as u8).as_bytes())
        .unwrap(); // Production mode
    params
        .write_all(
            format!(
                "{}\n",
                pstats.intersection_node_count / 2 - pstats.triple_node_count
            )
            .as_bytes(),
        )
        .unwrap();
    params
        .write_all(format!("{}\n", domain_size[0]).as_bytes())
        .unwrap();
    params
        .write_all(format!("{}\n", domain_size[1]).as_bytes())
        .unwrap();
    params
        .write_all(format!("{}\n", domain_size[2]).as_bytes())
        .unwrap();
}

/// Writes radii.dat (Radii Data)
///
/// # Arguments
///
/// * `final_fractures` - Vector array of indices of fractures left after isolated fracture removal
/// * `accepted_poly` - Vector array of all accetped fractures
/// * `output` - Path to output folder
fn write_radii_file(final_fractures: &[usize], accepted_poly: &[Poly], output: &str) {
    info!("Writing Radii File (radii.dat)");
    let file = format!("{}/radii.dat", output);
    let mut radii = File::create(file).unwrap();
    radii.write_all("Format: xRadius yRadius Family# Removed (-2 = userPolygon, -1 = userRectangle, 0 = userEllipse, > 0 is family in order of famProb)\n".as_bytes()).unwrap();
    let final_fract_limit = final_fractures.len() - 1;
    let size = accepted_poly.len();
    let mut cur_final_idx = 0;

    for (i, poly) in accepted_poly.iter().enumerate().take(size) {
        radii
            .write_all(
                format!(
                    "{:.8} {:.8} {:.8}",
                    poly.xradius,
                    poly.yradius,
                    poly.family_num + 1
                )
                .as_bytes(),
            )
            .unwrap();

        // If poly is not in finalFractures list
        // mark that is was removed
        if i != final_fractures[cur_final_idx] {
            radii.write_all(" R\n".as_bytes()).unwrap();
        } else {
            if cur_final_idx < final_fract_limit {
                cur_final_idx += 1;
            }

            radii.write_all("\n".as_bytes()).unwrap();
        }
    }
}

/// Wrtes translations file (translations.dat)
///
/// # Arguments
///
/// * `final_fractures` - Vector array of indices of fractures left after isolated fracture removal
/// * `accepted_poly` - Vector array of all accetped fractures
/// * `output` - Path to output folder
fn write_fracture_translations(final_fractures: &[usize], accepted_poly: &[Poly], output: &str) {
    info!("Writing Fracture Translations File (translations.dat)");
    let file_path = format!("{}/translations.dat", output);
    let mut file = File::create(file_path).unwrap();
    file.write_all(
        "Format: x y z  (R = removed from domain due to fracture isolation)\n".as_bytes(),
    )
    .unwrap();
    let final_fract_limit = final_fractures.len() - 1;
    let size = accepted_poly.len();
    let mut cur_final_idx = 0;

    for (i, poly) in accepted_poly.iter().enumerate().take(size) {
        file.write_all(
            format!(
                "{:.10} {:.10} {:.10}",
                poly.translation[0], poly.translation[1], poly.translation[2]
            )
            .as_bytes(),
        )
        .unwrap();

        // If poly is not in finalFractures list
        // mark that is was removed
        if i != final_fractures[cur_final_idx] {
            file.write_all(" R\n".as_bytes()).unwrap();
        } else {
            if cur_final_idx < final_fract_limit {
                cur_final_idx += 1;
            }

            file.write_all("\n".as_bytes()).unwrap();
        }
    }
}

/// Deprecated Function
/// Writes final radii file (after isoloated fractures have been removed)
///
/// # Arguments
///
/// * `final_fractures` - Vector array of indices of fractures left after isolated fracture removal
/// * `accepted_poly` - Vector array of all accetped fractures
/// * `output` - Path to output folder
fn write_final_poly_radii(final_fractures: &[usize], accepted_poly: &[Poly], output: &str) {
    let file = format!("{}/radii_Final.dat", output);
    let mut radii_final = File::create(file).unwrap();
    radii_final
        .write_all("Fracture Radii List After Isolated Fracture and Cluster Removal\n".as_bytes())
        .unwrap();
    radii_final.write_all("Format: xRadius yRadius Family# (-2 = userPolygon, -1 = userRectangle, 0 = userEllipse, > 0 is family in order of famProb)\n".as_bytes()).unwrap();
    let size = final_fractures.len();

    for i in 0..size {
        radii_final
            .write_all(
                format!(
                    "{} {} {}\n",
                    accepted_poly[final_fractures[i]].xradius,
                    accepted_poly[final_fractures[i]].yradius,
                    accepted_poly[final_fractures[i]].family_num + 1
                )
                .as_bytes(),
            )
            .unwrap();
    }
}

/// Deprecated Function
/// Writes final radii file (after isoloated fractures have been removed)
///
/// # Arguments
///
/// * `final_fractures` - Vector array of indices of fractures left after isolated fracture removal
/// * `accepted_poly` - Vector array of all accetped fractures
/// * `output` - Path to output folder
fn write_final_poly_area(final_fractures: &[usize], accepted_poly: &[Poly], output: &str) {
    let file = format!("{}/surface_area_Final.dat", output);
    let mut area_final = File::create(file).unwrap();
    area_final
        .write_all("Fracture Surface Area After Isolated Fracture and Cluster Removal\n".as_bytes())
        .unwrap();
    let size = final_fractures.len();

    for i in 0..size {
        area_final
            .write_all(format!("{}\n", accepted_poly[final_fractures[i]].area).as_bytes())
            .unwrap();
    }
}

/// Writes radii file (radii_AllAccepted_Fam_#.dat) for all accepted
/// fractures BEFORE isolated fracture removal (one file per family)
///
/// # Arguments
///
/// * `family_num` - Family number for which radii file will be written for
///     -2 - User Rectangles, -1 - User Ellipses, Family# >= 0 - Family in order of 'FamProb' in input file
/// * `accepted_poly` - Vector array of all accetped fractures
/// * `output` - Path to output folder
fn write_all_accepted_radii_of_family(family_num: isize, accepted_poly: &[Poly], output: &str) {
    let file_name = format!("{}/radii_AllAccepted_Fam_{}.dat", output, family_num + 1);
    let mut file = File::open(file_name).unwrap();
    file.write_all(
        format!(
            "Fracture Radii List Before Isolated Fracture and Cluster Removal (Family {})\n",
            family_num + 1
        )
        .as_bytes(),
    )
    .unwrap();
    file.write_all("Format: xRadius yRadius Distrubution# (-2 = userPolygon, -1 = userRectangle, 0 = userEllipse, > 0 is family in order of famProb)\n".to_string().as_bytes()).unwrap();
    let size = accepted_poly.len();

    for poly in accepted_poly.iter().take(size) {
        if poly.family_num == family_num {
            file.write_all(
                format!(
                    "{} {} {}\n",
                    poly.xradius,
                    poly.yradius,
                    poly.family_num + 1
                )
                .as_bytes(),
            )
            .unwrap();
        }
    }
}

/// Writes radii file (radii_AllAccepted_Fam_#.dat) for all accepted
/// fractures AFTER isolated fracture removal (one file per family)
///
/// # Arguments
///
/// * `final_fractures` - Vector array of indices of fractures left after isolated fracture removal
/// * `family_num` - Family number for which radii file will be written for
///    -2 - User Rectangles, -1 - User Ellipses, Family# >= 0 - Family in order of 'FamProb' in input file
/// * `accepted_poly` - Vector array of all accetped fractures
/// * `output` - Path to output folder
fn write_final_radii_of_family(
    final_fractures: &[usize],
    family_num: isize,
    accepted_poly: &[Poly],
    output: &str,
) {
    let file_name = format!("{}/radii_Final_Fam_{}.dat", output, family_num + 1);
    let mut file = File::open(file_name).unwrap();
    file.write_all(
        format!(
            "Fracture Radii List After Isolated Fracture and Cluster Removal (Family {})\n",
            family_num + 1
        )
        .as_bytes(),
    )
    .unwrap();
    file.write_all("Format: xRadius yRadius Distrubution# (-2 = userPolygon, -1 = userRectangle, 0 = userEllipse, > 0 is family in order of famProb)\n".to_string().as_bytes()).unwrap();
    let size = final_fractures.len();

    for i in 0..size {
        if accepted_poly[final_fractures[i]].family_num == family_num {
            file.write_all(
                format!(
                    "{} {} {}\n",
                    accepted_poly[final_fractures[i]].xradius,
                    accepted_poly[final_fractures[i]].yradius,
                    accepted_poly[final_fractures[i]].family_num + 1
                )
                .as_bytes(),
            )
            .unwrap();
        }
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
/// * `output` - Path to output folder
fn write_triple_pts(
    triple_points: &[Point3<f64>],
    final_fractures: &[usize],
    accepted_poly: &[Poly],
    int_pts: &[IntersectionPoints],
    output: &str,
) {
    let file_name = format!("{}/triple_points.dat", output);
    let mut file = File::create(file_name).unwrap();
    // Save triple point indices from final list of fractures to temp array
    // There will be duplicates  which need to be removed
    // before writing to file
    let mut triple_pts_list = Vec::new();
    let mut size = final_fractures.len();

    for i in 0..size {
        let intersect_count = accepted_poly[final_fractures[i]].intersection_index.len();

        for j in 0..intersect_count {
            let triple_size = int_pts[accepted_poly[final_fractures[i]].intersection_index[j]]
                .triple_points_idx
                .len();

            for k in 0..triple_size {
                let triple_pt_idx = int_pts
                    [accepted_poly[final_fractures[i]].intersection_index[j]]
                    .triple_points_idx[k];
                triple_pts_list.push(triple_pt_idx);
            }
        }
    }

    size = triple_pts_list.len();

    if size > 0 {
        //remove duplicates
        triple_pts_list.sort();
        let mut final_pts = Vec::new();
        let mut prev_pt = triple_pts_list[0];
        final_pts.push(prev_pt);

        for triple_pts in triple_pts_list.iter().take(size).skip(1) {
            let cur_pt = triple_pts;

            if cur_pt != &prev_pt {
                final_pts.push(*cur_pt);
            }

            prev_pt = *cur_pt;
        }

        triple_pts_list.clear();
        size = final_pts.len();

        for i in 0..size {
            file.write_all(
                format!(
                    "{:.17} {:.17} {:.17}\n",
                    triple_points[final_pts[i]].x,
                    triple_points[final_pts[i]].y,
                    triple_points[final_pts[i]].z
                )
                .as_bytes(),
            )
            .unwrap();
        }
    }
}

/// Write rejections.dat, rejection statistics
///
/// # Arguments
///
/// * `pstats` - Stats structure of program statistics
/// * `output` - Path to output folder
fn write_rejection_stats(pstats: &Stats, output: &str) {
    info!("Writing Rejection Statistics File (rejections.dat)");
    let _ = std::fs::create_dir_all(output);
    let file_name = format!("{}/rejections.dat", output);
    let mut file = File::create(file_name).unwrap();
    file.write_all(
        format!(
            "Short Intersection: {}\n",
            pstats.rejection_reasons.short_intersection
        )
        .as_bytes(),
    )
    .unwrap();
    file.write_all(
        format!(
            "Close to Node: {}\n",
            pstats.rejection_reasons.close_to_node
        )
        .as_bytes(),
    )
    .unwrap();
    file.write_all(
        format!(
            "Close to Edge: {}\n",
            pstats.rejection_reasons.close_to_edge
        )
        .as_bytes(),
    )
    .unwrap();
    file.write_all(
        format!(
            "Vertex Close to Edge: {}\n",
            pstats.rejection_reasons.close_point_to_edge
        )
        .as_bytes(),
    )
    .unwrap();
    file.write_all(format!("Outside of Domain: {}\n", pstats.rejection_reasons.outside).as_bytes())
        .unwrap();
    file.write_all(
        format!("Triple Intersection: {}\n", pstats.rejection_reasons.triple).as_bytes(),
    )
    .unwrap();
    file.write_all(
        format!(
            "Intersections Too Close: {}\n",
            pstats.rejection_reasons.inter_close_to_inter
        )
        .as_bytes(),
    )
    .unwrap();
}

/// Write rejections.dat, rejection statistics
///
/// # Arguments
///
/// * `pstats` - Stats structure of program statistics
/// * `output` - Path to output folder
fn write_user_rejected_fracture_information(pstats: &Stats, output: &str) {
    if !pstats.rejected_user_fracture.is_empty() {
        info!("Writing User Fracture Rejection File (userFractureRejections.dat)");
        let file_name = format!("{}/userFractureRejections.dat", output);
        let mut file = File::open(file_name).unwrap();
        file.write_all("Fracture id,User Fracture Type\n".as_bytes())
            .unwrap();

        for i in 0..pstats.rejected_user_fracture.len() {
            file.write_all(
                format!(
                    "{},{}\n",
                    pstats.rejected_user_fracture[i].id,
                    pstats.rejected_user_fracture[i].user_fracture_type
                )
                .as_bytes(),
            )
            .unwrap();
        }
    }
}

/// Writes families.dat, Shape families definition file
///
/// # Arguments
///
/// * `user_ellipses_on_off` - If true, user defined ellipses are used
/// * `user_rectangles_on_off` - If true, user defined rectangles are used
/// * `user_polygon_by_coord` - If true, user defined polygons are used
/// * `stop_condition` - Stop condition for fracture generation (0 - nPoly, 1 - P32)
/// * `orientation_option` - Orientation option for fracture generation
/// * `n_fam_ell` - Number of ellipse families
/// * `layers` - Vector array of layers
/// * `regions` - Vector array of regions
/// * `fam_prob_original` - Vector array of family probabilities
/// * `frac_families` - Vector array of all fracture shape families
/// * `output` - Path to output folder
#[allow(clippy::too_many_arguments)]
fn write_frac_fams(
    user_ellipses_on_off: bool,
    user_rectangles_on_off: bool,
    user_polygon_by_coord: bool,
    stop_condition: u8,
    n_fam_ell: usize,
    layers: &[f64],
    regions: &[f64],
    frac_families: &FractureFamilyOption,
    output: &str,
) {
    let rad_to_deg = 180. / std::f64::consts::PI;
    info!("Writing Family Definitions File (families.dat)");
    let file_name = format!("{}/families.dat", output);
    let mut file = File::create(file_name).unwrap();

    //TODO: add stub code in families.dat for userDefined fractures, IF there are user defined fractures

    if user_ellipses_on_off {
        file.write_all("UserDefined Ellipse Family: 0\n\n".as_bytes())
            .unwrap();
    }

    if user_rectangles_on_off {
        file.write_all("UserDefined Rectangle Family: -1\n\n".as_bytes())
            .unwrap();
    }

    if user_polygon_by_coord {
        file.write_all("UserDefined Polygon Family: -2\n\n".as_bytes())
            .unwrap();
    }

    for (i, shape) in frac_families.families.iter().enumerate() {
        //name(rect or ell) and number of family
        file.write_all(
            format!(
                "{} Family: {}\n",
                shape.shape,
                get_family_number(n_fam_ell, i as isize, shape.shape)
            )
            .as_bytes(),
        )
        .unwrap();
        file.write_all(format!("Global Family: {}\n", i + 1).as_bytes())
            .unwrap();

        // Print vertice number
        match shape.shape {
            Shape::Ellipse(n) => {
                file.write_all(format!("Number of Vertices: {}\n", n).as_bytes())
                    .unwrap();
            }
            Shape::Rectangle => {
                file.write_all("Number of Vertices: 4\n".as_bytes())
                    .unwrap();
            }
        }

        // aspect ratio
        file.write_all(format!("Aspect Ratio: {}\n", shape.aspect_ratio).as_bytes())
            .unwrap();

        // p32 target
        if stop_condition == 1 {
            file.write_all(
                format!("P32 (Fracture Intensity) Target: {}\n", shape.p32_target).as_bytes(),
            )
            .unwrap();
        }

        // beta distribution, rotation around normal vector
        if !shape.beta_distribution {
            file.write_all(
                "Beta Distribution (Rotation Around Normal Vector): Uniform\n".as_bytes(),
            )
            .unwrap();
        } else {
            file.write_all(
                format!("Beta (Rotation Around Normal Vector)-rad: {}\n", shape.beta).as_bytes(),
            )
            .unwrap();
            file.write_all(
                format!(
                    "Beta (Rotation Around Normal Vector)-deg: {}\n",
                    shape.beta * rad_to_deg
                )
                .as_bytes(),
            )
            .unwrap();
        }

        file.write_all(shape.orientation.to_string().as_bytes())
            .unwrap();

        // Print layer family belongs to
        if shape.layer == 0 {
            file.write_all("Layer: Entire domain\n".as_bytes()).unwrap();
        } else {
            let idx = (shape.layer - 1) * 2;
            file.write_all(format!("Layer Number: {}\n", shape.layer).as_bytes())
                .unwrap();
            file.write_all(format!("Layer: {{{}, {}}}\n", layers[idx], layers[idx + 1]).as_bytes())
                .unwrap();
        }

        // Print layer family belongs to
        if shape.region == 0 {
            file.write_all("Region: Entire domain\n".as_bytes())
                .unwrap();
        } else {
            let idx = (shape.region - 1) * 6;
            file.write_all(format!("Region Number: {}\n", shape.region).as_bytes())
                .unwrap();
            file.write_all(
                format!(
                    "Region: {{{}, {}, {}, {}, {}, {}}}\n",
                    regions[idx],
                    regions[idx + 1],
                    regions[idx + 2],
                    regions[idx + 3],
                    regions[idx + 4],
                    regions[idx + 5]
                )
                .as_bytes(),
            )
            .unwrap();
        }

        // Print distribution data
        file.write_all(shape.radius.to_string().as_bytes()).unwrap();

        file.write_all(
            format!(
                "Family Insertion Probability: {}\n\n",
                frac_families.original_probabilities[i]
            )
            .as_bytes(),
        )
        .unwrap();
    }
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
) {
    info!("Writing Connectivity Data (connectivity.dat)");
    let file_name = format!("{}/connectivity.dat", output);
    let mut file = File::create(file_name).unwrap();

    // Format for connectivity: Line number = poly number, integers on line are intersecting fractures
    //  eg:
    // 2 3
    // 1 3
    // 1 2
    // Fracture 1 intersects frac 2 and 3, fract 2 intersects 1 and 3 etc

    for i in 0..final_fractures.len() {
        let size = accepted_poly[final_fractures[i]].intersection_index.len();

        for j in 0..size {
            let int_idx = accepted_poly[final_fractures[i]].intersection_index[j];

            // Don't write fractures own fracture ID

            // NOTE: At this stage in the program, the fracture ID's have been changed to
            // negative. See 'adjustIntFractIds()' for more details.
            // Use '-' to make them positive again.
            if -int_pts[int_idx].fract1 == i as isize + 1 {
                file.write_all(format!("{} ", -int_pts[int_idx].fract2).as_bytes())
                    .unwrap();
            } else {
                file.write_all(format!("{} ", -int_pts[int_idx].fract1).as_bytes())
                    .unwrap();
            }
        }

        file.write_all("\n".as_bytes()).unwrap();
    }
}

/// Writes poly_info.dat
/// Writes fracture rotation data. Also includes shape families each fracture belongs to.
///
/// # Arguments
///
/// * `eps` - Epsilon value for floating point comparisons
/// * `domain_size` - domain size in x, y, z
/// * `accepted_poly` - Array off all polygons in DFN before isolated fracture removal
/// * `final_fractures` - Array of indices to polys in 'accepted_poly' which are left after isolated fracture removal
/// * `frac_families` - Array of all fracture shape families
/// * `output` - Path to output folder
fn write_rotation_data(
    eps: f64,
    domain_size: &Vector3<f64>,
    accepted_poly: &[Poly],
    final_fractures: &[usize],
    frac_families: &FractureFamilyOption,
    output: &str,
) {
    let file_output_file = format!("{}/../poly_info.dat", output);
    let mut file = File::create(file_output_file).unwrap();
    info!("Writing Rotation Data File (poly_info.dat)");
    let mut max_domain_size = domain_size[0];

    if max_domain_size < domain_size[1] {
        max_domain_size = domain_size[1];
    }

    if max_domain_size < domain_size[2] {
        max_domain_size = domain_size[2];
    }

    max_domain_size *= 10.;

    for i in 0..final_fractures.len() {
        // poly's normal is already normalized at this point
        let normal = accepted_poly[final_fractures[i]].normal;
        let e3 = Vector3::new(0., 0., 1.);
        // Rotation angle in radians
        let mut theta = normal.dot(&e3).acos();
        // rad to deg
        theta *= 180.0 / std::f64::consts::PI;
        // Rotation into xy plane
        let mut v = e3.cross(&normal);

        if !(v[0].abs() < eps && v[1].abs() < eps && v[2].abs() < eps) {
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
        let fam_num = if accepted_poly[final_fractures[i]].family_num == -1 {
            frac_families.families.len() as isize + 1
        } else if accepted_poly[final_fractures[i]].family_num == -2 {
            frac_families.families.len() as isize + 2
        } else {
            accepted_poly[final_fractures[i]].family_num + 1
        };

        // Format: fracture#, x0, y0, z0, x1, y1, z1, family#
        file.write_all(
            format!(
                "{} {} {:.15} {:.15} {:.15} {:.15} {:.15} {:.15} {:.15}\n",
                i + 1,
                fam_num,
                theta,
                x0,
                y0,
                z0,
                x1,
                y1,
                z1
            )
            .as_bytes(),
        )
        .unwrap();
    }
}

/// Writes normal_vectors.dat
///
/// Writes fracture rotation data. Also includes shape families each fracture belongs to.
///
/// # Arguments
///
/// * `accepted_poly` - Array off all polygons in DFN before isolated fracture removal
/// * `final_fractures` - Array of indices to polys in 'accepted_poly' which are left after isolated fracture removal
/// * `output` - Path to output folder
fn write_normal_vectors(accepted_poly: &[Poly], final_fractures: &[usize], output: &str) {
    let file_output_file = format!("{}/normal_vectors.dat", output);
    let mut file = File::create(file_output_file).unwrap();
    info!("Writing Normal Vectors into File (normal_vectors.dat)");

    for i in 0..final_fractures.len() {
        // poly's normal is already normalized at this point
        // Format: nx, ny, nz
        //std::cout <<  std::setprecision(15) << acceptedPoly[finalFractures[i]].normal[0] << " "
        //<< acceptedPoly[finalFractures[i]].normal[1] << " " << acceptedPoly[finalFractures[i]].normal[2]  << "\n";
        file.write_all(
            format!(
                " {:.15} {:.15} {:.15}\n",
                accepted_poly[final_fractures[i]].normal[0],
                accepted_poly[final_fractures[i]].normal[1],
                accepted_poly[final_fractures[i]].normal[2]
            )
            .as_bytes(),
        )
        .unwrap();
    }
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
fn write_rejects_per_attempt(pstats: &Stats, output: &str) {
    let file_output_file = format!("{}/rejectsPerAttempt.dat", output);
    let mut file = File::create(file_output_file).unwrap();
    info!("Writing Rotation Data File (rejectsPerAttempt.dat)");

    for i in 0..pstats.rejects_per_attempt.len() {
        file.write_all(format!("{}\n", pstats.rejects_per_attempt[i]).as_bytes())
            .unwrap();
    }
}

/// Writes graph data files to intersections_list.dat and fracture_info.dat
///
/// # Arguments
///
/// * `eps` - Epsilon value for floating point comparisons
/// * `domain_size` - domain size in x, y, z
/// * `final_fractures` - Vector array of indices of fractures left after isolated fracture removal
/// * `accepted_poly` - Vector array of all accetped fractures (before isolated fracture removal)
/// * `int_pts` - Vector array of all intersections
/// * `output` - Path to output folder
fn write_graph_data(
    eps: f64,
    domain_size: &Vector3<f64>,
    final_fractures: &[usize],
    accepted_poly: &[Poly],
    int_pts: &[IntersectionPoints],
    output: &str,
) {
    let domain_x = domain_size[0] * 0.5;
    let domain_y = domain_size[1] * 0.5;
    let domain_z = domain_size[2] * 0.5;
    // Keeps track of current un-rotated points we are working with
    info!("Writing Graph Data Files");
    //adjustIntFractIDs(finalFractures,acceptedPoly, intPts);
    // Make new intersection file in intersections folder
    let file = format!("{}/intersection_list.dat", output);
    let mut int_file = File::create(file).unwrap();
    int_file
        .write_all("f1 f2 x y z length\n".as_bytes())
        .unwrap();
    // Make new intersection file in intersections folder
    let file2 = format!("{}/fracture_info.dat", output);
    let mut fract_file = File::create(file2).unwrap();
    fract_file
        .write_all("num_connections perm aperture\n".as_bytes())
        .unwrap();

    for i in 0..final_fractures.len() {
        let mut num_conn = 0;
        // Go through each final fracture's intersections and write to output
        let size = accepted_poly[final_fractures[i]].intersection_index.len();

        for j in 0..size {
            // Used to measure current length before rotation (used to calculate number of points
            // to discretize) based on original intersection. This fixes any precision errors we
            // when calculating length after rotation, both rotations for the same intersection will
            // always have the same step size and same number of discretized points.
            let poly_int_idx = accepted_poly[final_fractures[i]].intersection_index[j];
            // fracture 1 is i
            // fracture 2 is the other intersecting fracture
            let fract1;
            let fract2;

            if -int_pts[poly_int_idx].fract1 == i as isize + 1 {
                fract1 = -int_pts[poly_int_idx].fract1;
                fract2 = -int_pts[poly_int_idx].fract2;
            } else {
                fract2 = -int_pts[poly_int_idx].fract1;
                fract1 = -int_pts[poly_int_idx].fract2;
            }

            if fract1 < fract2 {
                write_mid_point(
                    &mut int_file,
                    fract1,
                    fract2,
                    int_pts[poly_int_idx].p1.x,
                    int_pts[poly_int_idx].p1.y,
                    int_pts[poly_int_idx].p1.z,
                    int_pts[poly_int_idx].p2.x,
                    int_pts[poly_int_idx].p2.y,
                    int_pts[poly_int_idx].p2.z,
                );
                num_conn += 1;
            }
        }

        // Find intersections with domain boundaries
        // bools are so that lines of intersection are only found and written to file once
        // once an intersection with a domain boundary is found the bool = true
        let mut right = false;
        let mut left = false;
        let mut top = false;
        let mut bottom = false;
        let mut front = false;
        let mut back = false;
        let temp = accepted_poly[final_fractures[i]].number_of_nodes;

        // fracture 2 index are based on FEHM format
        // top / Z+ / 1
        // bottom / Z- / 2
        // left / X- / 3
        // front / Y+ / 4
        // right / X- / 5
        // back / Y- / 6
        for k in 0..temp {
            // Check boundary faces
            // Update which boundaries newPoly touches
            let idx = k as usize * 3;

            // if fracture if on right boundary
            if accepted_poly[final_fractures[i]].vertices[idx] >= domain_x - eps {
                for kk in 0..temp {
                    let iidx = kk as usize * 3;

                    if accepted_poly[final_fractures[i]].vertices[iidx] >= domain_x - eps
                        && idx != iidx
                        && !right
                    {
                        write_mid_point(
                            &mut int_file,
                            i as isize + 1,
                            -5,
                            accepted_poly[final_fractures[i]].vertices[idx],
                            accepted_poly[final_fractures[i]].vertices[idx + 1],
                            accepted_poly[final_fractures[i]].vertices[idx + 2],
                            accepted_poly[final_fractures[i]].vertices[iidx],
                            accepted_poly[final_fractures[i]].vertices[iidx + 1],
                            accepted_poly[final_fractures[i]].vertices[iidx + 2],
                        );
                        right = true;
                        num_conn += 1;
                        break;
                    }
                }
            }
            // if fracture if on left boundary
            else if accepted_poly[final_fractures[i]].vertices[idx] <= -domain_x + eps {
                for kk in 0..temp {
                    let iidx = kk as usize * 3;

                    if accepted_poly[final_fractures[i]].vertices[iidx] <= -domain_x + eps
                        && idx != iidx
                        && !left
                    {
                        write_mid_point(
                            &mut int_file,
                            i as isize + 1,
                            -3,
                            accepted_poly[final_fractures[i]].vertices[idx],
                            accepted_poly[final_fractures[i]].vertices[idx + 1],
                            accepted_poly[final_fractures[i]].vertices[idx + 2],
                            accepted_poly[final_fractures[i]].vertices[iidx],
                            accepted_poly[final_fractures[i]].vertices[iidx + 1],
                            accepted_poly[final_fractures[i]].vertices[iidx + 2],
                        );
                        left = true;
                        num_conn += 1;
                        break;
                    }
                }
            }

            // if fracture if on front boundary
            if accepted_poly[final_fractures[i]].vertices[idx + 1] >= domain_y - eps {
                for kk in 0..temp {
                    let iidx = kk as usize * 3;

                    if accepted_poly[final_fractures[i]].vertices[iidx + 1] >= domain_y - eps
                        && idx != iidx
                        && !front
                    {
                        write_mid_point(
                            &mut int_file,
                            i as isize + 1,
                            -4,
                            accepted_poly[final_fractures[i]].vertices[idx],
                            accepted_poly[final_fractures[i]].vertices[idx + 1],
                            accepted_poly[final_fractures[i]].vertices[idx + 2],
                            accepted_poly[final_fractures[i]].vertices[iidx],
                            accepted_poly[final_fractures[i]].vertices[iidx + 1],
                            accepted_poly[final_fractures[i]].vertices[iidx + 2],
                        );
                        front = true;
                        num_conn += 1;
                        break;
                    }
                }
            }
            // if fracture if on back boundary
            else if accepted_poly[final_fractures[i]].vertices[idx + 1] <= -domain_y + eps {
                for kk in 0..temp {
                    let iidx = kk as usize * 3;

                    if accepted_poly[final_fractures[i]].vertices[iidx + 1] <= -domain_y + eps
                        && idx != iidx
                        && !back
                    {
                        write_mid_point(
                            &mut int_file,
                            i as isize + 1,
                            -6,
                            accepted_poly[final_fractures[i]].vertices[idx],
                            accepted_poly[final_fractures[i]].vertices[idx + 1],
                            accepted_poly[final_fractures[i]].vertices[idx + 2],
                            accepted_poly[final_fractures[i]].vertices[iidx],
                            accepted_poly[final_fractures[i]].vertices[iidx + 1],
                            accepted_poly[final_fractures[i]].vertices[iidx + 2],
                        );
                        back = true;
                        num_conn += 1;
                        break;
                    }
                }
            }

            if accepted_poly[final_fractures[i]].vertices[idx + 2] >= domain_z - eps {
                for kk in 0..temp {
                    let iidx = kk as usize * 3;

                    if accepted_poly[final_fractures[i]].vertices[iidx + 2] >= domain_z - eps
                        && idx != iidx
                        && !top
                    {
                        write_mid_point(
                            &mut int_file,
                            i as isize + 1,
                            -1,
                            accepted_poly[final_fractures[i]].vertices[idx],
                            accepted_poly[final_fractures[i]].vertices[idx + 1],
                            accepted_poly[final_fractures[i]].vertices[idx + 2],
                            accepted_poly[final_fractures[i]].vertices[iidx],
                            accepted_poly[final_fractures[i]].vertices[iidx + 1],
                            accepted_poly[final_fractures[i]].vertices[iidx + 2],
                        );
                        top = true;
                        num_conn += 1;
                        break;
                    }
                }
            } else if accepted_poly[final_fractures[i]].vertices[idx + 2] <= -domain_z + eps {
                for kk in 0..temp {
                    let iidx = kk as usize * 3;

                    if accepted_poly[final_fractures[i]].vertices[iidx + 2] <= -domain_z + eps
                        && idx != iidx
                        && !bottom
                    {
                        write_mid_point(
                            &mut int_file,
                            i as isize + 1,
                            -2,
                            accepted_poly[final_fractures[i]].vertices[idx],
                            accepted_poly[final_fractures[i]].vertices[idx + 1],
                            accepted_poly[final_fractures[i]].vertices[idx + 2],
                            accepted_poly[final_fractures[i]].vertices[iidx],
                            accepted_poly[final_fractures[i]].vertices[iidx + 1],
                            accepted_poly[final_fractures[i]].vertices[iidx + 2],
                        );
                        bottom = true;
                        num_conn += 1;
                        break;
                    }
                }
            }
        }

        // fractFile << num_conn << " " << std::setprecision(10) << acceptedPoly[finalFractures[i]].permeability << " " << acceptedPoly[finalFractures[i]].aperture << "\n";
        fract_file
            .write_all(format!("{} 0 0\n", num_conn).as_bytes())
            .unwrap();
    }
}

/// Writes mid point and length of line defined by x1,y1,z1 and x2,y2,z2 into file fp
///
/// # Arguments
///
/// * `fp` - File to write information into
/// * `fract1` - Fracture 1
/// * `fract2` - Fracture 2
/// * `x1` - x coordinate of first endpoint
/// * `y1` - y coordinate of first endpoint
/// * `z1` - z coordinate of first endpoint
/// * `x2` - x coordinate of second endpoint
/// * `y2` - y coordinate of second endpoint
/// * `z2` - z coordinate of second endpoint
#[allow(clippy::too_many_arguments)]
fn write_mid_point(
    fp: &mut File,
    fract1: isize,
    fract2: isize,
    x1: f64,
    y1: f64,
    z1: f64,
    x2: f64,
    y2: f64,
    z2: f64,
) {
    // Keeps track of current un-rotated points we are working with
    let temp_point1 = Point3::new(x1, y1, z1);
    let temp_point2 = Point3::new(x2, y2, z2);
    let temp_point3 = Point3::new(
        0.5 * (temp_point1.x + temp_point2.x),
        0.5 * (temp_point1.y + temp_point2.y),
        0.5 * (temp_point1.z + temp_point2.z),
    );
    let cur_length = distance(&temp_point1, &temp_point2);
    fp.write_all(
        format!(
            "{} {} {:.10} {:.10} {:.10} {:.10}\n",
            fract1, fract2, temp_point3.x, temp_point3.y, temp_point3.z, cur_length
        )
        .as_bytes(),
    )
    .unwrap();
}

/// Writes fracture numbers into ASCII files corresponding to which boundary they touch
///
/// # Arguments
///
/// * `final_fractures` - Vector array of indices of fractures left after isolated fracture removal
/// * `accepted_poly` - Vector array of all accepted fractures (before isolated fracture removal)
fn write_boundary_files(final_fractures: &[usize], accepted_poly: &[Poly], output: &str) {
    info!("Writing Boundary Files");
    let left_file_name = format!("{}/left.dat", output);
    let mut left_file = File::create(left_file_name).unwrap();
    let right_file_name = format!("{}/right.dat", output);
    let mut right_file = File::create(right_file_name).unwrap();
    let front_file_name = format!("{}/front.dat", output);
    let mut front_file = File::create(front_file_name).unwrap();
    let back_file_name = format!("{}/back.dat", output);
    let mut back_file = File::create(back_file_name).unwrap();
    let top_file_name = format!("{}/top.dat", output);
    let mut top_file = File::create(top_file_name).unwrap();
    let bottom_file_name = format!("{}/bottom.dat", output);
    let mut bottom_file = File::create(bottom_file_name).unwrap();

    for i in 0..final_fractures.len() {
        // If touching X max
        if accepted_poly[final_fractures[i]].faces[0] {
            right_file
                .write_all(format!("{}\n", i + 1).as_bytes())
                .unwrap();
        }

        // If touching X min
        if accepted_poly[final_fractures[i]].faces[1] {
            left_file
                .write_all(format!("{}\n", i + 1).as_bytes())
                .unwrap();
        }

        // If touching Y max
        if accepted_poly[final_fractures[i]].faces[2] {
            front_file
                .write_all(format!("{}\n", i + 1).as_bytes())
                .unwrap();
        }

        // If touching Y min
        if accepted_poly[final_fractures[i]].faces[3] {
            back_file
                .write_all(format!("{}\n", i + 1).as_bytes())
                .unwrap();
        }

        // If touching Z max
        if accepted_poly[final_fractures[i]].faces[4] {
            top_file
                .write_all(format!("{}\n", i + 1).as_bytes())
                .unwrap();
        }

        // If touching Z min
        if accepted_poly[final_fractures[i]].faces[5] {
            bottom_file
                .write_all(format!("{}\n", i + 1).as_bytes())
                .unwrap();
        }
    }
}
