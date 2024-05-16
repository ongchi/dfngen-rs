use parry3d_f64::na::{distance, Point3, Vector3};

use super::domain::domain_truncation;
use super::insert_shape::print_reject_reason;

use crate::{
    computational_geometry::{create_bounding_box, intersection_checking},
    io::input::Input,
    math_functions::get_area,
    structures::{IntersectionPoints, Poly, RejectedUserFracture, Stats},
};

// ****************************************************************
// ***********  Insert User Ellipses By Coord  ********************
//     Inserts user ellipses using defined coordinates
//     provided by the user (see input file).
//     Intersection checking, FRAM, and rejection/accptance are all contained
//     within this function.
//     Arg 1: Array for all accepted polygons
//     Arg 2: Array for all accepted intersections
//     Arg 3: Program statistics structure
//     Arg 4: Array of all triple intersection points
pub fn insert_user_ell_by_coord(
    input: &mut Input,
    accepted_poly: &mut Vec<Poly>,
    intpts: &mut Vec<IntersectionPoints>,
    pstats: &mut Stats,
    triple_points: &mut Vec<Point3<f64>>,
) {
    println!(
        "{} User Ellipses By Coordinates Defined",
        &input.nEllByCoord
    );

    for i in 0..input.nEllByCoord {
        let mut new_poly = Poly::default();
        let mut rejected_user_fracture = RejectedUserFracture::default();
        new_poly.family_num = -1; // Using -1 for all user specified ellipses
        new_poly.vertices.reserve(input.nEllNodes * 3); // 3 * number of nodes
                                                        // Set number of nodes  - needed for rotations
        new_poly.number_of_nodes = input.nEllNodes as isize;
        let poly_vert_idx = i * 3 * input.nEllNodes; // Each polygon has nEllNodes * 3 vertices

        // Initialize vertices
        for j in 0..input.nEllNodes {
            let v_idx = j * 3;
            new_poly.vertices[v_idx] = input.userEllCoordVertices[poly_vert_idx + v_idx];
            new_poly.vertices[v_idx + 1] = input.userEllCoordVertices[poly_vert_idx + 1 + v_idx];
            new_poly.vertices[v_idx + 2] = input.userEllCoordVertices[poly_vert_idx + 2 + v_idx];
        }

        // Get a normal vector
        // Vector from fist node to node accross middle of polygon
        let mid_pt_idx = 3 * (input.nEllNodes / 2);
        let v1 = Vector3::new(
            new_poly.vertices[mid_pt_idx] - new_poly.vertices[0],
            new_poly.vertices[mid_pt_idx + 1] - new_poly.vertices[1],
            new_poly.vertices[mid_pt_idx + 2] - new_poly.vertices[2],
        );
        // Vector from first node to 2nd node
        let v2 = Vector3::new(
            new_poly.vertices[3] - new_poly.vertices[0],
            new_poly.vertices[4] - new_poly.vertices[1],
            new_poly.vertices[5] - new_poly.vertices[2],
        );
        let x_prod1 = v2.cross(&v1).normalize();
        // Set normal vector
        new_poly.normal = x_prod1;
        // Estimate radius
        new_poly.xradius = 0.5 * v2.magnitude(); // across middle if even number of nodes
        let temp_idx1 = 3 * (mid_pt_idx / 2); // Get idx for node 1/4 around polygon
        let temp_idx2 = 3 * (temp_idx1 + mid_pt_idx); // Get idx for node 3/4 around polygon
                                                      // across middle close to perpendicular to xradius magnitude calculation
        new_poly.yradius = 0.5
            * distance(
                &Point3::from_slice(&new_poly.vertices[temp_idx1..temp_idx1 + 3]),
                &Point3::from_slice(&new_poly.vertices[temp_idx2..temp_idx2 + 3]),
            );
        new_poly.aspect_ratio = new_poly.yradius / new_poly.xradius;
        // Estimate translation (middle of poly)
        // Use midpoint between 1st and and half way around polygon
        // Note: For polygons defined by coordinates, the coordinates
        // themselves provide the translation. We need to estimate the center
        // of the polygon and init. the translation array
        new_poly.translation[0] = 0.5 * (new_poly.vertices[0] + new_poly.vertices[mid_pt_idx]);
        new_poly.translation[1] = 0.5 * (new_poly.vertices[1] + new_poly.vertices[mid_pt_idx + 1]);
        new_poly.translation[2] = 0.5 * (new_poly.vertices[2] + new_poly.vertices[mid_pt_idx + 2]);

        if domain_truncation(input, &mut new_poly, &input.domainSize) {
            // Poly completely outside domain
            new_poly.vertices.clear();
            pstats.rejection_reasons.outside += 1;
            pstats.rejected_poly_count += 1;
            println!("User Ellipse (defined by coordinates) {} was rejected for being outside the defined domain.", i + 1);
            rejected_user_fracture.id = (i + 1) as isize;
            rejected_user_fracture.user_fracture_type = -1;
            pstats.rejected_user_fracture.push(rejected_user_fracture);
            continue; // Go to next poly (go to next iteration of for loop)
        }

        create_bounding_box(&mut new_poly);
        // Line of intersection and FRAM
        let reject_code = intersection_checking(
            input,
            &mut new_poly,
            accepted_poly,
            intpts,
            pstats,
            triple_points,
        );

        if reject_code == 0 {
            // If intersection is ok (FRAM passed all tests)
            if new_poly.truncated {
                pstats.truncated += 1;
            }

            // Incriment counter of accepted polys
            pstats.accepted_poly_count += 1;
            // Calculate poly's area
            new_poly.area = get_area(&new_poly);
            // Add new rejectsPerAttempt counter
            pstats.rejects_per_attempt.push(0);
            println!(
                "User Defined Elliptical Fracture (Defined By Coordinates) {} Acceptedn",
                i + 1
            );
            accepted_poly.push(new_poly); // Save newPoly to accepted polys list
        } else {
            new_poly.vertices.clear(); // Need to delete manually, created with new[]
            pstats.rejects_per_attempt[pstats.accepted_poly_count] += 1;
            pstats.rejected_poly_count += 1;
            println!(
                "Rejected Eser Defined Elliptical Fracture (Defined By Coordinates) {}",
                i + 1
            );
            print_reject_reason(reject_code, &new_poly);
            rejected_user_fracture.id = (i + 1) as isize;
            rejected_user_fracture.user_fracture_type = -1;
            pstats.rejected_user_fracture.push(rejected_user_fracture);
        }
    } // End loop

    input.userEllCoordVertices.clear();
}
