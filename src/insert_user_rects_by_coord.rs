use crate::{
    computational_geometry::{create_bounding_box, intersection_checking},
    domain::domain_truncation,
    insert_shape::print_reject_reason,
    math_functions::get_area,
    read_input::Input,
    structures::{IntPoints, Point, Poly, RejectedUserFracture, Stats},
    vector_functions::{cross_product, magnitude, normalize},
};

// *************************************************************
// ***********  Insert User Rects By Coord  ********************
//     Inserts user rectangles using defined coordinates
//     provided by the user (see input file).
//     Intersection checking, FRAM, and rejection/accptance are all contained
//     within this function.
//     Arg 1: Array for all accepted polygons
//     Arg 2: Array for all accepted intersections
//     Arg 3: Program statistics structure
//     Arg 4: Array of all triple intersection points
pub fn insert_user_rects_by_coord(
    input: &mut Input,
    accepted_poly: &mut Vec<Poly>,
    intpts: &mut Vec<IntPoints>,
    pstats: &mut Stats,
    triple_points: &mut Vec<Point>,
) {
    println!(
        "{} User Rectangles By Coordinates Defined",
        input.nRectByCoord
    );

    for i in 0..input.nRectByCoord {
        let mut new_poly = Poly::default();
        let mut rejected_user_fracture = RejectedUserFracture::default();
        new_poly.family_num = -2; // Using -2 for all user specified rectangles
        new_poly.vertices.reserve(12); // 4 * {x,y,z}
                                       // Set number of nodes  - needed for rotations
        new_poly.number_of_nodes = 4;
        let poly_vert_idx = i * 12; // Each polygon has 4 vertices (12 elements, 4*{x,y,z}))

        // Initialize vertices
        for j in 0..4 {
            let v_idx = j * 3;
            new_poly.vertices[v_idx] = input.userRectCoordVertices[poly_vert_idx + v_idx];
            new_poly.vertices[v_idx + 1] = input.userRectCoordVertices[poly_vert_idx + 1 + v_idx];
            new_poly.vertices[v_idx + 2] = input.userRectCoordVertices[poly_vert_idx + 2 + v_idx];
        }

        // Check that rectangle lays one a single plane:
        // let xProd1 = cross Product vector (1st node to 2nd node) with vector(1st node to 3rd node)
        // and xProd2 = cross product vector (1st node to 3th node) with vector (1st node to 4th node)
        // Then, cross product xProd1 and xProd2, if this produces zero vector, all coords are on the same plane
        // v1 is vector from first vertice to third vertice
        // Vector from fist node to 3rd node (vector through middle of sqare)
        let v1 = [
            new_poly.vertices[6] - new_poly.vertices[0],
            new_poly.vertices[7] - new_poly.vertices[1],
            new_poly.vertices[8] - new_poly.vertices[2],
        ];
        // Vector from first node to 2nd node
        let v2 = [
            new_poly.vertices[3] - new_poly.vertices[0],
            new_poly.vertices[4] - new_poly.vertices[1],
            new_poly.vertices[5] - new_poly.vertices[2],
        ];
        let x_prod1 = cross_product(&v2, &v1);
        // Vector from fist node to 4th node
        let v3 = [
            new_poly.vertices[9] - new_poly.vertices[0],
            new_poly.vertices[10] - new_poly.vertices[1],
            new_poly.vertices[11] - new_poly.vertices[2],
        ];
        // let x_prod2 = crossProduct(&v3, &v1);
        // let x_prod3 = crossProduct(&x_prod1, &x_prod2);
        //will be zero vector if all vertices are on the same plane
        //TODO: Error check below is too sensitive. Adjust it.
        // Error check for points not on the same plane
        //        if (std::abs(magnitude(xProd3[0],xProd3[1],xProd3[2])) > eps) { //points do not lay on the same plane. reject poly else meshing will fail
        // if (!(std::abs(xProd3[0]) < eps && std::abs(xProd3[1]) < eps && std::abs(xProd3[2]) < eps)) {
        //     delete[] newPoly.vertices;
        //     pstats.rejectedPolyCount++;
        //     std::cout << "\nUser Rectangle (defined by coordinates) " << i+1 << " was rejected. The defined vertices are not co-planar.\n";
        //     std::cout << "Please check user defined coordinates for rectanle " << i+1 << " in input file\n";
        //     delete[] xProd1;
        //     delete[] xProd2;
        //     delete[] xProd3;
        //     continue; //go to next poly
        // }

        // Set normal vector
        new_poly.normal[0] = x_prod1[0]; //x
        new_poly.normal[1] = x_prod1[1]; //y
        new_poly.normal[2] = x_prod1[2]; //z
        normalize(&mut new_poly.normal);
        // Set radius (x and y radii might be switched based on order of users coordinates)
        new_poly.xradius = 0.5 * magnitude(v2[0], v2[1], v2[2]);
        new_poly.yradius = 0.5 * magnitude(v3[0], v3[1], v3[2]);
        new_poly.aspect_ratio = new_poly.yradius / new_poly.xradius;
        // Estimate translation
        // Use midpoint between 1st and 3rd vertices
        // Note: For polygons defined by coordinates, the coordinates
        // themselves provide the translation. We are just filling the
        // translation array for completeness even though the translation
        // array might not be used
        new_poly.translation[0] = 0.5 * (new_poly.vertices[0] + new_poly.vertices[6]);
        new_poly.translation[1] = 0.5 * (new_poly.vertices[1] + new_poly.vertices[7]);
        new_poly.translation[2] = 0.5 * (new_poly.vertices[2] + new_poly.vertices[8]);

        if domain_truncation(input, &mut new_poly, &input.domainSize) {
            // Poly completely outside domain
            new_poly.vertices.clear();
            pstats.rejection_reasons.outside += 1;
            pstats.rejected_poly_count += 1;
            println!("User Rectangle (defined by coordinates) {} was rejected for being outside the defined domain.", i + 1);
            rejected_user_fracture.id = i as isize + 1;
            rejected_user_fracture.user_fracture_type = -2;
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
                "User Defined Rectangular Fracture (Defined By Coordinates) {} Accepted",
                i + 1
            );
            accepted_poly.push(new_poly); // Save newPoly to accepted polys list
        } else {
            new_poly.vertices.clear(); // Need to delete manually, created with new[]
            pstats.rejects_per_attempt[pstats.accepted_poly_count] += 1;
            pstats.rejected_poly_count += 1;
            println!(
                "Rejected User Defined Rectangular Fracture (Defined By Coordinates) {}",
                i + 1
            );
            print_reject_reason(reject_code, &new_poly);
            rejected_user_fracture.id = i as isize + 1;
            rejected_user_fracture.user_fracture_type = -2;
            pstats.rejected_user_fracture.push(rejected_user_fracture);
        }
    } // End loop

    input.userRectCoordVertices.clear();
}
