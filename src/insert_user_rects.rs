use crate::{
    computational_geometry::{
        apply_rotation2_d, apply_rotation3_d, create_bounding_box, intersection_checking, translate,
    },
    domain::domain_truncation,
    insert_shape::{initialize_rect_vertices, print_reject_reason},
    math_functions::get_area,
    read_input::Input,
    structures::{IntPoints, Point, Poly, RejectedUserFracture, Stats},
    vector_functions::normalize,
};

// ***********************************************************************
// ********************  Insert User Rectangles  *************************
//     Inserts a user defined rectangle into the domain
//     Intersection checking, FRAM, and rejection/accptance is all contained
//     within this function.
//     Arg 1: Array for all accepted polygons
//     Arg 2: Array for all accepted intersections
//     Arg 3: Program statistics structure
//     Arg 4: Array of all triple intersection points
pub fn insert_user_rects(
    input: &mut Input,
    accepted_poly: &mut Vec<Poly>,
    intpts: &mut Vec<IntPoints>,
    pstats: &mut Stats,
    triple_points: &mut Vec<Point>,
) {
    println!("{} User Rectangles Defined", input.nUserRect);

    for i in 0..input.nUserRect {
        let mut new_poly = Poly::default();
        let mut rejected_user_fracture = RejectedUserFracture::default();
        new_poly.family_num = -2; // Using -2 for all user specified rectangles
        new_poly.vertices.reserve(12); // 4*{x,y,z}
                                       // Set number of nodes. Needed for rotations.
        new_poly.number_of_nodes = 4;
        let index = i * 3; // Index to start of vertices/nodes
                           // initializeRectVertices() sets newpoly.xradius, newpoly.yradius, newpoly.aperture
        initialize_rect_vertices(&mut new_poly, input.urRadii[i], input.uraspect[i]);

        // Convert angle to rad if necessary
        let angle = if input.urAngleOption {
            input.urBeta[i] * std::f64::consts::PI / 180.
        } else {
            input.urBeta[i]
        };

        // Initialize normal to {0,0,1}. need initialized for 3D rotation
        new_poly.normal[0] = 0.; //x
        new_poly.normal[1] = 0.; //y
        new_poly.normal[2] = 1.; //z
                                 // Apply 2d rotation matrix, twist around origin
                                 // Assumes polygon on x-y plane
                                 // Angle must be in rad
        apply_rotation2_d(&mut new_poly, angle);
        // Rotate into 3D from poly.normal to "urnormal", new normal
        normalize(&mut input.urnormal[index..index + 3].try_into().unwrap());
        // Rotate vertices to urnormal[index] (new normal)
        apply_rotation3_d(
            input,
            &mut new_poly,
            &input.urnormal[index..index + 3].try_into().unwrap(),
        );
        // Save newPoly's new normal vector
        new_poly.normal[0] = input.urnormal[index];
        new_poly.normal[1] = input.urnormal[index + 1];
        new_poly.normal[2] = input.urnormal[index + 2];
        // Translate newPoly to urtranslation
        translate(
            &mut new_poly,
            &input.urtranslation[index..index + 3].try_into().unwrap(),
        );

        if domain_truncation(input, &mut new_poly, &input.domainSize) {
            //poly completely outside domain
            new_poly.vertices.clear();
            pstats.rejection_reasons.outside += 1;
            pstats.rejected_poly_count += 1;
            println!(
                "User Rectangle {} was rejected for being outside the defined domain.",
                i + 1
            );
            rejected_user_fracture.id = (i + 1) as isize;
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
            // If intersection is ok
            if new_poly.truncated {
                pstats.truncated += 1;
            }

            // Incriment counter of accepted polys
            pstats.accepted_poly_count += 1;
            // Calculate poly's area
            new_poly.area = get_area(&new_poly);
            // Add new rejectsPerAttempt counter
            pstats.rejects_per_attempt.push(0);
            println!("User Defined Rectangular Fracture {} Accepted", i + 1);
            accepted_poly.push(new_poly); // Save newPoly to accepted polys list
        } else {
            new_poly.vertices.clear(); // Delete manually, created with new[]
            pstats.rejected_poly_count += 1;
            pstats.rejects_per_attempt[pstats.accepted_poly_count] += 1;
            println!("Rejected user defined rectangular fracture {}", i + 1);
            print_reject_reason(reject_code, &new_poly);
            rejected_user_fracture.id = (i + 1) as isize;
            rejected_user_fracture.user_fracture_type = -2;
            pstats.rejected_user_fracture.push(rejected_user_fracture);
        }
    }
}
