use crate::{
    computational_geometry::{
        apply_rotation2_d, apply_rotation3_d, create_bounding_box, intersection_checking, translate,
    },
    domain::domain_truncation,
    generating_points::generate_theta,
    insert_shape::{initialize_ell_vertices, print_reject_reason},
    math_functions::get_area,
    read_input::Input,
    structures::{IntPoints, Point, Poly, RejectedUserFracture, Stats},
    vector_functions::normalize,
};

// ***********************************************************************
// *********************  Insert User Ellipse  ***************************
// Inserts a user defined ellipse into the domain.
// Intersection checking, FRAM, and rejection/accptance are all called
// within this function.
// Arg 1: Array for all accepted polygons
// Arg 2: Array for all accepted intersections
// Arg 3: Program statistics structure
// Arg 4: Array of all triple intersection points
pub fn insert_user_ell(
    input: &mut Input,
    accepted_poly: &mut Vec<Poly>,
    intpts: &mut Vec<IntPoints>,
    pstats: &mut Stats,
    triple_points: &mut Vec<Point>,
) {
    println!("\n{} User Ellipses Defined\n", input.nUserEll);

    for i in 0..input.nUserEll {
        let index = i * 3; // Index to start of vertices/nodes
        let mut new_poly = Poly::default(); // New poly/fracture to be tested
        let mut rejected_user_fracture = RejectedUserFracture::default();
        new_poly.family_num = -1; // Using -1 for all user specified ellipses
        new_poly.vertices = Vec::with_capacity(input.uenumPoints[i] * 3);
        // Set number of nodes  - needed for rotations
        new_poly.number_of_nodes = input.uenumPoints[i] as isize;
        // Initialize translation data
        new_poly.translation[0] = input.uetranslation[index];
        new_poly.translation[1] = input.uetranslation[index + 1];
        new_poly.translation[2] = input.uetranslation[index + 2];
        // Generate theta array used to place vertices
        let mut theta_ary = Vec::new();
        generate_theta(&mut theta_ary, input.ueaspect[i], input.uenumPoints[i]);
        // Initialize vertices on x-y plane
        initialize_ell_vertices(
            &mut new_poly,
            input.ueRadii[i],
            input.ueaspect[i],
            &theta_ary,
            input.uenumPoints[i],
        );
        // Convert angle to rad if necessary
        let angle = if input.ueAngleOption {
            input.ueBeta[i] * std::f64::consts::PI / 180.
        } else {
            input.ueBeta[i]
        };

        // Initialize normal to {0,0,1}. need initialized for 3D rotation
        new_poly.normal[0] = 0.; //x
        new_poly.normal[1] = 0.; //y
        new_poly.normal[2] = 1.; //z
                                 // Apply 2d rotation matrix, twist around origin
                                 // Assumes polygon on x-y plane
                                 // Angle must be in rad
        apply_rotation2_d(&mut new_poly, angle);
        // Normalize user denined normal vector
        let mut tmp_uenormal = input.uenormal[index..index + 3].try_into().unwrap();
        normalize(&mut tmp_uenormal);
        for (i, j) in (index..index + 3).enumerate() {
            input.uenormal[j] = tmp_uenormal[i];
        }
        // Rotate vertices to uenormal[index] (new normal)
        apply_rotation3_d(input, &mut new_poly, &tmp_uenormal);
        // Save newPoly's new normal vector
        new_poly.normal[0] = input.uenormal[index];
        new_poly.normal[1] = input.uenormal[index + 1];
        new_poly.normal[2] = input.uenormal[index + 2];
        // Translate newPoly to uetranslation
        translate(
            &mut new_poly,
            input.uetranslation[index..index + 3].try_into().unwrap(),
        );

        if domain_truncation(input, &mut new_poly, &input.domainSize) {
            // Poly completely outside domain
            new_poly.vertices.clear();
            pstats.rejection_reasons.outside += 1;
            pstats.rejected_poly_count += 1;
            println!(
                "\nUser Ellipse {} was rejected for being outside the defined domain.",
                i + 1
            );
            rejected_user_fracture.id = i as isize + 1;
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
            //if intersection is ok
            if new_poly.truncated {
                pstats.truncated += 1;
            }

            // Incriment counter of accepted polys
            pstats.accepted_poly_count += 1;
            // Calculate poly's area
            new_poly.area = get_area(&new_poly);
            // Add new rejectsPerAttempt counter
            pstats.rejects_per_attempt.push(0);
            println!("User Defined Elliptical Fracture {} Accepted", i + 1);
            accepted_poly.push(new_poly); // Save newPoly to accepted polys list
        } else {
            new_poly.vertices.clear(); // Need to delete manually, created with new[]
            pstats.rejects_per_attempt[pstats.accepted_poly_count] += 1;
            pstats.rejected_poly_count += 1;
            println!("\nRejected User Defined Elliptical Fracture {}", i + 1);
            print_reject_reason(reject_code, &new_poly);
            rejected_user_fracture.id = i as isize + 1;
            rejected_user_fracture.user_fracture_type = -1;
            pstats.rejected_user_fracture.push(rejected_user_fracture);

            #[cfg(feature = "testing")]
            {
                panic!()
            }
        }

        println!("\n");
    }
}
