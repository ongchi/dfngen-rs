use parry3d_f64::na::{Point3, Vector3};

use super::domain::domain_truncation;
use super::insert_shape::{initialize_ell_vertices, print_reject_reason};

use crate::{
    computational_geometry::{
        apply_rotation2_d, apply_rotation3_d, create_bounding_box, intersection_checking, translate,
    },
    distribution::generating_points::generate_theta,
    io::input::Input,
    math_functions::get_area,
    structures::{IntersectionPoints, Poly, RejectedUserFracture, Stats},
};

fn create_poly(input: &mut Input, idx: usize) -> Poly {
    let index = idx * 3;

    let mut new_poly = Poly {
        family_num: -1,
        // Set number of nodes. Needed for rotations.
        number_of_nodes: input.uenumPoints[idx] as isize,
        // Initialize normal to {0,0,1}. need initialized for 3D rotation
        normal: Vector3::new(0., 0., 1.),
        translation: Vector3::from_row_slice(&input.uetranslation[index..index + 3]),
        ..Default::default()
    };

    new_poly.vertices = Vec::with_capacity(input.uenumPoints[idx] * 3);

    // Generate theta array used to place vertices
    let mut theta_ary = Vec::new();
    generate_theta(&mut theta_ary, input.ueaspect[idx], input.uenumPoints[idx]);

    // Initialize vertices on x-y plane
    initialize_ell_vertices(
        &mut new_poly,
        input.ueRadii[idx],
        input.ueaspect[idx],
        &theta_ary,
        input.uenumPoints[idx],
    );

    // Convert angle to rad if necessary
    let angle = if input.ueAngleOption {
        input.ueBeta[idx] * std::f64::consts::PI / 180.
    } else {
        input.ueBeta[idx]
    };

    // Apply 2d rotation matrix, twist around origin
    // Assumes polygon on x-y plane
    // Angle must be in rad
    apply_rotation2_d(&mut new_poly, angle);
    // Normalize user denined normal vector
    let normal = Vector3::from_row_slice(&input.uenormal[index..index + 3]).normalize();

    // Rotate vertices to uenormal[index] (new normal)
    apply_rotation3_d(
        &mut new_poly,
        &Vector3::from_row_slice(&input.urnormal[index..index + 3]),
        input.eps,
    );

    // Save newPoly's new normal vector
    new_poly.normal = normal;
    input.uenormal[index] = normal.x;
    input.uenormal[index + 1] = normal.y;
    input.uenormal[index + 2] = normal.z;

    // Translate newPoly to uetranslation
    translate(
        &mut new_poly,
        Vector3::new(
            input.uetranslation[index],
            input.uetranslation[index + 1],
            input.uetranslation[index + 2],
        ),
    );

    new_poly
}

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
    intpts: &mut Vec<IntersectionPoints>,
    pstats: &mut Stats,
    triple_points: &mut Vec<Point3<f64>>,
) {
    let family_id = -1;
    let npoly = input.nUserEll;
    println!("{} User Ellipses Defined", npoly);

    for i in 0..npoly {
        let mut new_poly = create_poly(input, i);

        if domain_truncation(input.h, input.eps, &mut new_poly, &input.domainSize) {
            // Poly completely outside domain
            pstats.rejection_reasons.outside += 1;
            pstats.rejected_poly_count += 1;
            println!(
                "\nUser Ellipse {} was rejected for being outside the defined domain.",
                i + 1
            );
            pstats
                .rejected_user_fracture
                .push(RejectedUserFracture::new(i + 1, family_id));
            continue; // Go to next poly (go to next iteration of for loop)
        }

        create_bounding_box(&mut new_poly);
        // Line of intersection and FRAM
        let reject_code = intersection_checking(
            input.h,
            input.eps,
            input.rFram,
            input.disableFram,
            input.tripleIntersections,
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
            pstats.rejects_per_attempt[pstats.accepted_poly_count] += 1;
            pstats.rejected_poly_count += 1;
            println!("\nRejected User Defined Elliptical Fracture {}", i + 1);
            print_reject_reason(reject_code, &new_poly);
            pstats
                .rejected_user_fracture
                .push(RejectedUserFracture::new(i + 1, family_id));
        }
    }
}
