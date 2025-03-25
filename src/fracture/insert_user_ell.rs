use parry3d_f64::na::{Point3, Vector3};

use super::domain::domain_truncation;
use super::insert_shape::{initialize_ell_vertices, print_reject_reason};

use crate::{
    computational_geometry::{
        apply_rotation2_d, apply_rotation3_d, create_bounding_box, intersection_checking, translate,
    },
    distribution::generating_points::generate_theta,
    math_functions::get_area,
    structures::{IntersectionPoints, Poly, RejectedUserFracture, Stats},
};

#[allow(clippy::too_many_arguments)]
fn create_poly(
    eps: f64,
    ue_angle_option: bool,
    uenum_points: &[usize],
    uetranslation: &[f64],
    ueaspect: &[f64],
    ue_radii: &[f64],
    ue_beta: &[f64],
    uenormal: &mut [f64],
    urnormal: &[f64],
    idx: usize,
) -> Poly {
    let index = idx * 3;

    let mut new_poly = Poly {
        family_num: -1,
        // Set number of nodes. Needed for rotations.
        number_of_nodes: uenum_points[idx] as isize,
        // Initialize normal to {0,0,1}. need initialized for 3D rotation
        normal: Vector3::new(0., 0., 1.),
        translation: Vector3::from_row_slice(&uetranslation[index..index + 3]),
        ..Default::default()
    };

    new_poly.vertices = Vec::with_capacity(uenum_points[idx] * 3);

    // Generate theta array used to place vertices
    let mut theta_ary = Vec::new();
    generate_theta(&mut theta_ary, ueaspect[idx], uenum_points[idx]);

    // Initialize vertices on x-y plane
    initialize_ell_vertices(
        &mut new_poly,
        ue_radii[idx],
        ueaspect[idx],
        &theta_ary,
        uenum_points[idx],
    );

    // Convert angle to rad if necessary
    let angle = if ue_angle_option {
        ue_beta[idx] * std::f64::consts::PI / 180.
    } else {
        ue_beta[idx]
    };

    // Apply 2d rotation matrix, twist around origin
    // Assumes polygon on x-y plane
    // Angle must be in rad
    apply_rotation2_d(&mut new_poly, angle);
    // Normalize user denined normal vector
    let normal = Vector3::from_row_slice(&uenormal[index..index + 3]).normalize();

    // Rotate vertices to uenormal[index] (new normal)
    apply_rotation3_d(
        &mut new_poly,
        &Vector3::from_row_slice(&urnormal[index..index + 3]),
        eps,
    );

    // Save newPoly's new normal vector
    new_poly.normal = normal;
    uenormal[index] = normal.x;
    uenormal[index + 1] = normal.y;
    uenormal[index + 2] = normal.z;

    // Translate newPoly to uetranslation
    translate(
        &mut new_poly,
        Vector3::new(
            uetranslation[index],
            uetranslation[index + 1],
            uetranslation[index + 2],
        ),
    );

    new_poly
}

/// Insert User Ellipse
///
/// Inserts a user defined ellipse into the domain.
/// Intersection checking, FRAM, and rejection/accptance are all called
/// within this function.
///
/// # Arguments
///
/// * `h` - Minimum feature size
/// * `eps` - Epsilon value for floating point comparisons
/// * `n_user_ell` - Number of user defined ellipses
/// * `ue_angle_option` - Option to use angle in degrees or radians
/// * `r_fram` - Uses a relaxed version of the FRAM algorithm. The mesh may not be perfectly conforming
/// * `disable_fram` - If true, FRAM is disabled
/// * `triple_intersections` - If true, triple intersections are allowed
/// * `domain_size` - Size of the domain
/// * `uenum_points` - User ellipses number of points per ellipse array.
/// * `uetranslation` - User ellipses translation array.
/// * `ueaspect` - User ellipses aspect ratio array.
/// * `ue_radii` - User ellipses radii array.
/// * `ue_beta` - User ellipses beta array.
/// * `uenormal` - User ellipses normal vector array.
/// * `urnormal` - User rectangles normal vector array.
/// * `accepted_poly` - Array for all accepted polygons
/// * `intpts` - Array for all accepted intersections
/// * `pstats` - Program statistics structure
/// * `triple_points` - Array of all triple intersection points
#[allow(clippy::too_many_arguments)]
pub fn insert_user_ell(
    h: f64,
    eps: f64,
    n_user_ell: usize,
    ue_angle_option: bool,
    r_fram: bool,
    disable_fram: bool,
    triple_intersections: bool,
    domain_size: &Vector3<f64>,
    uenum_points: &[usize],
    uetranslation: &[f64],
    ueaspect: &[f64],
    ue_radii: &[f64],
    ue_beta: &[f64],
    uenormal: &mut [f64],
    urnormal: &[f64],
    accepted_poly: &mut Vec<Poly>,
    intpts: &mut Vec<IntersectionPoints>,
    pstats: &mut Stats,
    triple_points: &mut Vec<Point3<f64>>,
) {
    let family_id = -1;
    let npoly = n_user_ell;
    println!("{} User Ellipses Defined", npoly);

    for i in 0..npoly {
        let mut new_poly = create_poly(
            eps,
            ue_angle_option,
            uenum_points,
            uetranslation,
            ueaspect,
            ue_radii,
            ue_beta,
            uenormal,
            urnormal,
            i,
        );

        if domain_truncation(h, eps, &mut new_poly, domain_size) {
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
            h,
            eps,
            r_fram,
            disable_fram,
            triple_intersections,
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
