use parry3d_f64::na::{Point3, Vector3};

use super::domain::domain_truncation;
use super::insert_shape::{initialize_rect_vertices, print_reject_reason};

use crate::{
    computational_geometry::{
        apply_rotation2_d, apply_rotation3_d, create_bounding_box, intersection_checking, translate,
    },
    math_functions::get_area,
    structures::{IntersectionPoints, Poly, RejectedUserFracture, Stats},
};

#[allow(clippy::too_many_arguments)]
fn create_poly(
    eps: f64,
    ur_radii: &[f64],
    uraspect: &[f64],
    ur_beta: &[f64],
    urnormal: &mut [f64],
    urtranslation: &[f64],
    ur_angle_option: bool,
    idx: usize,
) -> Poly {
    let mut new_poly = Poly {
        family_num: -2,
        // Set number of nodes. Needed for rotations.
        number_of_nodes: 4,
        // Initialize normal to {0,0,1}. need initialized for 3D rotation
        normal: Vector3::new(0., 0., 1.),
        ..Default::default()
    };

    new_poly.vertices.reserve(12); // 4*{x,y,z}

    // Index to start of vertices/nodes
    let index = idx * 3;

    // initializeRectVertices() sets newpoly.xradius, newpoly.yradius, newpoly.aperture
    initialize_rect_vertices(&mut new_poly, ur_radii[idx], uraspect[idx]);

    // Convert angle to rad if necessary
    let angle = if ur_angle_option {
        ur_beta[idx] * std::f64::consts::PI / 180.
    } else {
        ur_beta[idx]
    };

    // Apply 2d rotation matrix, twist around origin
    // Assumes polygon on x-y plane
    // Angle must be in rad
    apply_rotation2_d(&mut new_poly, angle);

    // Rotate vertices to urnormal[index] (new normal)
    let normal = Vector3::from_row_slice(&urnormal[index..index + 3]).normalize();
    apply_rotation3_d(&mut new_poly, &normal, eps);

    // Update normal vector
    new_poly.normal = normal;
    urnormal[index] = normal.x;
    urnormal[index + 1] = normal.y;
    urnormal[index + 2] = normal.z;

    // Translate newPoly to urtranslation
    translate(
        &mut new_poly,
        Vector3::new(
            urtranslation[index],
            urtranslation[index + 1],
            urtranslation[index + 2],
        ),
    );

    new_poly
}

/// Insert User Rectangles
///
/// Inserts a user defined rectangle into the domain
/// Intersection checking, FRAM, and rejection/accptance is all contained
/// within this function.
///
/// # Arguments
///
/// * `h` - Minimum feature size
/// * `eps` - Epsilon value for floating point comparisons
/// * `disable_fram` - If true, FRAM is disabled
/// * `triple_intersections` - If true, triple intersections are accepted
/// * `n_user_rect` - Number of user defined rectangles
/// * `domain_size` - Size of the domain
/// * `ur_radii` - Array of radii for user defined rectangles
/// * `uraspect` - Array of aspect ratios for user defined rectangles
/// * `ur_beta` - User defined rectangles beta array
/// * `urnormal` - User defined rectangles normal array
/// * `urtranslation` - User defined rectangles translation array
/// * `ur_angle_option` - If true, angle is in degrees else in radians
/// * `accepted_poly` - Array for all accepted polygons
/// * `intpts` - Array for all accepted intersections
/// * `pstats` - Program statistics structure
/// * `triple_points` - Array of all triple intersection points
#[allow(clippy::too_many_arguments)]
pub fn insert_user_rects(
    h: f64,
    eps: f64,
    r_fram: bool,
    disable_fram: bool,
    triple_intersections: bool,
    n_user_rect: usize,
    domain_size: &Vector3<f64>,
    ur_radii: &[f64],
    uraspect: &[f64],
    ur_beta: &[f64],
    urnormal: &mut [f64],
    urtranslation: &[f64],
    ur_angle_option: bool,
    accepted_poly: &mut Vec<Poly>,
    intpts: &mut Vec<IntersectionPoints>,
    pstats: &mut Stats,
    triple_points: &mut Vec<Point3<f64>>,
) {
    let npoly = n_user_rect;
    let family_id = -2;
    println!("{} User Rectangles Defined", npoly);

    for i in 0..npoly {
        let mut new_poly = create_poly(
            eps,
            ur_radii,
            uraspect,
            ur_beta,
            urnormal,
            urtranslation,
            ur_angle_option,
            i,
        );

        if domain_truncation(h, eps, &mut new_poly, domain_size) {
            //poly completely outside domain
            pstats.rejection_reasons.outside += 1;
            pstats.rejected_poly_count += 1;
            println!(
                "User Rectangle {} was rejected for being outside the defined domain.",
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
            pstats.rejected_poly_count += 1;
            pstats.rejects_per_attempt[pstats.accepted_poly_count] += 1;
            println!("Rejected user defined rectangular fracture {}", i + 1);
            print_reject_reason(reject_code, &new_poly);
            pstats
                .rejected_user_fracture
                .push(RejectedUserFracture::new(i + 1, family_id));
        }
    }
}
