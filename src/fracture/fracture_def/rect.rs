use parry3d_f64::na::{Point3, Vector3};
use tracing::{info, warn};

use crate::{
    computational_geometry::{create_bounding_box, intersection_checking},
    fracture::{domain::domain_truncation, insert_shape::print_reject_reason},
    math_functions::get_area,
    structures::{IntersectionPoints, Poly, RejectedUserFracture, Stats},
};

use super::UserDefinedFractures;

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
    domain_size: &Vector3<f64>,
    accepted_poly: &mut Vec<Poly>,
    intpts: &mut Vec<IntersectionPoints>,
    pstats: &mut Stats,
    triple_points: &mut Vec<Point3<f64>>,
    user_defined_rects: &UserDefinedFractures,
) {
    let npoly = user_defined_rects.n_frac;
    let family_id = -2;
    info!("{} User Rectangles Defined", npoly);

    let new_polys = user_defined_rects.create_polys(eps);

    for (i, mut new_poly) in new_polys.into_iter().enumerate() {
        if domain_truncation(h, eps, &mut new_poly, domain_size) {
            //poly completely outside domain
            pstats.rejection_reasons.outside += 1;
            pstats.rejected_poly_count += 1;
            warn!(
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
            info!("User Defined Rectangular Fracture {} Accepted", i + 1);
            accepted_poly.push(new_poly); // Save newPoly to accepted polys list
        } else {
            pstats.rejected_poly_count += 1;
            pstats.rejects_per_attempt[pstats.accepted_poly_count] += 1;
            info!("Rejected user defined rectangular fracture {}", i + 1);
            print_reject_reason(reject_code, &new_poly);
            pstats
                .rejected_user_fracture
                .push(RejectedUserFracture::new(i + 1, family_id));
        }
    }
}
