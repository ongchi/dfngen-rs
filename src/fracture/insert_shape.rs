use std::cell::RefCell;
use std::rc::Rc;

use rand_mt::Mt64;
use tracing::warn;

use crate::fracture::fracture_family::FractureFamily;

use super::poly::Poly;

/// Retranslate Polygon
///
/// Translates the polygon to a new position, preserving its size, shape,
/// and normal vector
/// This helps hit target distributions since we reject less
///
/// # Arguments
///
/// * `new_poly` - Polygon
/// * `frac_fam` - Fracture family structure which Polygon belongs to
/// * `generator` - Random Generator
pub fn re_translate_poly(
    new_poly: &mut Poly,
    frac_fam: &FractureFamily,
    generator: Rc<RefCell<Mt64>>,
) {
    new_poly.group_num = 0; // Clear cluster group information
    new_poly.intersection_index.clear(); // Clear any saved intersections

    // Move poly back to origin
    for vertex in &mut new_poly.vertices {
        *vertex -= new_poly.translation;
    }

    // Translate to new position
    let t = frac_fam.sample_position(generator.clone());
    new_poly.translate(t);
}

/// Print Rejection Reson
///
/// Function prints fracture rejection reasons to user based on reject code
/// Currtenly only used for user defined fractures.
///
/// # Arguments
///
/// * `reject_code` - Rejection code
/// * `new_poly` - Poly which was rejected
pub fn print_reject_reason(reject_code: i32, new_poly: &Poly) {
    if new_poly.family_num >= 0 {
        warn!(
            "Attempted fracture from family {} was rejected:",
            new_poly.family_num
        );
    }

    match reject_code {
        -2 => {
            warn!("\trejectCode = -2: Intersection of length < h.");
        }

        -1 => {
            warn!("\trejectCode = -1: Fracture too close to a node.")
        }

        -6 => {
            warn!("\trejectCode = -6: Fracture too close to another fracture's edge.");
        }

        -7 => {
            warn!("\trejectCode = -7: Fractures intersecting on same plane");
        }

        -10 => {
            warn!("\trejectCode = -10: Rejected triple intersection due to triple intersections being turned off in input file.")
        }

        -11 => {
            warn!("\trejectCode = -11: Fracture's intersection landed too close to a previous intersection.");
        }

        -12 => {
            warn!("\trejectCode = -12: Fracture created a triple intersection with an angle too small.");
        }

        -13 => {
            warn!("\trejectCode = -13: Fracture created a triple intersection with the triple intersection point too close to an intersection's endpoint.");
        }

        -14 => {
            warn!("\trejectCode = -14: Fracture created a triple intersection with the triple intersection point too close to another triple intersection point.");
        }

        _ => {
            warn!("\trejectCode = {}", reject_code);
        }
    }
}
