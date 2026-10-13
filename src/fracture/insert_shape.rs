use rand_mt::Mt64;
use tracing::warn;

use crate::fracture::fracture_family::FractureFamily;
use crate::structures::RejectReason;

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
pub fn re_translate_poly(new_poly: &mut Poly, frac_fam: &FractureFamily, generator: &mut Mt64) {
    new_poly.group_num = 0; // Clear cluster group information
    new_poly.intersection_index.clear(); // Clear any saved intersections

    // Move poly back to origin
    for vertex in &mut new_poly.vertices {
        *vertex -= new_poly.translation;
    }

    // Translate to new position
    let t = frac_fam.sample_position(generator);
    new_poly.translate(t);
}

/// Print Rejection Reson
///
/// Function prints the fracture reject reason to the user
///
/// # Arguments
///
/// * `reject_reason` - Reason why FRAM rejected the poly
/// * `new_poly` - Poly which was rejected
pub fn print_reject_reason(reject_reason: RejectReason, new_poly: &Poly) {
    if let Some(family_index) = new_poly.origin.family_index() {
        warn!(
            "Attempted fracture from family {} was rejected:",
            family_index
        );
    }

    warn!(
        "\trejectCode = {}: {}",
        reject_reason.code(),
        reject_reason.description()
    );
}
