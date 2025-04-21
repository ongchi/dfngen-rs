use parry3d_f64::na::Point3;
use tracing::{info, warn};

use crate::{
    computational_geometry::intersection_checking,
    structures::{IntersectionPoints, Poly, Stats},
};

/// Remove Fractures Smaller Than Minimum Size
///
/// Function is designed to be used AFTER DFN generation
/// Originally created to compare the difference in distributions
/// if small fractures were removed after a DFN was generated
/// as opposed limiting their insertion during DFN generation.
///
/// The minimum size options in the input files will still be used.
/// If the user wishes to also removeFractures after DFN generation,
/// 'minSize' here must be larger than the minimum size fractures in the
/// input file. Fratrues with radii less than 'minSize' will be removed
/// after DFN has been created.
///
/// # Arguments
///
/// * `h` - Minimum feature size
/// * `eps` - Epsilon value for floating point comparisons
/// * `r_fram` - Uses a relaxed version of the FRAM algorithm. The mesh may not be perfectly conforming
/// * `disable_fram` - If true, FRAM is disabled
/// * `triple_intersections` - If true, triple intersections are accepted
/// * `min_size` - Minimum size. All fractures smaller than 'minSize' will be removed.
/// * `accepted_polys` - Array of accepted polygons
/// * `int_pts` - Array of accepted intersections
/// * `triple_points` - Array of all triple intersection points
/// * `pstats` - Stats structure (DFN Statisctics)
///
/// NOTE: Must be executed before getCluster()
///       This funciton rebuilds the DFN. Using getCluster() before this
///       funciton executes causes undefined behavior.
#[allow(clippy::too_many_arguments)]
pub fn remove_fractures(
    h: f64,
    eps: f64,
    r_fram: bool,
    disable_fram: bool,
    triple_intersections: bool,
    min_size: f64,
    accepted_polys: &mut Vec<Poly>,
    int_pts: &mut Vec<IntersectionPoints>,
    triple_points: &mut Vec<Point3<f64>>,
    pstats: &mut Stats,
) {
    let mut final_poly_list = Vec::new();
    // Clear GroupData
    pstats.group_data.clear();
    // Clear FractGroup
    pstats.fract_group.clear();
    // Clear Triple Points
    triple_points.clear();
    // Clear IntPoints
    int_pts.clear();
    // Re-init nextGroupNum
    pstats.next_group_num = 1;

    for poly in accepted_polys.iter_mut() {
        if poly.xradius < min_size {
            poly.vertices.clear();
            continue;
        }

        let mut new_poly = poly.clone();
        new_poly.group_num = 0; // Reset cluster group number
        new_poly.intersection_index.clear(); // Remove ref to old intersections
                                             // Find line of intersection and FRAM check
        let reject_code = intersection_checking(
            h,
            eps,
            r_fram,
            disable_fram,
            triple_intersections,
            &mut new_poly,
            &mut final_poly_list,
            int_pts,
            pstats,
            triple_points,
        );

        // IF POLY ACCEPTED:
        if reject_code == 0 {
            // Intersections are ok
            // SAVING POLYGON (intersection and triple points saved witchin intersectionChecking())
            final_poly_list.push(new_poly); // SAVE newPoly to accepted polys list
        } else {
            // Poly rejected
            warn!("Error rebuilding dfn, previously accepted fracture was rejected during DFN rebuild.");
        }
    }

    info!("Rebuilding DFN complete.");
    accepted_polys.clear();
    accepted_polys.extend(final_poly_list);
}
