use tracing::{info, warn};

use crate::{
    computational_geometry::intersection_checking,
    structures::{DFNGen, PolyOptions},
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
pub fn remove_fractures(min_size: f64, opts: PolyOptions, dfngen: &mut DFNGen) {
    let mut final_poly_list = Vec::new();
    // Clear GroupData
    dfngen.pstats.group_data.clear();
    // Clear FractGroup
    dfngen.pstats.fract_group.clear();
    // Clear Triple Points
    dfngen.triple_points.clear();
    // Clear IntPoints
    dfngen.intpts.clear();
    // Re-init nextGroupNum
    dfngen.pstats.next_group_num = 1;

    for poly in dfngen.accepted_poly.iter_mut() {
        if poly.xradius < min_size {
            poly.vertices.clear();
            continue;
        }

        let mut new_poly = poly.clone();
        new_poly.group_num = 0; // Reset cluster group number
        new_poly.intersection_index.clear(); // Remove ref to old intersections
                                             // Find line of intersection and FRAM check
        let reject_code = intersection_checking(
            opts.h,
            opts.eps,
            opts.r_fram,
            opts.disable_fram,
            opts.triple_intersections,
            &mut new_poly,
            &mut final_poly_list,
            &mut dfngen.intpts,
            &mut dfngen.pstats,
            &mut dfngen.triple_points,
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
    dfngen.accepted_poly.clear();
    dfngen.accepted_poly.extend(final_poly_list);
}
