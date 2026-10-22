use rand::distr::Uniform;
use rand::Rng;
use rand_mt::Mt64;
use tracing::{debug, info};

use crate::dfn::{Dfn, PolyOptions, RejectedFracture};
use crate::distribution::RadiusFunction;
use crate::error::DfngenError;
use crate::fracture::family_probability::{
    cumulative_sum, index_from_probability, index_from_probability_and_p32_status,
    p32_progress_fraction, remove_family_probability,
};
use crate::fracture::fracture_family::{FractureFamily, FractureFamilyCollection};
use crate::fracture::insert_user_defined_fractures;
use crate::fracture::poly::{FractureOrigin, Poly};
use crate::fram::{check_intersections, log_reject_reason};
use crate::io::input::{Input, StopCondition};
use crate::mesh::{mesh_network, FractureMesh};
use crate::region_boundaries::boundary_lines;
use crate::stats::{FractureCount, Stats};

/// Result of a generation
pub struct Network {
    /// Input options
    pub input: Input,
    /// Stochastic fracture families, with their final P32
    pub families: FractureFamilyCollection,
    /// Accepted and rejected fractures, intersections, triple points and statistics
    pub dfn: Dfn,
    /// Sorted indices into `dfn.accepted_polys` of the fractures which remain
    /// after the cluster analysis
    pub final_fractures: Vec<usize>,
    /// Triangle mesh of the final fractures. `None` without `[mesh]`.
    pub mesh: Option<FractureMesh>,
}

/// Progress of the insertion of the stochastic fractures
pub struct ProgressUpdate {
    /// Number of accepted stochastic fractures
    pub accepted_count: usize,
    /// Part of the stop condition which is complete, from 0 to 1
    pub fraction: f64,
}

/// Generates a fracture network
///
/// # Arguments
///
/// * `input` - Input options. `seed` is the seed of the random generator, also if it is 0.
/// * `family_collection` - Stochastic fracture families
/// * `on_progress` - Called after each accepted stochastic fracture
pub fn generate(
    input: Input,
    mut family_collection: FractureFamilyCollection,
    on_progress: &mut dyn FnMut(&ProgressUpdate),
) -> Result<Network, DfngenError> {
    let mut dfn = Dfn::new();
    let total_families = family_collection.families.len();

    let options = PolyOptions {
        min_feature_size: input.min_feature_size,
        eps: input.eps,
        domain_size: input.domain_size,
        domain_mesh: input.domain_mesh.clone(),
        fram: input.fram.clone(),
        region_boundaries: input.region_boundaries.clone(),
    };

    let mut generator = Mt64::new(input.seed);
    let domain_volume = input.domain_volume;

    if total_families > 0 {
        prepare_radii_lists(
            &input,
            &mut family_collection,
            &mut dfn.stats,
            &mut generator,
        )?;
    }

    if total_families == 0 && input.user_fractures.is_empty() {
        return Err(DfngenError::NoFractureDefined);
    }

    insert_user_defined_fractures(&input.user_fractures, &options, &mut dfn);

    if total_families > 0 {
        insert_stochastic_fractures(
            &input,
            &options,
            &mut family_collection,
            &mut dfn,
            &mut generator,
            on_progress,
        )?;
    }

    // Remove last element off the rejects per attempt counter.
    // It will have one extra item due to how each element is initialized.
    dfn.stats.rejects_per_attempt.pop();

    // The close to node check is inside of the close to edge check function
    // for optimization (only need do intersection close to node on one condition)
    // On close to node rejections, close to edge is counted as well.
    // To get the correct number we must subtract the close to node count
    // (they were counted in `close_to_edge` and `close_to_node`)
    dfn.stats.rejection_reasons.close_to_edge -= dfn.stats.rejection_reasons.close_to_node;

    info!("Network Generation Complete");
    info!("Version of dfngen-rs: {}", env!("CARGO_PKG_VERSION"));

    if input.stop_condition == StopCondition::P32 && total_families > 0 {
        info!("Final p32 values per family:");

        for (i, family) in family_collection.families.iter().enumerate() {
            info!(
                "Family {} target P32 = {}, Final P32 = {}",
                i + 1,
                family.p32_target,
                family.current_p32
            );
        }
    }

    report_statistics_before_removal(&input, &family_collection, &mut dfn, domain_volume);

    if input.remove_fractures_less_than > 0. {
        info!(
            "Removing fractures with radius less than {} and rebuilding DFN",
            input.remove_fractures_less_than
        );
        let size = dfn.accepted_polys.len();
        dfn.remove_small_fractures(&options, input.remove_fractures_less_than);
        info!(
            "Removed {} fractures with radius less than {}",
            size - dfn.accepted_polys.len(),
            input.remove_fractures_less_than
        );
    }

    // Remove any isolated fractures and return
    // a list of polygon indices. If input option
    // `keep_only_largest_cluster` is true, return largest
    // cluster. DFN will keep all fractures with intersections.
    let mut final_fractures = dfn.clusters.final_fractures(
        input.keep_isolated_fractures,
        input.keep_only_largest_cluster,
        dfn.stats.accepted_count,
    );
    // Sort fracture indices to retain order by acceptance
    final_fractures.sort();

    if final_fractures.is_empty() {
        return Err(DfngenError::NoIntersectingFractures);
    }

    report_statistics_after_removal(
        &input,
        &family_collection,
        &mut dfn,
        &final_fractures,
        domain_volume,
    );

    report_intersection_statistics(&dfn);
    report_rejection_statistics(&input, &dfn.stats);

    if total_families > 0 {
        report_estimation_statistics(&family_collection, &dfn.stats);
    }

    info!("Seed: {}", input.seed);

    let mesh = input
        .mesh
        .as_ref()
        .map(|options| mesh_network(&input, options, &dfn, &final_fractures))
        .transpose()?;

    Ok(Network {
        input,
        families: family_collection,
        dfn,
        final_fractures,
        mesh,
    })
}

/// Pre-generates the sorted radii list of each family and initializes the
/// per-family counters in `stats`.
fn prepare_radii_lists(
    input: &Input,
    family_collection: &mut FractureFamilyCollection,
    stats: &mut Stats,
    generator: &mut Mt64,
) -> Result<(), DfngenError> {
    let total_families = family_collection.families.len();

    match input.stop_condition {
        StopCondition::NPoly(n_poly) => {
            // Npoly Option
            // Estimate fractures, generate radii lists for n_poly option
            info!("Building radii lists for nPoly option...");
            family_collection.generate_radii(input.force_large_fractures, n_poly, generator)?;
            info!("Building radii lists for nPoly option Complete");
        }
        StopCondition::P32 => {
            info!("Building radii lists for the P32 target of each family...");
            family_collection.generate_p32_radii(
                input.force_large_fractures,
                input.domain_volume,
                &input.region_volumes,
                generator,
            )?;
        }
    }

    if input.fram.is_some() {
        for fracture_family in &family_collection.families {
            info!(
                "Estimated {} fractures for family {}",
                fracture_family.radii_list.len(),
                fracture_family.id
            );
        }
    }

    family_collection.sort_radii();

    // Keep count of accepted & rejected fractures by family
    stats.accepted_per_family = vec![0; total_families];
    stats.rejected_per_family = vec![0; total_families];

    // Save sizes of pre-generated radii lists per family.
    // Print as part of statistics to user
    stats.expected_per_family.extend(
        family_collection
            .families
            .iter()
            .map(|fracture_family| fracture_family.radii_list.len()),
    );

    Ok(())
}

/// Inserts fractures from the stochastic families until the stop condition
/// (`n_poly` count or P32 target of every family) is met.
fn insert_stochastic_fractures(
    input: &Input,
    options: &PolyOptions,
    family_collection: &mut FractureFamilyCollection,
    dfn: &mut Dfn,
    generator: &mut Mt64,
    on_progress: &mut dyn FnMut(&ProgressUpdate),
) -> Result<(), DfngenError> {
    let total_families = family_collection.families.len();
    let domain_volume = input.domain_volume;

    // Probabilities of the families which have not met their P32 target (P32 option).
    // This is a copy: the family probabilities must stay in their original state.
    let mut probabilities = family_collection.probabilities.clone();
    let mut cdf = cumulative_sum(&probabilities);

    // Initialize uniform distribution on [0,1]
    let uniform_dist = Uniform::new(0., 1.)?;

    // Used with the P32 stop condition.
    // Values are set to true once the families p32 requirement is met.
    // Once all elements have values all set to true, all families have had their
    // P32 requirement
    let mut p32_status = vec![false; total_families];

    // `n_poly` does not include the user-defined fractures
    let mut stochastic_accepted_count = 0;

    while match input.stop_condition {
        StopCondition::NPoly(n_poly) => stochastic_accepted_count < n_poly,
        StopCondition::P32 => p32_status.iter().any(|p| !p),
    } {
        let roll = generator.sample(uniform_dist);

        // cdf_index is the index to the CDF array of the family being inserted
        let (family_index, cdf_index) = match input.stop_condition {
            StopCondition::NPoly(_) => {
                // n_poly Option
                // Choose a family based purely on the family probabilities
                let family_index = index_from_probability(&cdf, roll);
                (family_index, family_index)
            }
            StopCondition::P32 => {
                // P32 Option
                // Choose a family based on probabiliyis AND their target p32 completion status.
                // If a family has already met is fracture intinisty req. (p32) don't choose that family anymore
                index_from_probability_and_p32_status(&p32_status, &cdf, roll)
            }
        };

        let mut new_poly = family_collection.families[family_index].create_poly(
            input.eps,
            family_index,
            generator,
        )?;

        let region_ids = boundary_region_ids(input, &family_collection.families[family_index]);
        new_poly.region_ids = region_ids.clone();
        let mut reject_counter = 0;

        // Re-translate the polygon until it is accepted or reaches its reject limit
        loop {
            if new_poly.is_outside_box_domain(input.min_feature_size, input.eps, &input.domain_size)
            {
                // Poly was completely outside domain. The center of a fracture
                // is in the mesh, thus this does not occur.
                dfn.stats.rejection_reasons.outside += 1;

                // Retranslate poly and try again, preserving normal, size, and family
                family_collection.families[family_index].retranslate(&mut new_poly, generator);
                continue;
            }

            new_poly.assign_bounding_box();

            if options.region_boundaries.is_some() {
                new_poly.boundary_lines = boundary_lines(
                    &new_poly.vertices,
                    &new_poly.normal,
                    &options.domain_mesh,
                    &region_ids,
                );
            }

            // Find line of intersection and FRAM check
            let fram_result = check_intersections(options, &mut new_poly, dfn);

            // A reject comes only from FRAM, thus FRAM has its options
            let (Err(reject_reason), Some(fram_options)) = (fram_result, &input.fram) else {
                dfn.stats.accepted_count += 1;
                stochastic_accepted_count += 1;
                dfn.stats.accepted_per_family[family_index] += 1;
                // Make new rejection counter for next fracture attempt
                dfn.stats.rejects_per_attempt.push(0);

                new_poly.assign_area();

                let family = &mut family_collection.families[family_index];
                family.add_area_to_p32(new_poly.area, domain_volume, &input.region_volumes);

                // If the last inserted pologon met the p32 reqirement, set that familiy to no longer
                // insert any more fractures. ajust the CDF and familiy probabilites to account for this
                if input.stop_condition == StopCondition::P32
                    && family.current_p32 >= family.p32_target
                {
                    p32_status[family_index] = true;
                    info!("P32 For Family {} Completed", family_index + 1);

                    // Remove the completed family's element in the CDF and the probabilities
                    // Distribute the removed family probability evenly among the others and rebuild the CDF
                    if probabilities.len() > 1 {
                        cdf = remove_family_probability(&mut probabilities, cdf_index);
                    }
                }

                on_progress(&ProgressUpdate {
                    accepted_count: stochastic_accepted_count,
                    fraction: progress_fraction(
                        input,
                        stochastic_accepted_count,
                        &family_collection.families,
                    ),
                });

                // Intersections and triple points are saved within check_intersections()
                dfn.accepted_polys.push(new_poly);
                break;
            };

            reject_counter += 1;
            // Inc reject counter for current attempt
            // (number of rejects until next fracture accepted)
            dfn.stats.rejects_per_attempt[dfn.stats.accepted_count] += 1;

            log_reject_reason(reject_reason, family_index);

            if reject_counter >= fram_options.rejects_per_fracture {
                dfn.stats.rejected_count += 1;
                dfn.stats.rejected_per_family[family_index] += 1;

                if fram_options.write_rejected_fractures {
                    dfn.rejected_polys.push(RejectedFracture {
                        poly: new_poly,
                        reason: reject_reason,
                    });
                }

                // Stop retranslating polygon if its reached its reject limit
                break;
            }

            debug!("Translating rejected fracture to new position");

            dfn.stats.retranslated_count += 1;
            family_collection.families[family_index].retranslate(&mut new_poly, generator);
        }
    }

    Ok(())
}

/// Regions of a stochastic fracture: its boundary lines are at their
/// boundary, and its mesh has no triangle outside them. No region: all
/// cells of the mesh.
fn boundary_region_ids(input: &Input, family: &FractureFamily) -> Vec<usize> {
    match &input.region_boundaries {
        Some(options) if options.use_family_region && family.region_id > 0 => {
            vec![family.region_id]
        }
        _ => Vec::new(),
    }
}

/// Part of the stop condition which is complete, from 0 to 1
fn progress_fraction(input: &Input, accepted_count: usize, families: &[FractureFamily]) -> f64 {
    match input.stop_condition {
        StopCondition::NPoly(n_poly) => accepted_count as f64 / n_poly as f64,
        StopCondition::P32 => {
            let current_and_target_p32: Vec<(f64, f64)> = families
                .iter()
                .map(|family| (family.current_p32, family.p32_target))
                .collect();

            p32_progress_fraction(&current_and_target_p32)
        }
    }
}

/// Fracture area (one side) summed over a set of fractures.
struct AreaSummary {
    total: f64,
    /// One element per stochastic family
    per_family: Vec<f64>,
    user_defined: f64,
}

fn sum_areas<'a>(polys: impl Iterator<Item = &'a Poly>, total_families: usize) -> AreaSummary {
    let mut summary = AreaSummary {
        total: 0.,
        per_family: vec![0.; total_families],
        user_defined: 0.,
    };

    for poly in polys {
        summary.total += poly.area;

        match poly.origin.family_index() {
            Some(family_index) => summary.per_family[family_index] += poly.area,
            None => summary.user_defined += poly.area,
        }
    }

    summary
}

/// Number of fractures and polygons. The pieces of a user-defined fracture
/// (same `fracture_id`) are one fracture.
fn count_fractures<'a>(polys: impl IntoIterator<Item = &'a Poly>) -> FractureCount {
    FractureCount::from_fracture_ids(polys.into_iter().map(|poly| poly.user_fracture_id))
}

fn count_user_defined_fractures<'a>(polys: impl IntoIterator<Item = &'a Poly>) -> FractureCount {
    count_fractures(
        polys
            .into_iter()
            .filter(|poly| poly.origin == FractureOrigin::UserDefined),
    )
}

fn log_total_area_statistics(total_area: f64, fracture_count: usize, domain_volume: f64) {
    info!("Total Surface Area: {} m^2", total_area * 2.);
    info!(
        "Total Fracture Density (P30): {}",
        fracture_count as f64 / domain_volume
    );
    info!(
        "Total Fracture Intensity (P32): {}",
        (total_area * 2.) / domain_volume
    );
}

/// Logs the region a family belongs to.
fn log_family_location(family: &FractureFamily, input: &Input) {
    match family.region_id {
        0 => info!("Region: Whole Domain"),
        region_id => info!(
            "Region: {} of the domain mesh",
            input.domain_mesh.regions()[region_id - 1].name
        ),
    }
}

fn report_statistics_before_removal(
    input: &Input,
    family_collection: &FractureFamilyCollection,
    dfn: &mut Dfn,
    domain_volume: f64,
) {
    let accepted_count = count_fractures(&dfn.accepted_polys);
    info!("Statistics Before Isolated Fractures Removed:",);
    info!("Fractures: {}", accepted_count);

    let areas = sum_areas(dfn.accepted_polys.iter(), family_collection.families.len());
    dfn.stats.area_before_removal = areas.total;

    log_total_area_statistics(areas.total, accepted_count.fractures, domain_volume);

    for (i, family) in family_collection.families.iter().enumerate() {
        info!("Family: {}", i + 1);
        info!("Accepted: {}", dfn.stats.accepted_per_family[i]);
        info!("Rejected: {}", dfn.stats.rejected_per_family[i]);
        log_family_location(family, input);
        info!("Surface Area: {} m^2", areas.per_family[i] * 2.);
        info!("Fracture Intensity (P32): {}", family.current_p32);
    }

    if !input.user_fractures.is_empty() {
        info!("User Defined:");
        info!(
            "Accepted: {}",
            count_user_defined_fractures(&dfn.accepted_polys)
        );
        info!("Rejected: {}", dfn.stats.rejected_user_fracture_count);
        info!("Surface Area: {} m^2", areas.user_defined * 2.);
        info!(
            "Fracture Intensity (P32): {}",
            areas.user_defined * 2. / domain_volume
        );
    }
}

fn report_statistics_after_removal(
    input: &Input,
    family_collection: &FractureFamilyCollection,
    dfn: &mut Dfn,
    final_fractures: &[usize],
    domain_volume: f64,
) {
    let total_families = family_collection.families.len();

    let final_polys: Vec<&Poly> = final_fractures
        .iter()
        .map(|&i| &dfn.accepted_polys[i])
        .collect();
    let accepted_count = count_fractures(&dfn.accepted_polys);
    let final_count = count_fractures(final_polys.iter().copied());
    let final_user_defined_count = count_user_defined_fractures(final_polys.iter().copied());
    let areas = sum_areas(final_polys.iter().copied(), total_families);

    // Re-count number of accepted fracture per family after isloated fractures were removed
    let mut final_count_per_family = vec![0; total_families];
    for poly in &final_polys {
        if let Some(family_index) = poly.origin.family_index() {
            final_count_per_family[family_index] += 1;
        }
    }

    dfn.stats.area_after_removal = areas.total;

    info!("Statistics After Isolated Fractures Removed:");
    info!("Final Number of Fractures: {}", final_count);
    info!(
        "Isolated Fractures Removed: {}",
        accepted_count - final_count
    );
    info!(
        "Fractures before isolated fractures removed: {}",
        accepted_count
    );

    log_total_area_statistics(areas.total, final_count.fractures, domain_volume);

    for (i, family) in family_collection.families.iter().enumerate() {
        info!("Family: {}", i + 1);
        info!(
            "Fractures After Isolated Fracture Removal: {}",
            final_count_per_family[i]
        );
        info!(
            "Isolated Fractures Removed: {}",
            dfn.stats.accepted_per_family[i] - final_count_per_family[i]
        );
        info!("Accepted: {}", dfn.stats.accepted_per_family[i]);
        info!("Rejected: {}", dfn.stats.rejected_per_family[i]);
        log_family_location(family, input);
        info!("Surface Area: {} m^2", areas.per_family[i] * 2.);
        info!(
            "Fracture Intensity (P32): {}",
            areas.per_family[i] * 2. / domain_volume
        );
    }

    if !input.user_fractures.is_empty() {
        info!("User Defined Shapes:");
        info!("Final: {}", final_user_defined_count);
        info!("Surface Area: {} m^2", areas.user_defined * 2.);
        info!(
            "Fracture Intensity (P32): {}",
            areas.user_defined * 2. / domain_volume
        );
    }

    let stochastic_rejected_count: usize = dfn.stats.rejected_per_family.iter().sum();
    let rejected_count = FractureCount::of_single_polygons(stochastic_rejected_count)
        + dfn.stats.rejected_user_fracture_count;

    info!(
        "Fractures Accepted (Before Isolated Fracture Removal): {}",
        accepted_count
    );
    info!(
        "Final Fractures (After Isolated Fracture Removal): {}",
        final_count
    );
    info!("Total Fractures Rejected: {}", rejected_count);
    info!(
        "Total Fractures Re-translated: {}",
        dfn.stats.retranslated_count
    );
}

fn report_intersection_statistics(dfn: &Dfn) {
    let stats = &dfn.stats;

    info!(
        "Number of Triple Intersection Points (Before Isolated Fracture Removal): {}",
        dfn.triple_points.len()
    );
    info!("Intersection Statistics:");
    info!("Number of Intersections: {}", dfn.intersections.len());
    info!("Intersections Shortened: {}", stats.intersections_shortened);
    info!(
        "Original Intersection (Before Intersection Shrinking) Length: {} m",
        stats.original_length
    );
    info!(
        "Intersection Length Discarded: {} m",
        stats.discarded_length
    );
    info!(
        "Final Intersection Length: {} m",
        stats.original_length - stats.discarded_length
    );
}

fn report_rejection_statistics(input: &Input, stats: &Stats) {
    let reasons = &stats.rejection_reasons;

    info!("Rejection Statistics:");
    info!("{} Short Intersections", reasons.short_intersection);
    info!("{} Close to Node", reasons.close_to_node);
    info!("{} Close to Edge", reasons.close_to_edge);
    info!("{} Vertex Close to Edge", reasons.close_point_to_edge);
    info!("{} Outside of Domain", reasons.outside);
    info!("{} Triple intersection Rejections", reasons.triple);
    info!(
        "{} Intersections Close to Other Intersections",
        reasons.inter_close_to_inter
    );

    if input.region_boundaries.is_some() {
        info!("{} Boundary Lines of the Regions", reasons.boundary_line);
    }
}

/// Compares the estimated and the actual number of fractures per family.
fn report_estimation_statistics(family_collection: &FractureFamilyCollection, stats: &Stats) {
    info!("Fracture Estimation statistics:");
    info!("NOTE: If estimation and actual are very different, expected family distributions might not be accurate.");
    for (i, family) in family_collection.families.iter().enumerate() {
        match family.radius.function {
            RadiusFunction::Constant(_) => {
                info!("Family {} Using constant size", family.id);
            }
            _ => {
                info!(
                    "Family {} Estimated: {}",
                    family.id, stats.expected_per_family[i]
                );

                info!(
                    "Actual: {}",
                    stats.accepted_per_family[i] + stats.rejected_per_family[i]
                );
            }
        }
    }
}
