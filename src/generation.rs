use rand::distr::Uniform;
use rand::Rng;
use rand_mt::Mt64;
use tracing::{info, warn};

use crate::cg::intersection_checking;
use crate::error::DfngenError;
use crate::fracture::fracture_estimating::dry_run;
use crate::fracture::fracture_family::{FractureFamily, FractureFamilyCollection, RadiusOption};
use crate::fracture::insert_shape::{print_reject_reason, re_translate_poly};
use crate::fracture::insert_user_defined_fractures;
use crate::fracture::poly::Poly;
use crate::io::input::{Input, StopCondition};
use crate::math_functions::{
    cumsum, index_from_prob, index_from_prob_and_p32_status, p32_progress_fraction,
    remove_family_probability,
};
use crate::structures::{DFNGen, PolyOptions, RadiusFunction, RejectedFracture, Stats};

/// Result of a generation
pub struct Network {
    /// Input options. `stop_condition` and `n_poly` are those which the generation used.
    pub input: Input,
    /// Stochastic fracture families, with their final P32
    pub families: FractureFamilyCollection,
    /// Accepted and rejected fractures, intersections, triple points and statistics
    pub dfn: DFNGen,
    /// Sorted indices into `dfn.accepted_poly` of the fractures which remain
    /// after the cluster analysis
    pub final_fractures: Vec<usize>,
}

/// Progress of the insertion of the stochastic fractures
pub struct ProgressUpdate {
    /// Number of accepted fractures
    pub accepted_count: usize,
    /// Part of the stop condition which is complete, from 0 to 1
    pub fraction: f64,
}

/// Generates a fracture network
///
/// # Arguments
///
/// * `input` - Input options. `seed` is the seed of the random generator, also if it is 0.
/// * `frac_fam_col` - Stochastic fracture families
/// * `on_progress` - Called after each accepted stochastic fracture
pub fn generate(
    mut input: Input,
    mut frac_fam_col: FractureFamilyCollection,
    on_progress: &mut dyn FnMut(&ProgressUpdate),
) -> Result<Network, DfngenError> {
    let mut dfngen = DFNGen::new();
    let total_families = frac_fam_col.families.len();

    let poly_opts = PolyOptions {
        min_feature_size: input.min_feature_size,
        eps: input.eps,
        domain_size: input.domain_size,
        domain_mesh: input.domain_mesh.clone(),
        r_fram: input.r_fram,
        is_fram_enabled: input.is_fram_enabled,
        triple_intersections: input.triple_intersections,
        check_user_fractures: input.check_user_fractures,
        keep_rejected_fractures: input.output_rejected_fractures,
    };

    let mut generator = Mt64::new(input.seed);
    let dom_vol = input.domain_volume;

    if total_families > 0 {
        prepare_radii_lists(
            &input,
            &mut frac_fam_col,
            &mut dfngen.pstats,
            &mut generator,
        )?;
    }

    let num_user_defined_fractures =
        insert_user_defined_fractures(&input.user_fractures, &poly_opts, &mut dfngen);

    if total_families == 0 && input.stop_condition != StopCondition::NPoly {
        // If no stochastic shapes, use nPoly option with npoly = number of user polygons
        warn!("You have defined stop_condition = \"p32\" but have no stochastic fracture families defined.");
        warn!("Automatically setting stop_condition to \"n_poly\" for use with user defined fractures.");
        input.stop_condition = StopCondition::NPoly;

        if num_user_defined_fractures == 0 {
            return Err(DfngenError::NoFractureDefined);
        }

        // Set nPoly to the amount of user defined polygons
        input.n_poly = num_user_defined_fractures;
    }

    if total_families > 0 {
        insert_stochastic_fractures(
            &input,
            &poly_opts,
            &mut frac_fam_col,
            &mut dfngen,
            &mut generator,
            on_progress,
        )?;
    }

    // Remove last element off the rejects per attempt counter.
    // It will have one extra item due to how each element is initialized.
    dfngen.pstats.rejects_per_attempt.pop();

    // The close to node check is inside of the close to edge check function
    // for optimization (only need do intersection close to node on one condition)
    // On close to node rejections, close to edge is counted as well.
    // To get the correct number we must subtract the close to node count
    // (they were counted in closeToEdge AND closeToNode)
    dfngen.pstats.rejection_reasons.close_to_edge -= dfngen.pstats.rejection_reasons.close_to_node;

    info!("Network Generation Complete");
    info!("Version of dfngen-rs: {}", env!("CARGO_PKG_VERSION"));

    if input.stop_condition == StopCondition::P32 {
        info!("Final p32 values per family:");

        for (i, shape) in frac_fam_col.families.iter().enumerate() {
            info!(
                "Family {} target P32 = {}, Final P32 = {}",
                i + 1,
                shape.p32_target,
                shape.current_p32
            );
        }
    }

    report_statistics_before_removal(&input, &frac_fam_col, &mut dfngen, dom_vol);

    if input.remove_fractures_less_than > 0. {
        info!(
            "Removing fractures with radius less than {} and rebuilding DFN",
            input.remove_fractures_less_than
        );
        let size = dfngen.accepted_poly.len();
        dfngen.remove_small_fractures(&poly_opts, input.remove_fractures_less_than);
        info!(
            "Removed {} fractures with radius less than {}",
            size - dfngen.accepted_poly.len(),
            input.remove_fractures_less_than
        );
    }

    // Remove any isolated fractures and return
    // a list of polygon indices. If input option
    // keepOnlyLargestCluster = 1, return largest
    // cluster. DFN will keep all fractures with intersections.
    let mut final_fractures = dfngen.clusters.final_fractures(
        input.keep_isolated_fractures,
        input.keep_only_largest_cluster,
        dfngen.pstats.accepted_poly_count,
    );
    // Sort fracture indices to retain order by acceptance
    final_fractures.sort();

    if final_fractures.is_empty() {
        return Err(DfngenError::NoIntersectingFractures);
    }

    report_statistics_after_removal(
        &input,
        &frac_fam_col,
        &mut dfngen,
        &final_fractures,
        dom_vol,
    );

    report_intersection_statistics(&dfngen);
    report_rejection_statistics(&dfngen.pstats);

    if total_families > 0 {
        report_estimation_statistics(&frac_fam_col, &dfngen.pstats);
    }

    info!("Seed: {}", input.seed);

    Ok(Network {
        input,
        families: frac_fam_col,
        dfn: dfngen,
        final_fractures,
    })
}

/// Pre-generates the sorted radii list of each family and initializes the
/// per-family counters in `pstats`.
fn prepare_radii_lists(
    input: &Input,
    frac_fam_col: &mut FractureFamilyCollection,
    pstats: &mut Stats,
    generator: &mut Mt64,
) -> Result<(), DfngenError> {
    let total_families = frac_fam_col.families.len();

    if input.stop_condition == StopCondition::NPoly {
        // Npoly Option
        // Estimate fractures, generate radii lists for nPoly option
        info!("Building radii lists for nPoly option...");
        frac_fam_col.generate_radii(input.force_large_fractures, input.n_poly, generator)?;
        info!("Building radii lists for nPoly option Complete");
    } else if input.is_fram_enabled {
        // P32 Option
        // ESTIMATE # FRACTURES NEEDED
        info!("Estimating number of fractures needed...");
        dry_run(input, frac_fam_col, generator)?;
    }

    // Add a percentage more radii to each radii
    // list using families' distribution.
    // First arg is percentage, eg: 0.1 will add 10% more fractures
    // to the radii list for each family
    if input.is_fram_enabled {
        frac_fam_col.generate_radii(
            false,
            (input.n_poly as f64 * input.radii_list_increase).ceil() as usize,
            generator,
        )?;

        for fracture_family in &frac_fam_col.families {
            info!(
                "Estimated {} fractures for family {}",
                fracture_family.radii_list.len(),
                fracture_family.id
            );
        }
    }

    frac_fam_col.sort_radii();

    // Keep count of accepted & rejected fractures by family
    pstats.accepted_from_fam = vec![0; total_families];
    pstats.rejected_from_fam = vec![0; total_families];

    // Save sizes of pre-generated radii lists per family.
    // Print as part of statistics to user
    pstats.expected_from_fam.extend(
        frac_fam_col
            .families
            .iter()
            .map(|fracture_family| fracture_family.radii_list.len()),
    );

    Ok(())
}

/// Inserts fractures from the stochastic families until the stop condition
/// (`nPoly` count or P32 target of every family) is met.
fn insert_stochastic_fractures(
    input: &Input,
    poly_opts: &PolyOptions,
    frac_fam_col: &mut FractureFamilyCollection,
    dfngen: &mut DFNGen,
    generator: &mut Mt64,
    on_progress: &mut dyn FnMut(&ProgressUpdate),
) -> Result<(), DfngenError> {
    let total_families = frac_fam_col.families.len();
    let dom_vol = input.domain_volume;

    // Probabilities of the families which have not met their P32 target (P32 option).
    // This is a copy: the family probabilities must stay in their original state.
    let mut probabilities = frac_fam_col.probabilities.clone();
    let mut cdf = cumsum(&probabilities);

    // Initialize uniform distribution on [0,1]
    let uniform_dist = Uniform::new(0., 1.)?;

    // Used with stopCondition = 1, P32 option.
    // Values are set to true once the families p32 requirement is met.
    // Once all elements have values all set to true, all families have had their
    // P32 requirement
    let mut p32_status = vec![false; total_families];

    while (input.stop_condition == StopCondition::NPoly
        && dfngen.pstats.accepted_poly_count < input.n_poly)
        || (input.stop_condition == StopCondition::P32 && p32_status.iter().any(|p| !p))
    {
        // cdfIdx holds the index to the CDF array for the current shape family being inserted
        let roll = generator.sample(uniform_dist);

        // cdf_idx is the index to the CDF array of the family being inserted
        let (family_index, cdf_idx) = if input.stop_condition == StopCondition::NPoly {
            // nPoly Option
            // Choose a family based purely on famProb probabilities
            let family_index = index_from_prob(&cdf, roll);
            (family_index, family_index)
        } else {
            // P32 Option
            // Choose a family based on probabiliyis AND their target p32 completion status.
            // If a family has already met is fracture intinisty req. (p32) don't choose that family anymore
            index_from_prob_and_p32_status(&p32_status, &cdf, roll)
        };

        let mut new_poly = frac_fam_col.families[family_index].create_poly(
            input.eps,
            family_index,
            RadiusOption::FromCacheOrRng,
            generator,
        )?;

        let mut reject_counter = 0;

        // Re-translate the polygon until it is accepted or reaches its reject limit
        loop {
            if new_poly.check_box_domain(input.min_feature_size, input.eps, &input.domain_size) {
                // Poly was completely outside domain
                dfngen.pstats.rejection_reasons.outside += 1;

                if reject_counter >= input.rejects_per_fracture {
                    // Reject poly, generate new polygon
                    break;
                }

                // Retranslate poly and try again, preserving normal, size, and shape
                re_translate_poly(
                    &mut new_poly,
                    &frac_fam_col.families[family_index],
                    generator,
                );
                continue;
            }

            new_poly.assign_bounding_box();
            // Find line of intersection and FRAM check
            let fram_result = intersection_checking(poly_opts, &mut new_poly, dfngen);

            let Err(reject_reason) = fram_result else {
                dfngen.pstats.accepted_poly_count += 1;
                dfngen.pstats.accepted_from_fam[family_index] += 1;
                // Make new rejection counter for next fracture attempt
                dfngen.pstats.rejects_per_attempt.push(0);

                new_poly.assign_area();

                let family = &mut frac_fam_col.families[family_index];
                family.add_area_to_p32(new_poly.area, dom_vol, &input.region_volumes);

                // If the last inserted pologon met the p32 reqirement, set that familiy to no longer
                // insert any more fractures. ajust the CDF and familiy probabilites to account for this
                if input.stop_condition == StopCondition::P32
                    && family.current_p32 >= family.p32_target
                {
                    p32_status[family_index] = true;
                    info!("P32 For Family {} Completed", family_index + 1);

                    // Remove the completed family's element in 'CDF[]' and 'famProb[]'
                    // Distribute the removed family probability evenly among the others and rebuild the CDF
                    if probabilities.len() > 1 {
                        cdf = remove_family_probability(&mut probabilities, cdf_idx);
                    }
                }

                on_progress(&ProgressUpdate {
                    accepted_count: dfngen.pstats.accepted_poly_count,
                    fraction: progress_fraction(
                        input,
                        dfngen.pstats.accepted_poly_count,
                        &frac_fam_col.families,
                    ),
                });

                // Intersections and triple points are saved within intersection_checking()
                dfngen.accepted_poly.push(new_poly);
                break;
            };

            reject_counter += 1;
            // Inc reject counter for current attempt
            // (number of rejects until next fracture accepted)
            dfngen.pstats.rejects_per_attempt[dfngen.pstats.accepted_poly_count] += 1;

            if input.print_reject_reasons {
                print_reject_reason(reject_reason, &new_poly);
            }

            if reject_counter >= input.rejects_per_fracture {
                dfngen.pstats.rejected_poly_count += 1;
                dfngen.pstats.rejected_from_fam[family_index] += 1;

                if poly_opts.keep_rejected_fractures {
                    dfngen.rejected_poly.push(RejectedFracture {
                        poly: new_poly,
                        reason: reject_reason,
                    });
                }

                // Stop retranslating polygon if its reached its reject limit
                break;
            }

            if input.print_reject_reasons {
                info!("Translating rejected fracture to new position");
            }

            dfngen.pstats.retranslated_poly_count += 1;
            re_translate_poly(
                &mut new_poly,
                &frac_fam_col.families[family_index],
                generator,
            );
        }
    }

    Ok(())
}

/// Part of the stop condition which is complete, from 0 to 1
fn progress_fraction(input: &Input, accepted_count: usize, families: &[FractureFamily]) -> f64 {
    match input.stop_condition {
        StopCondition::NPoly => accepted_count as f64 / input.n_poly as f64,
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

fn log_total_area_statistics(total_area: f64, fracture_count: usize, dom_vol: f64) {
    info!("Total Surface Area: {} m^2", total_area * 2.);
    info!(
        "Total Fracture Density (P30): {}",
        fracture_count as f64 / dom_vol
    );
    info!(
        "Total Fracture Intensity (P32): {}",
        (total_area * 2.) / dom_vol
    );
}

/// Logs the region a family belongs to.
fn log_family_location(frac_fam: &FractureFamily, input: &Input) {
    match frac_fam.region_id {
        0 => info!("Region: Whole Domain"),
        region_id => info!(
            "Region: {} of the domain mesh",
            input.domain_mesh.regions()[region_id - 1].name
        ),
    }
}

fn report_statistics_before_removal(
    input: &Input,
    frac_fam_col: &FractureFamilyCollection,
    dfngen: &mut DFNGen,
    dom_vol: f64,
) {
    info!("Statistics Before Isolated Fractures Removed:",);
    info!("Fractures: {}", dfngen.accepted_poly.len());

    let areas = sum_areas(dfngen.accepted_poly.iter(), frac_fam_col.families.len());
    dfngen.pstats.area_before_removal = areas.total;

    log_total_area_statistics(areas.total, dfngen.accepted_poly.len(), dom_vol);

    for (i, frac_fam) in frac_fam_col.families.iter().enumerate() {
        info!("Family: {}", i + 1);
        info!("Accepted: {}", dfngen.pstats.accepted_from_fam[i]);
        info!("Rejected: {}", dfngen.pstats.rejected_from_fam[i]);
        log_family_location(frac_fam, input);
        info!("Surface Area: {} m^2", areas.per_family[i] * 2.);
        info!("Fracture Intensity (P32): {}", frac_fam.current_p32);
    }

    if areas.user_defined > 0. {
        info!("User Defined:");
        info!("Surface Area: {} m^2", areas.user_defined * 2.);
        info!(
            "Fracture Intensity (P32): {}",
            areas.user_defined * 2. / dom_vol
        );
    }
}

fn report_statistics_after_removal(
    input: &Input,
    frac_fam_col: &FractureFamilyCollection,
    dfngen: &mut DFNGen,
    final_fractures: &[usize],
    dom_vol: f64,
) {
    let total_families = frac_fam_col.families.len();

    info!("Statistics After Isolated Fractures Removed:");
    info!("Final Number of Fractures: {}", final_fractures.len());
    info!(
        "Isolated Fractures Removed: {}",
        dfngen.accepted_poly.len() - final_fractures.len()
    );
    info!(
        "Fractures before isolated fractures removed: {}",
        dfngen.accepted_poly.len()
    );

    let areas = sum_areas(
        final_fractures.iter().map(|&i| &dfngen.accepted_poly[i]),
        total_families,
    );
    dfngen.pstats.area_after_removal = areas.total;

    // Re-count number of accepted fracture per family after isloated fractures were removed
    let mut final_count_per_family = vec![0; total_families];
    for &i in final_fractures {
        if let Some(family_index) = dfngen.accepted_poly[i].origin.family_index() {
            final_count_per_family[family_index] += 1;
        }
    }

    log_total_area_statistics(areas.total, final_fractures.len(), dom_vol);

    for (i, frac_fam) in frac_fam_col.families.iter().enumerate() {
        info!("Family: {}", i + 1);
        info!(
            "Fractures After Isolated Fracture Removal: {}",
            final_count_per_family[i]
        );
        info!(
            "Isolated Fractures Removed: {}",
            dfngen.pstats.accepted_from_fam[i] - final_count_per_family[i]
        );
        info!("Accepted: {}", dfngen.pstats.accepted_from_fam[i]);
        info!("Rejected: {}", dfngen.pstats.rejected_from_fam[i]);
        log_family_location(frac_fam, input);
        info!("Surface Area: {} m^2", areas.per_family[i] * 2.);
        info!(
            "Fracture Intensity (P32): {}",
            areas.per_family[i] * 2. / dom_vol
        );
    }

    if areas.user_defined > 0. {
        info!("User Defined Shapes:");
        info!("Surface Area: {} m^2", areas.user_defined * 2.);
        info!(
            "Fracture Intensity (P32): {}",
            areas.user_defined * 2. / dom_vol
        );
    }

    info!(
        "{} Fractures Accepted (Before Isolated Fracture Removal)",
        dfngen.accepted_poly.len()
    );
    info!(
        "{} Final Fractures (After Isolated Fracture Removal)",
        final_fractures.len()
    );
    info!(
        "Total Fractures Rejected: {}",
        dfngen.pstats.rejected_poly_count
    );
    info!(
        "Total Fractures Re-translated: {}",
        dfngen.pstats.retranslated_poly_count
    );
}

fn report_intersection_statistics(dfngen: &DFNGen) {
    let pstats = &dfngen.pstats;

    info!(
        "Number of Triple Intersection Points (Before Isolated Fracture Removal): {}",
        dfngen.triple_points.len()
    );
    info!("Intersection Statistics:");
    info!("Number of Intersections: {}", dfngen.intpts.len());
    info!(
        "Intersections Shortened: {}",
        pstats.intersections_shortened
    );
    info!(
        "Original Intersection (Before Intersection Shrinking) Length: {} m",
        pstats.original_length
    );
    info!(
        "Intersection Length Discarded: {} m",
        pstats.discarded_length
    );
    info!(
        "Final Intersection Length: {} m",
        pstats.original_length - pstats.discarded_length
    );
}

fn report_rejection_statistics(pstats: &Stats) {
    let reasons = &pstats.rejection_reasons;

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
}

/// Compares the estimated and the actual number of fractures per family.
fn report_estimation_statistics(frac_fam_col: &FractureFamilyCollection, pstats: &Stats) {
    info!("Fracture Estimation statistics:");
    info!("NOTE: If estimation and actual are very different, expected family distributions might not be accurate.");
    info!("If this is the case, try increasing or decreasing the 'generation.radii_list_increase' option in the input file.");

    for (i, shape) in frac_fam_col.families.iter().enumerate() {
        match shape.radius.function {
            RadiusFunction::Constant(_) => {
                info!("Family {} Using constant size", shape.id);
            }
            _ => {
                info!(
                    "Family {} Estimated: {}",
                    shape.id, pstats.expected_from_fam[i]
                );

                info!(
                    "Actual: {}",
                    pstats.accepted_from_fam[i] + pstats.rejected_from_fam[i]
                );
            }
        }
    }
}
