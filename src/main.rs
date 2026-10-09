use std::cell::RefCell;
use std::fs::File;
use std::io::Write;
use std::rc::Rc;
use std::time::{SystemTime, UNIX_EPOCH};

use clap::Parser;
use rand::distr::Uniform;
use rand::Rng;
use rand_mt::Mt64;
use tracing::{error, info, warn};
use tracing_subscriber::{layer::SubscriberExt, util::SubscriberInitExt};

use crate::cg::intersection_checking;
use crate::error::DfngenError;
use crate::fracture::cluster_groups::get_cluster;
use crate::fracture::fracture_estimating::dry_run;
use crate::fracture::fracture_family::{FractureFamily, FractureFamilyCollection, RadiusOption};
use crate::fracture::insert_shape::{print_reject_reason, re_translate_poly};
use crate::fracture::insert_user_defined_fractures;
use crate::fracture::poly::Poly;
use crate::io::input::{read_input, Input};
use crate::io::output::write_output;
use crate::math_functions::{
    adjust_cdf_and_fam_prob, cumsum, index_from_prob, index_from_prob_and_p32_status,
};
use crate::structures::{DFNGen, PolyOptions, RadiusFunction, Stats};

mod cg;
mod distribution;
mod error;
mod fracture;
mod io;
mod math_functions;
mod structures;

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    /// Path to input file
    input_file: String,

    /// Path to output folder
    output_folder: String,
}

fn main() -> Result<(), DfngenError> {
    let cli = Cli::parse();

    std::fs::create_dir_all(&cli.output_folder)?;
    init_tracing(File::create(format!(
        "{}/DFN_output.txt",
        cli.output_folder
    ))?);

    info!("Starting dfngen-rs");

    let mut dfngen = DFNGen::new();

    let (mut input, mut frac_fam_col) = read_input(&cli.input_file);
    let total_families = frac_fam_col.families.len();

    let poly_opts = PolyOptions {
        h: input.h,
        eps: input.eps,
        domain_size: input.domainSize,
        r_fram: input.rFram,
        disable_fram: input.disableFram,
        triple_intersections: input.tripleIntersections,
    };

    let generator = Rc::new(RefCell::new(Mt64::new(match input.seed {
        0 => SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .map(|t| t.as_secs())?,
        seed => seed,
    })));
    let dom_vol = domain_volume(&input);

    if total_families > 0 {
        prepare_radii_lists(
            &input,
            &mut frac_fam_col,
            &mut dfngen.pstats,
            generator.clone(),
        )?;
    }

    let num_user_defined_fractures = insert_user_defined_fractures(
        input.insertUserRectanglesFirst,
        &input.ext_fracture_files,
        &poly_opts,
        &mut dfngen,
    );

    if total_families == 0 && input.stopCondition != 0 {
        // If no stochastic shapes, use nPoly option with npoly = number of user polygons
        warn!("You have defined stopCondition = 1 (P32 program stopping condition) but have no stochastic shape families defined.");
        warn!("Automatically setting stopCondition to 0 for use with user defined polygons and nPoly.");
        input.stopCondition = 0;

        if num_user_defined_fractures == 0 {
            panic!("ERROR: No fracture defined, please check input file for errors.");
        }

        // Set nPoly to the amount of user defined polygons
        input.nPoly = num_user_defined_fractures;
    }

    if total_families > 0 {
        insert_stochastic_fractures(
            &input,
            &cli.output_folder,
            &poly_opts,
            &mut frac_fam_col,
            &mut dfngen,
            generator,
        )?;
    }

    // The close to node check is inside of the close to edge check function
    // for optimization (only need do intersection close to node on one condition)
    // On close to node rejections, close to edge is counted as well.
    // To get the correct number we must subtract the close to node count
    // (they were counted in closeToEdge AND closeToNode)
    dfngen.pstats.rejection_reasons.close_to_edge -= dfngen.pstats.rejection_reasons.close_to_node;

    info!("Network Generation Complete");
    info!("Version of dfngen-rs: {}", env!("CARGO_PKG_VERSION"));

    if input.stopCondition == 1 {
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

    if input.removeFracturesLessThan > 0. {
        info!(
            "Removing fractures with radius less than {} and rebuilding DFN",
            input.removeFracturesLessThan
        );
        let size = dfngen.accepted_poly.len();
        dfngen.remove_small_fractures(&poly_opts, input.removeFracturesLessThan);
        info!(
            "Removed {} fractures with radius less than {}",
            size - dfngen.accepted_poly.len(),
            input.removeFracturesLessThan
        );
    }

    if input.polygonBoundaryFlag {
        info!("Extracting fractures from a polygon boundary domain");
        let size = dfngen.accepted_poly.len();
        dfngen.remove_fractures_beyond_domains(&poly_opts, &input.domainVertices);
        info!(
            "Removed {} fractures outside subdomain",
            size - dfngen.accepted_poly.len()
        );
    }

    // Remove any isolated fractures and return
    // a list of polygon indices matching the users
    // boundaryFaces option. If input option
    // keepOnlyLargestCluster = 1, return largest
    // cluster matching users boundaryFaces option
    // If ignoreBoundaryFaces input option is on,
    // DFN will keep all fractures with intersections.
    let mut final_fractures = get_cluster(
        input.keepIsolatedFractures,
        input.keepOnlyLargestCluster,
        input.ignoreBoundaryFaces,
        &input.boundaryFaces,
        &dfngen.pstats,
    );
    // Sort fracture indices to retain order by acceptance
    final_fractures.sort();
    // Error check for no boundary connection
    let mut print_connectivity_error = false;

    if final_fractures.is_empty() && !input.ignoreBoundaryFaces {
        print_connectivity_error = true;
        //if there is no fracture network connected users defined boundary faces
        //switch to ignore boundary faces option with notice to user that there is no connectivity
        final_fractures = get_cluster(
            input.keepIsolatedFractures,
            input.keepOnlyLargestCluster,
            input.ignoreBoundaryFaces,
            &input.boundaryFaces,
            &dfngen.pstats,
        );
        //if still no fractures, there is no fracture network
    }

    if final_fractures.is_empty() {
        error!("DFN Generation has finished, however there are no intersecting fractures. Please adjust input parameters.");
        error!("Try increasing the fracture density, or shrinking the domain.",);
        std::process::exit(0)
    }

    report_statistics_after_removal(
        &input,
        &frac_fam_col,
        &mut dfngen,
        &final_fractures,
        dom_vol,
    );

    if print_connectivity_error {
        error!("DFN generation has finished but the formed",);
        error!("fracture network does not make a connection between",);
        error!("the user's specified boundary faces.");
        error!("Try increasing the fracture density, shrinking the domain",);
        error!("or consider using the 'ignoreBoundaryFaces' option.",);
        panic!();
    }

    report_intersection_statistics(&dfngen);
    report_rejection_statistics(&dfngen.pstats);

    if total_families > 0 {
        report_estimation_statistics(&frac_fam_col, &dfngen.pstats);
    }

    info!("Seed: {}", input.seed);

    // Write all output files
    write_output(
        &input,
        &cli.output_folder,
        &mut dfngen.accepted_poly,
        &mut dfngen.intpts,
        &dfngen.triple_points,
        &mut dfngen.pstats,
        &final_fractures,
        &frac_fam_col,
    )?;

    // Duplicate node counters are set in write_output(), so this must come after it.
    info!(
        "Lagrit Should Remove {} Nodes ({}/2 - {})",
        dfngen.pstats.intersection_node_count / 2 - dfngen.pstats.triple_node_count,
        dfngen.pstats.intersection_node_count,
        dfngen.pstats.triple_node_count
    );

    info!("dfngen-rs - Complete");

    Ok(())
}

/// Logs to stdout (compact) and to the given file (JSON).
fn init_tracing(log_file: File) {
    tracing_subscriber::registry()
        .with(tracing_subscriber::EnvFilter::from_default_env())
        .with(
            tracing_subscriber::fmt::layer()
                .compact()
                .with_writer(std::io::stdout),
        )
        .with(
            tracing_subscriber::fmt::layer()
                .json()
                .with_writer(log_file),
        )
        .init();
}

fn domain_volume(input: &Input) -> f64 {
    input.domainSize[0] * input.domainSize[1] * input.domainSize[2]
}

/// Pre-generates the sorted radii list of each family and initializes the
/// per-family counters in `pstats`.
fn prepare_radii_lists(
    input: &Input,
    frac_fam_col: &mut FractureFamilyCollection,
    pstats: &mut Stats,
    generator: Rc<RefCell<Mt64>>,
) -> Result<(), DfngenError> {
    let total_families = frac_fam_col.families.len();

    if input.stopCondition == 0 {
        // Npoly Option
        // Estimate fractures, generate radii lists for nPoly option
        info!("Building radii lists for nPoly option...");
        frac_fam_col.generate_radii(input.forceLargeFractures, input.nPoly, generator.clone())?;
        info!("Building radii lists for nPoly option Complete");
    } else if !input.disableFram {
        // P32 Option
        // ESTIMATE # FRACTURES NEEDED
        info!("Estimating number of fractures needed...");
        dry_run(input, frac_fam_col, generator.clone());
    }

    // Add a percentage more radii to each radii
    // list using families' distribution.
    // First arg is percentage, eg: 0.1 will add 10% more fractures
    // to the radii list for each family
    if !input.disableFram {
        frac_fam_col.generate_radii(
            false,
            (input.nPoly as f64 * input.radiiListIncrease).ceil() as usize,
            generator,
        )?;

        for fracture_family in &frac_fam_col.families {
            info!(
                "Estimated {} fractures for {} family {}",
                fracture_family.radii_list.len(),
                fracture_family.shape,
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

    // Init first rejects per insertion attempt counter
    pstats.rejects_per_attempt.push(0);

    Ok(())
}

/// Inserts fractures from the stochastic families until the stop condition
/// (`nPoly` count or P32 target of every family) is met.
fn insert_stochastic_fractures(
    input: &Input,
    output_folder: &str,
    poly_opts: &PolyOptions,
    frac_fam_col: &mut FractureFamilyCollection,
    dfngen: &mut DFNGen,
    generator: Rc<RefCell<Mt64>>,
) -> Result<(), DfngenError> {
    let total_families = frac_fam_col.families.len();
    let dom_vol = domain_volume(input);

    // 'cdf' size will shrink along when used with fracture intensity (P32) option
    let mut cdf = cumsum(&frac_fam_col.probabilities);
    let mut cdf_size = total_families;

    // Initialize uniform distribution on [0,1]
    let uniform_dist = Uniform::new(0., 1.)?;

    let mut radii_all = if input.outputAllRadii {
        let mut file = File::open(format!("{}/radii/radii_All.dat", output_folder))?;
        let _ = file.write_all("Format: xRadius yRadius Family# (-1 = userRectangle, 0 = userEllipse, > 0 is family in order of famProb)\n".as_bytes());
        Some(file)
    } else {
        None
    };

    // Used with stopCondition = 1, P32 option.
    // Values are set to true once the families p32 requirement is met.
    // Once all elements have values all set to true, all families have had their
    // P32 requirement
    let mut p32_status = vec![false; total_families];

    while (input.stopCondition == 0 && dfngen.pstats.accepted_poly_count < input.nPoly)
        || (input.stopCondition == 1 && p32_status.iter().any(|p| !p))
    {
        // cdfIdx holds the index to the CDF array for the current shape family being inserted
        let mut cdf_idx = 0;
        let roll = generator.borrow_mut().sample(uniform_dist);

        let family_index = if input.stopCondition == 0 {
            // nPoly Option
            // Choose a family based purely on famProb probabilities
            index_from_prob(&cdf, roll)
        } else {
            // P32 Option
            // Choose a family based on probabiliyis AND their target p32 completion status.
            // If a family has already met is fracture intinisty req. (p32) don't choose that family anymore
            index_from_prob_and_p32_status(&p32_status, &cdf, roll, total_families, &mut cdf_idx)
        };

        let mut new_poly = frac_fam_col.families[family_index].create_poly(
            input.eps,
            family_index,
            RadiusOption::FromCacheOrRng,
            generator.clone(),
        );

        if let Some(file) = radii_all.as_mut() {
            let _ = writeln!(
                file,
                "{:.8} {:.8} {:.8}",
                new_poly.xradius,
                new_poly.yradius,
                new_poly.family_num + 1
            );
        }

        let mut reject_counter = 0;

        // Re-translate the polygon until it is accepted or reaches its reject limit
        loop {
            // Truncate poly if needed
            if new_poly.domain_truncation(input.h, input.eps, &input.domainSize) {
                // Poly was completely outside domain, or was truncated to less than
                // 3 vertices due to vertices being too close together
                dfngen.pstats.rejection_reasons.outside += 1;

                if reject_counter >= input.rejectsPerFracture {
                    // Reject poly, generate new polygon
                    break;
                }

                // Retranslate poly and try again, preserving normal, size, and shape
                re_translate_poly(
                    input.eps,
                    &mut new_poly,
                    &frac_fam_col.families[family_index],
                    generator.clone(),
                );
                continue;
            }

            new_poly.assign_bounding_box();
            // Find line of intersection and FRAM check
            let reject_code = intersection_checking(
                poly_opts,
                &mut new_poly,
                &mut dfngen.accepted_poly,
                &mut dfngen.intpts,
                &mut dfngen.pstats,
                &mut dfngen.triple_points,
            );

            if reject_code == 0 {
                dfngen.pstats.accepted_poly_count += 1;
                dfngen.pstats.accepted_from_fam[family_index] += 1;
                // Make new rejection counter for next fracture attempt
                dfngen.pstats.rejects_per_attempt.push(0);

                if new_poly.truncated {
                    dfngen.pstats.truncated += 1;
                }

                new_poly.assign_area();

                let family = &mut frac_fam_col.families[family_index];
                family.add_area_to_p32(new_poly.area, dom_vol, &input.layerVol, &input.regionVol);

                // If the last inserted pologon met the p32 reqirement, set that familiy to no longer
                // insert any more fractures. ajust the CDF and familiy probabilites to account for this
                if input.stopCondition == 1 && family.current_p32 >= family.p32_target {
                    p32_status[family_index] = true;
                    info!("P32 For Family {} Completed", family_index + 1);

                    // Remove the completed family's element in 'CDF[]' and 'famProb[]'
                    // Distribute the removed family probability evenly among the others and rebuild the CDF
                    if cdf_size > 1 {
                        adjust_cdf_and_fam_prob(
                            &mut cdf,
                            &mut frac_fam_col.probabilities,
                            &mut cdf_size,
                            cdf_idx,
                        );
                    }
                }

                if dfngen.pstats.accepted_poly_count.is_multiple_of(200) {
                    log_insertion_progress(input, &frac_fam_col.families, &dfngen.pstats);
                }

                // Intersections and triple points are saved within intersection_checking()
                dfngen.accepted_poly.push(new_poly);
                break;
            }

            reject_counter += 1;
            // Inc reject counter for current attempt
            // (number of rejects until next fracture accepted)
            dfngen.pstats.rejects_per_attempt[dfngen.pstats.accepted_poly_count] += 1;

            if input.printRejectReasons {
                print_reject_reason(reject_code, &new_poly);
            }

            if reject_counter >= input.rejectsPerFracture {
                dfngen.pstats.rejected_poly_count += 1;
                dfngen.pstats.rejected_from_fam[family_index] += 1;
                // Stop retranslating polygon if its reached its reject limit
                break;
            }

            if input.printRejectReasons {
                info!("Translating rejected fracture to new position");
            }

            dfngen.pstats.retranslated_poly_count += 1;
            re_translate_poly(
                input.eps,
                &mut new_poly,
                &frac_fam_col.families[family_index],
                generator.clone(),
            );
        }
    }

    // Remove last element off the rejects per attempt counter.
    // It will have one extra item due to how each element is initialized.
    dfngen.pstats.rejects_per_attempt.pop();

    Ok(())
}

fn log_insertion_progress(input: &Input, families: &[FractureFamily], pstats: &Stats) {
    info!("Accepted {} fractures", pstats.accepted_poly_count);
    info!("Rejected {} fractures", pstats.rejected_poly_count);
    info!("Re-translated {} fractures", pstats.retranslated_poly_count);
    info!("Current p32 values per family:");

    for shape in families {
        if input.stopCondition == 0 {
            info!(
                "{} family {} Current P32 = {:.8}",
                shape.shape, shape.id, shape.current_p32
            );
        } else {
            info!(
                "{} family {} target P32 = {:.8}, Current P32 = {}",
                shape.shape, shape.id, shape.p32_target, shape.current_p32
            );
        }

        if input.stopCondition == 1 && shape.p32_target <= shape.current_p32 {
            info!("...Done");
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

        if poly.family_num >= 0 {
            summary.per_family[poly.family_num as usize] += poly.area;
        } else {
            summary.user_defined += poly.area;
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

/// Logs the layer and the region a family belongs to.
fn log_family_location(frac_fam: &FractureFamily, input: &Input) {
    if frac_fam.layer_id > 0 {
        let idx = (frac_fam.layer_id - 1) * 2;
        info!("Layer: {}", frac_fam.layer_id);
        info!(
            "Layer {{-z, +z}}: {{{}, {}}}",
            input.layers[idx],
            input.layers[idx + 1]
        );
    } else {
        info!("Layer: Whole Domain");
    }

    if frac_fam.region_id > 0 {
        let idx = (frac_fam.region_id - 1) * 6;
        info!("Region: {}", frac_fam.region_id);
        info!("{{-x,+x,-y,+y,-z,+z}}: {:?}", &input.regions[idx..idx + 6]);
    } else {
        info!("Region: Whole Domain");
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
    info!("Truncated: {}", dfngen.pstats.truncated);

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
        let fam_num = dfngen.accepted_poly[i].family_num;
        if fam_num >= 0 {
            final_count_per_family[fam_num as usize] += 1;
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
    info!("If this is the case, try increasing or decreasing the 'radiiListIncrease' option in the input file.");

    for (i, shape) in frac_fam_col.families.iter().enumerate() {
        match shape.radius.function {
            RadiusFunction::Constant(_) => {
                info!("{} Family {} Using constant size", shape.shape, shape.id);
            }
            _ => {
                info!(
                    "{} Family {} Estimated: {}",
                    shape.shape, shape.id, pstats.expected_from_fam[i]
                );

                info!(
                    "Actual: {}",
                    pstats.accepted_from_fam[i] + pstats.rejected_from_fam[i]
                );
            }
        }
    }
}
