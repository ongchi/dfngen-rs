use std::cell::RefCell;
use std::fs::File;
use std::io::Write;
use std::rc::Rc;
use std::time::{SystemTime, UNIX_EPOCH};

use clap::Parser;
use fracture::fracture_family::RadiusOption;
use itertools::zip_eq;
use rand::distr::Uniform;
use rand::Rng;
use rand_mt::Mt64;
use tracing::{error, info, warn};
use tracing_subscriber::{layer::SubscriberExt, util::SubscriberInitExt};

use crate::computational_geometry::{
    create_bounding_box, domain_truncation, intersection_checking, polygon_boundary,
    remove_fractures,
};
use crate::error::DfngenError;
use crate::fracture::cluster_groups::get_cluster;
use crate::fracture::fracture_estimating::dry_run;
use crate::fracture::insert_shape::{get_family_number, print_reject_reason, re_translate_poly};
use crate::fracture::insert_user_defined_fractures;
use crate::io::input::read_input;
use crate::io::output::write_output;
use crate::math_functions::{
    adjust_cdf_and_fam_prob, cumsum, get_area, index_from_prob, index_from_prob_and_p32_status,
};
use crate::structures::{DFNGen, PolyOptions, RadiusFunction, Shape};

mod computational_geometry;
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
    let output_file = File::create(format!("{}/DFN_output.txt", cli.output_folder))?;

    // Setup tracing
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
                .with_writer(output_file),
        )
        .init();

    info!("Starting dfngen-rs");

    /************* Initialize Arrays/Vectors and Structures **************/
    let mut dfngen = DFNGen::new();

    // Read input variables. Most input variables are global
    let (mut input, mut frac_fam_opt) = read_input(&cli.input_file);
    let total_families = frac_fam_opt.families.len();

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
    let dom_vol = input.domainSize[0] * input.domainSize[1] * input.domainSize[2];

    if !frac_fam_opt.families.is_empty() {
        if input.stopCondition == 0 {
            // Npoly Option
            // Estimate fractures, generate radii lists for nPoly option
            info!("Building radii lists for nPoly option...");
            frac_fam_opt.generate_radii(
                input.forceLargeFractures,
                input.nPoly,
                generator.clone(),
            )?;
            info!("Building radii lists for nPoly option Complete");
        } else {
            // P32 Option
            // ESTIMATE # FRACTURES NEEDED
            if !input.disableFram {
                info!("Estimating number of fractures needed...");
                dry_run(&mut input, &mut frac_fam_opt, generator.clone());
            }
        }

        // Add a percentage more radii to each radii
        // list using families' distribution.
        // First arg is percentage, eg: 0.1 will add 10% more fractures
        // to the radii list for each family
        if !input.disableFram {
            frac_fam_opt.generate_radii(
                false,
                (input.nPoly as f64 * input.radiiListIncrease).ceil() as usize,
                generator.clone(),
            )?;

            let mut ell_id = 0;
            let mut rect_id = 0;
            for fracture_family in frac_fam_opt.families.iter() {
                let fam_id = match fracture_family.shape {
                    Shape::Ellipse(_) => {
                        ell_id += 1;
                        ell_id
                    }
                    Shape::Rectangle => {
                        rect_id += 1;
                        rect_id
                    }
                };

                info!(
                    "Estimated {} fractures for {} family {}",
                    fracture_family.radii_list.len(),
                    fracture_family.shape,
                    fam_id
                );
            }
        }

        frac_fam_opt.sort_radii();

        // Keep count of accepted & rejected fractures by family
        dfngen.pstats.accepted_from_fam = vec![0; total_families];
        dfngen.pstats.rejected_from_fam = vec![0; total_families];

        for fracture_family in &frac_fam_opt.families {
            // Save sizes of pre-generated radii lists per family.
            // Print as part of statistics to user
            dfngen
                .pstats
                .expected_from_fam
                .push(fracture_family.radii_list.len());
        }

        // Init first rejects per insertion attempt counter
        dfngen.pstats.rejects_per_attempt.push(0);
    }

    // ********************* User Defined Fractures Insertion ************************

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

    if !frac_fam_opt.families.is_empty() {
        // 'cdf' size will shrink along when used with fracture intensity (P32) option
        let mut cdf = cumsum(&frac_fam_opt.probabilities);
        let mut cdf_size = total_families;

        // Initialize uniform distribution on [0,1]
        let uniform_dist = Uniform::new(0., 1.)?;

        let radii_folder = format!("{}/radii", cli.output_folder);
        let mut radii_all = None;

        if input.outputAllRadii {
            let file = radii_folder + "/radii_All.dat";
            radii_all = Some(File::open(file)?);
            radii_all.as_mut().map(|f| f.write_all("Format: xRadius yRadius Family# (-1 = userRectangle, 0 = userEllipse, > 0 is family in order of famProb)\n".as_bytes()));
        }

        // Used with stopCondition = 1, P32 option.
        // Values are set to true once the families p32 requirement is met.
        // Once all elements have values all set to true, all families have had their
        // P32 requirement
        let mut p32_status = vec![false; total_families];

        // ********* Begin stochastic fracture insertion ***********
        while (input.stopCondition == 0 && dfngen.pstats.accepted_poly_count < input.nPoly)
            || (input.stopCondition == 1 && p32_status.iter().any(|p| !p))
        {
            // cdfIdx holds the index to the CDF array for the current shape family being inserted
            let mut cdf_idx = 0;

            let family_index = if input.stopCondition == 0 {
                // nPoly Option
                // Choose a family based purely on famProb probabilities
                index_from_prob(&cdf, generator.clone().borrow_mut().sample(uniform_dist))
            } else {
                // P32 Option
                // Choose a family based on probabiliyis AND their target p32 completion status.
                // If a family has already met is fracture intinisty req. (p32) don't choose that family anymore
                index_from_prob_and_p32_status(
                    &mut p32_status,
                    &cdf,
                    generator.clone().borrow_mut().sample(uniform_dist),
                    total_families,
                    &mut cdf_idx,
                )
            };

            let mut new_poly = frac_fam_opt.families[family_index].create_poly(
                input.eps,
                input.nFamEll,
                family_index,
                RadiusOption::FromCacheOrRng,
                &input.domainSize,
                &input.domainSizeIncrease,
                &input.layers,
                &input.regions,
                generator.clone(),
            );

            if input.outputAllRadii {
                // Output all radii
                radii_all.as_mut().map(|f| {
                    f.write_all(
                        format!(
                            "{:.8} {:.8} {:.8}\n",
                            &new_poly.xradius,
                            &new_poly.yradius,
                            &new_poly.family_num + 1
                        )
                        .as_bytes(),
                    )
                });
            }

            let mut reject_code = 1;
            let mut reject_counter = 0;

            while reject_code != 0 {
                // Truncate poly if needed
                // 1 if poly is outside of domain or has less than 3 vertices
                if domain_truncation(input.h, input.eps, &mut new_poly, &input.domainSize) {
                    // Poly was completely outside domain, or was truncated to less than
                    // 3 vertices due to vertices being too close together
                    dfngen.pstats.rejection_reasons.outside += 1;

                    // Test if newPoly has reached its limit of insertion attempts
                    if reject_counter >= input.rejectsPerFracture {
                        // Created with new, delete manually
                        new_poly.vertices.clear();
                        // reject_counter += 1;
                        // Reject poly, generate new polygon
                        break;
                    } else {
                        // Retranslate poly and try again, preserving normal, size, and shape
                        re_translate_poly(
                            input.eps,
                            &input.domainSize,
                            &input.domainSizeIncrease,
                            &input.layers,
                            &input.regions,
                            &mut new_poly,
                            &frac_fam_opt.families[family_index],
                            generator.clone(),
                        );
                        continue; // Go to next iteration of while loop, test new translation
                    }
                }

                // Create/assign bounding box
                create_bounding_box(&mut new_poly);
                // Find line of intersection and FRAM check
                // rejectCode = intersectionChecking(newPoly, acceptedPoly, intPts, pstats, triplePoints);
                // Find line of intersection and FRAM check
                reject_code = intersection_checking(
                    input.h,
                    input.eps,
                    input.rFram,
                    input.disableFram,
                    input.tripleIntersections,
                    &mut new_poly,
                    &mut dfngen.accepted_poly,
                    &mut dfngen.intpts,
                    &mut dfngen.pstats,
                    &mut dfngen.triple_points,
                );

                // IF POLY ACCEPTED:
                if reject_code == 0 {
                    // Intersections are ok
                    // Incriment counter of accepted polys
                    dfngen.pstats.accepted_poly_count += 1;
                    dfngen.pstats.accepted_from_fam[family_index] += 1;
                    // Make new rejection counter for next fracture attempt
                    dfngen.pstats.rejects_per_attempt.push(0);

                    if new_poly.truncated {
                        dfngen.pstats.truncated += 1;
                    }

                    // Calculate poly's area
                    new_poly.area = get_area(&new_poly);

                    // Update P32
                    if frac_fam_opt.families[family_index].layer == 0
                        && frac_fam_opt.families[family_index].region == 0
                    {
                        // Whole domain
                        frac_fam_opt.families[family_index].current_p32 +=
                            new_poly.area * 2. / dom_vol;
                    } else if frac_fam_opt.families[family_index].layer > 0
                        && frac_fam_opt.families[family_index].region == 0
                    {
                        // Layer
                        frac_fam_opt.families[family_index].current_p32 += new_poly.area * 2.
                            / input.layerVol[frac_fam_opt.families[family_index].layer - 1];
                    } else if frac_fam_opt.families[family_index].layer == 0
                        && frac_fam_opt.families[family_index].region > 0
                    {
                        // Region
                        frac_fam_opt.families[family_index].current_p32 += new_poly.area * 2.
                            / input.regionVol[frac_fam_opt.families[family_index].region - 1];
                    }

                    if input.stopCondition == 1 {
                        // If the last inserted pologon met the p32 reqirement, set that familiy to no longer
                        // insert any more fractures. ajust the CDF and familiy probabilites to account for this
                        if frac_fam_opt.families[family_index].current_p32
                            >= frac_fam_opt.families[family_index].p32_target
                        {
                            p32_status[family_index] = true; // Mark family as having its p32 requirement met
                            info!("P32 For Family {} Completed", family_index + 1);

                            // Adjust CDF, PDF. Reduce their size by 1.
                            // Remove the completed family's element in 'CDF[]' and 'famProb[]'
                            // Distribute the removed family probability evenly among the others and rebuild the CDF
                            // familyIndex = index of family's probability
                            // cdfIdx = index of the completed family's correspongding CDF index
                            if cdf_size > 1 {
                                // If there are still more families to insert
                                // Remove completed family from CDF and famProb
                                adjust_cdf_and_fam_prob(
                                    &mut cdf,
                                    &mut frac_fam_opt.probabilities,
                                    &mut cdf_size,
                                    cdf_idx,
                                );
                            }
                        }
                    }

                    // Output to user: print running program status to user
                    if dfngen.pstats.accepted_poly_count % 200 == 0 {
                        info!("Accepted {} fractures", dfngen.pstats.accepted_poly_count);
                        info!("Rejected {} fractures", dfngen.pstats.rejected_poly_count);
                        info!(
                            "Re-translated {} fractures",
                            dfngen.pstats.retranslated_poly_count
                        );
                        info!("Current p32 values per family:");

                        for (i, shape) in frac_fam_opt.families.iter().enumerate() {
                            if input.stopCondition == 0 {
                                info!(
                                    "{} family {} Current P32 = {:.8}",
                                    shape.shape,
                                    get_family_number(input.nFamEll, i as isize, shape.shape),
                                    shape.current_p32
                                );
                            } else {
                                info!(
                                    "{} family {} target P32 = {:.8}, Current P32 = {}",
                                    shape.shape,
                                    get_family_number(input.nFamEll, i as isize, shape.shape),
                                    shape.p32_target,
                                    shape.current_p32
                                );
                            }

                            if input.stopCondition == 1 && shape.p32_target <= shape.current_p32 {
                                info!("...Done");
                            }
                        }
                    }

                    // SAVING POLYGON (intersection and triple points saved witchin intersectionChecking())
                    dfngen.accepted_poly.push(new_poly.clone()); // SAVE newPoly to accepted polys list
                } else {
                    // Poly rejected
                    // Inc reject counter for current poly
                    reject_counter += 1;
                    // Inc reject counter for current attempt
                    // (number of rejects until next fracture accepted)
                    dfngen.pstats.rejects_per_attempt[dfngen.pstats.accepted_poly_count] += 1;

                    if input.printRejectReasons {
                        print_reject_reason(reject_code, &new_poly);
                    }

                    if reject_counter >= input.rejectsPerFracture {
                        new_poly.vertices.clear(); // Delete manually, created with new[]
                        dfngen.pstats.rejected_poly_count += 1;
                        dfngen.pstats.rejected_from_fam[family_index] += 1;
                        // Stop retranslating polygon if its reached its reject limit
                        break; // Break will cause code to go to next poly
                    } else {
                        // Translate poly to new position
                        if input.printRejectReasons {
                            info!("Translating rejected fracture to new position");
                        }

                        dfngen.pstats.retranslated_poly_count += 1;
                        re_translate_poly(
                            input.eps,
                            &input.domainSize,
                            &input.domainSizeIncrease,
                            &input.layers,
                            &input.regions,
                            &mut new_poly,
                            &frac_fam_opt.families[family_index],
                            generator.clone(),
                        );
                    }
                } // End else poly rejected
            } // End loop while for re-translating polys option (reject == 1)
        } // !!!!  END MAIN LOOP !!!! end while loop for inserting polyons

        /************************** DFN GENERATION IS COMPLETE ***************************/

        // Remove last element off the rejects per attempt counter.
        // It will have one extra item due to how each element is initialized.
        if !dfngen.pstats.rejects_per_attempt.is_empty() {
            let _ = dfngen.pstats.rejects_per_attempt.pop();
        }
    }

    // The close to node check is inside of the close to edge check function
    // for optimization (only need do intersection close to node on one condition)
    // On close to node rejections, close to edge is counted as well.
    // To get the correct number we must subtract the close to node count
    // (they were counted in closeToEdge AND closeToNode)
    dfngen.pstats.rejection_reasons.close_to_edge -= dfngen.pstats.rejection_reasons.close_to_node;

    // // Assign apertures and permiability to accepted polygons
    // for (unsigned int i = 0; i < acceptedPoly.size(); i++) {
    //     assignAperture(acceptedPoly[i], generator);
    //     // NOTE: must assign aperture before permeability
    //     assignPermeability(acceptedPoly[i]);
    // }
    // Copy end of DFN generation stats to file, as well as print to screen
    info!("Network Generation Complete");
    info!("Version of dfngen-rs: {}", env!("CARGO_PKG_VERSION"));

    if input.stopCondition == 1 {
        info!("Final p32 values per family:");

        for (i, shape) in frac_fam_opt.families.iter().enumerate() {
            info!(
                "Family {} target P32 = {}, Final P32 = {}",
                i + 1,
                shape.p32_target,
                shape.current_p32
            );
            info!(
                "Family {} target P32 = {}, Final P32 = {}",
                i + 1,
                shape.p32_target,
                shape.current_p32,
            );
        }
    }

    // Calculate total area, volume
    let mut user_defined_shapes_area = 0.;
    let mut family_area = vec![0.; total_families];

    info!("Statistics Before Isolated Fractures Removed:",);
    info!("Fractures: {}", dfngen.accepted_poly.len());
    info!("Truncated: {}", dfngen.pstats.truncated);

    // Calculate total fracture area, and area per family
    for i in 0..dfngen.accepted_poly.len() {
        let area = dfngen.accepted_poly[i].area;
        dfngen.pstats.area_before_removal += area;

        if dfngen.accepted_poly[i].family_num >= 0 {
            family_area[dfngen.accepted_poly[i].family_num as usize] += area;
        } else {
            // User-defined polygon
            user_defined_shapes_area += area;
        }
    }

    info!(
        "Total Surface Area: {} m^2",
        dfngen.pstats.area_before_removal * 2.
    );
    info!(
        "Total Fracture Density (P30): {}",
        dfngen.accepted_poly.len() as f64 / dom_vol
    );
    info!(
        "Total Fracture Intensity (P32): {}",
        (dfngen.pstats.area_before_removal * 2.) / dom_vol
    );

    // Print family stats to user
    for (i, (frac_fam, fam_area)) in zip_eq(&frac_fam_opt.families, &family_area).enumerate() {
        info!("Family: {}", i + 1);
        info!("Accepted: {}", dfngen.pstats.accepted_from_fam[i],);
        info!("Rejected: {}", dfngen.pstats.rejected_from_fam[i]);

        if frac_fam.layer > 0 {
            let idx = (frac_fam.layer - 1) * 2;
            info!("Layer: {}", frac_fam.layer,);
            info!(
                "Layer {{-z, +z}}: {{{}, {}}}",
                input.layers[idx],
                input.layers[idx + 1]
            );
        } else {
            info!("Layer: Whole Domain");
        }

        if frac_fam.region > 0 {
            let idx = (frac_fam.region - 1) * 6;
            info!("Region: {}", frac_fam.region);
            info!("{{-x,+x,-y,+y,-z,+z}}: {:?}", &input.regions[idx..idx + 6]);
        } else {
            info!("Region: Whole Domain");
        }

        info!("Surface Area: {} m^2", fam_area * 2.,);
        info!("Fracture Intensity (P32): {}", frac_fam.current_p32);
    }

    if user_defined_shapes_area > 0. {
        info!("User Defined:");
        info!("Surface Area: {} m^2", user_defined_shapes_area * 2.,);
        info!(
            "Fracture Intensity (P32): {}",
            user_defined_shapes_area * 2. / dom_vol
        );
    }

    if input.removeFracturesLessThan > 0. {
        info!(
            "Removing fractures with radius less than {} and rebuilding DFN",
            input.removeFracturesLessThan
        );
        let size = dfngen.accepted_poly.len();
        remove_fractures(input.removeFracturesLessThan, &poly_opts, &mut dfngen);
        info!(
            "Removed {} fractures with radius less than {}",
            size - dfngen.accepted_poly.len(),
            input.removeFracturesLessThan
        );
    }

    if input.polygonBoundaryFlag {
        info!("Extracting fractures from a polygon boundary domain");
        let size = dfngen.accepted_poly.len();
        polygon_boundary(
            &poly_opts,
            input.numOfDomainVertices,
            &input.domainVertices,
            &mut dfngen,
        );
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

    // ************************* Print Statistics to User ***********************************
    info!("Statistics After Isolated Fractures Removed:");
    info!("Final Number of Fractures: {}", final_fractures.len());
    info!(
        "Isolated Fractures Removed: {}",
        dfngen.accepted_poly.len() - final_fractures.len()
    );
    info!(
        "Fractures before isolated fractures removed:: {}",
        dfngen.accepted_poly.len()
    );
    // Reset totalArea to 0
    user_defined_shapes_area = 0.;

    //zero out array.
    for area in family_area.iter_mut() {
        *area = 0.;
    }

    // Calculate total fracture area, and area per family
    for i in 0..final_fractures.len() {
        let area = dfngen.accepted_poly[final_fractures[i]].area;
        dfngen.pstats.area_after_removal += area;

        if dfngen.accepted_poly[final_fractures[i]].family_num >= 0 {
            family_area[dfngen.accepted_poly[final_fractures[i]].family_num as usize] += area;
        } else {
            // User-defined polygon
            user_defined_shapes_area += area;
        }
    }

    // Re-count number of accepted fracture per family after isloated fractures were removed
    let mut accepted_from_fam_counters = vec![0; total_families];

    if total_families > 0 {
        for i in &final_fractures {
            let fam_num = dfngen.accepted_poly[*i].family_num;
            if fam_num >= 0 {
                accepted_from_fam_counters[fam_num as usize] += 1;
            }
        }
    }

    info!(
        "Total Surface Area: {} m^2",
        dfngen.pstats.area_after_removal * 2.
    );
    info!(
        "Total Fracture Density (P30): {}",
        final_fractures.len() as f64 / dom_vol
    );
    info!(
        "Total Fracture Intensity (P32): {}",
        (dfngen.pstats.area_after_removal * 2.) / dom_vol
    );

    // Print family stats to user
    for i in 0..total_families {
        info!("Family: {}", i + 1);
        info!(
            "Fractures After Isolated Fracture Removal: {}",
            accepted_from_fam_counters[i]
        );
        info!(
            "Isolated Fractures Removed: {}",
            dfngen.pstats.accepted_from_fam[i] - accepted_from_fam_counters[i]
        );
        info!("Accepted: {}", dfngen.pstats.accepted_from_fam[i]);
        info!("Rejected: {}", dfngen.pstats.rejected_from_fam[i]);

        if frac_fam_opt.families[i].layer > 0 {
            let idx = (frac_fam_opt.families[i].layer - 1) * 2;
            info!("Layer: {}", frac_fam_opt.families[i].layer);
            info!(
                "Layer {{-z, +z}}: {{{}, {}}}",
                input.layers[idx],
                input.layers[idx + 1]
            );
        } else {
            info!("Layer: Whole Domain");
        }

        if frac_fam_opt.families[i].region > 0 {
            let idx = (frac_fam_opt.families[i].region - 1) * 6;
            info!("Region: {}", frac_fam_opt.families[i].region);
            info!("{{-x,+x,-y,+y,-z,+z}}: {:?}", &input.regions[idx..idx + 6]);
        } else {
            info!("Region: Whole Domain");
        }

        info!("Surface Area: {} m^2", family_area[i] * 2.);
        info!(
            "Fracture Intensity (P32): {}",
            family_area[i] * 2. / dom_vol
        );
    }

    if user_defined_shapes_area > 0. {
        info!("User Defined Shapes:");
        info!("Surface Area: {} m^2", user_defined_shapes_area * 2.);
        info!(
            "Fracture Intensity (P32): {}",
            user_defined_shapes_area * 2. / dom_vol
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

    if print_connectivity_error {
        error!("DFN generation has finished but the formed",);
        error!("fracture network does not make a connection between",);
        error!("the user's specified boundary faces.");
        error!("Try increasing the fracture density, shrinking the domain",);
        error!("or consider using the 'ignoreBoundaryFaces' option.",);
        panic!();
    }

    //************ Intersection Stats ***************
    info!(
        "Number of Triple Intersection Points (Before Isolated Fracture Removal): {}",
        dfngen.triple_points.len()
    );
    // Shrink intersection stats
    info!("Intersection Statistics:");
    info!("Number of Intersections: {}", dfngen.intpts.len());
    info!(
        "Intersections Shortened: {}",
        dfngen.pstats.intersections_shortened
    );
    info!(
        "Original Intersection (Before Intersection Shrinking) Length: {} m",
        dfngen.pstats.original_length
    );
    info!(
        "Intersection Length Discarded: {} m",
        dfngen.pstats.discarded_length
    );
    info!(
        "Final Intersection Length: {} m",
        dfngen.pstats.original_length - dfngen.pstats.discarded_length
    );
    // *********** Rejection Stats *******************
    info!("Rejection Statistics:");
    info!(
        "{} Short Intersections",
        dfngen.pstats.rejection_reasons.short_intersection
    );
    info!(
        "{} Close to Node",
        dfngen.pstats.rejection_reasons.close_to_node
    );
    info!(
        "{} Close to Edge",
        dfngen.pstats.rejection_reasons.close_to_edge
    );
    info!(
        "{} Vertex Close to Edge",
        dfngen.pstats.rejection_reasons.close_point_to_edge
    );
    info!(
        "{} Outside of Domain",
        dfngen.pstats.rejection_reasons.outside
    );
    info!(
        "{} Triple intersection Rejections",
        dfngen.pstats.rejection_reasons.triple
    );
    info!(
        "{} Intersections Close to Other Intersections",
        dfngen.pstats.rejection_reasons.inter_close_to_inter
    );

    if total_families > 0 {
        info!("Fracture Estimation statistics:");
        info!("NOTE: If estimation and actual are very different, expected family distributions might not be accurate.");
        info!("If this is the case, try increasing or decreasing the 'radiiListIncrease' option in the input file.");

        // Compare expected radii/poly size and actual
        for (i, shape) in frac_fam_opt.families.iter().enumerate() {
            match shape.radius.function {
                RadiusFunction::Constant(_) => {
                    info!(
                        "{} Family {} Using constant size",
                        shape.shape,
                        get_family_number(input.nFamEll, i as isize, shape.shape)
                    );
                }
                _ => {
                    info!(
                        "{} Family {} Estimated: {}",
                        shape.shape,
                        get_family_number(input.nFamEll, i as isize, shape.shape),
                        dfngen.pstats.expected_from_fam[i]
                    );

                    info!(
                        "Actual: {}",
                        dfngen.pstats.accepted_from_fam[i] + dfngen.pstats.rejected_from_fam[i]
                    );
                }
            }
        }
    }

    info!("Seed: {}", input.seed);

    // Write all output files
    write_output(
        &input,
        &cli.output_folder,
        &mut dfngen.accepted_poly,
        &mut dfngen.intpts,
        &mut dfngen.triple_points,
        &mut dfngen.pstats,
        &final_fractures,
        &frac_fam_opt,
    );

    // Duplicate node counters are set in writeOutput(). Write output must happen before
    // duplicate node prints
    // Print number of duplicate nodes (pstats.intersectionsNodeCount is set in writeOutpu() )
    info!(
        "Lagrit Should Remove {} Nodes ({}/2 - {})",
        dfngen.pstats.intersection_node_count / 2 - dfngen.pstats.triple_node_count,
        dfngen.pstats.intersection_node_count,
        dfngen.pstats.triple_node_count
    );

    info!("dfngen-rs - Complete");

    Ok(())
}
