use std::{cell::RefCell, fs::File, io::Write, path::Path, rc::Rc, time::SystemTime};

use clap::Parser;
use console::Term;
use rand::{distributions::Uniform, Rng};
use rand_mt::Mt19937GenRand64;

use crate::{
    cluster_groups::get_cluster,
    computational_geometry::{create_bounding_box, intersection_checking},
    distributions::Distributions,
    domain::domain_truncation,
    fracture_estimating::{
        add_radii_to_lists, dry_run, generate_radii_lists_n_poly_option, sort_radii,
    },
    insert_shape::{
        generate_poly, get_family_number, p32_complete, print_reject_reason, re_translate_poly,
        shape_type,
    },
    insert_user_ell::insert_user_ell,
    insert_user_ell_by_coord::insert_user_ell_by_coord,
    insert_user_polygon_by_coord::insert_user_polygon_by_coord,
    insert_user_rects::insert_user_rects,
    insert_user_rects_by_coord::insert_user_rects_by_coord,
    math_functions::{
        adjust_cdf_and_fam_prob, create_cdf, get_area, index_from_prob,
        index_from_prob_and_p32_status,
    },
    output::write_output,
    polygon_boundary::polygon_boundary,
    read_input::get_input,
    read_input_functions::get_time_based_seed,
    remove_fractures::remove_fractures,
    structures::{IntPoints, Point, Poly, Shape, Stats},
};

mod cluster_groups;
mod computational_geometry;
mod distributions;
mod domain;
mod exp_dist;
mod fracture_estimating;
mod generating_points;
mod insert_shape;
mod insert_user_ell;
mod insert_user_ell_by_coord;
mod insert_user_polygon_by_coord;
mod insert_user_rects;
mod insert_user_rects_by_coord;
mod math_functions;
mod output;
mod polygon_boundary;
mod read_input;
mod read_input_functions;
mod remove_fractures;
mod structures;
mod vector_functions;

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    /// Path to input file
    input_file: String,

    /// Path to output folder
    output_folder: String,
}

fn main() {
    let cli = Cli::parse();

    println!("Starting DFNGen");

    /************* Initialize Arrays/Vectors and Structures **************/
    // Vector to store accepted polygons/fractures
    let mut accepted_poly: Vec<Poly> = Vec::new();
    // Vector for storing intersections
    let mut int_pts: Vec<IntPoints> = Vec::new();
    // Vector for storing triple intersection points
    let mut triple_points: Vec<Point> = Vec::new();
    // Vector for shape families/ stochastic families
    let mut shape_families: Vec<Shape> = Vec::new();
    // Statistics structure:
    // Keeps track of DFN statistics (see definition in structures.h)
    let mut pstats = Stats::new();
    // *********************************************************************
    // Read Input File
    // Initialize input variables. Most input variables are global
    let mut input = get_input(&cli.input_file, &mut shape_families);
    // Set epsilon
    input.eps = input.h * 1e-8;
    println!("h: {}", input.h);
    let total_families = input.nFamEll + input.nFamRect;

    // Print shape families to screen
    // if totalFamilies > 0 {
    //     printShapeFams(shapeFamilies);
    // }

    // Initialize random generator with seed ( see c++ <random> )
    // Mersene Twister 19937 generator (64 bit)
    if input.seed == 0 {
        input.seed = get_time_based_seed();
    }

    let generator = Rc::new(RefCell::new(Mt19937GenRand64::new(input.seed)));
    let distributions = Rc::new(RefCell::new(Distributions::new(
        &input,
        generator.clone(),
        &mut shape_families,
    )));
    let dom_vol = input.domainSize[0] * input.domainSize[1] * input.domainSize[2];

    if total_families > 0 {
        if input.stopCondition == 0 {
            // Npoly Option
            // Estimate fractures, generate radii lists for nPoly option
            generate_radii_lists_n_poly_option(
                &input,
                &mut shape_families,
                &input.famProb,
                generator.clone(),
                distributions.clone(),
            );
        } else {
            // P32 Option
            // ESTIMATE # FRACTURES NEEDED
            if !input.disableFram {
                println!("Estimating number of fractures needed...");
                dry_run(
                    &mut input,
                    &mut shape_families,
                    generator.clone(),
                    distributions.clone(),
                );
            }
        }

        // Add a percentage more radii to each radii
        // list using families' distribution.
        // First arg is percentage, eg: 0.1 will add 10% more fractures
        // to the radii list for each family
        if !input.disableFram {
            add_radii_to_lists(
                &input,
                input.radiiListIncrease,
                &mut shape_families,
                generator.clone(),
                distributions.clone(),
            );

            for (j, shape) in shape_families.iter().enumerate() {
                if shape.distribution_type == 4 {
                    // Constant size
                    println!(
                        "{} family {} using constant size",
                        shape_type(shape),
                        get_family_number(&input, j as isize, shape.shape_family)
                    );
                } else {
                    println!(
                        "Estimated {} fractures for {} family {}",
                        shape.radii_list.len(),
                        shape_type(shape),
                        get_family_number(&input, j as isize, shape.shape_family)
                    );
                }
            }

            sort_radii(&mut shape_families);
        }

        // Keep count of accepted & rejected fractures by family
        pstats.accepted_from_fam.reserve(total_families);
        pstats.rejected_from_fam.reserve(total_families);
        // Save sizes of pre-generated radii lists per family.
        // Print as part of statistics to user
        pstats.expected_from_fam.reserve(total_families);

        // Zero arrays, init expectedFromFam array
        for (i, shape) in shape_families.iter().enumerate().take(total_families) {
            pstats.accepted_from_fam[i] = 0;
            pstats.rejected_from_fam[i] = 0;
            pstats.expected_from_fam[i] = shape.radii_list.len();
        }

        // Init first rejects per insertion attempt counter
        pstats.rejects_per_attempt.push(0);
    }

    // *********** SETUP HOT KEY *************
    let stdout = Term::buffered_stdout();
    let mut key = '\0';
    #[cfg(feature = "testing")]
    {
        // Set custom terminal settings for hotkey functionality
        set_conio_terminal_mode();
        atexit(reset_terminal_mode);
    }
    // *********  END SETUP HOT KEY **********

    // ********************* User Defined Shapes Insertion ************************
    // User Polygons are always inserted first
    if input.userPolygonByCoord {
        insert_user_polygon_by_coord(
            &input,
            &mut accepted_poly,
            &mut int_pts,
            &mut pstats,
            &mut triple_points,
        );
    }

    if input.insertUserRectanglesFirst {
        // Insert user rects first
        if input.userRectanglesOnOff {
            insert_user_rects(
                &mut input,
                &mut accepted_poly,
                &mut int_pts,
                &mut pstats,
                &mut triple_points,
            );
        }

        // Insert all user rectangles by coordinates
        if input.userRecByCoord {
            insert_user_rects_by_coord(
                &mut input,
                &mut accepted_poly,
                &mut int_pts,
                &mut pstats,
                &mut triple_points,
            );
        }

        // Insert all user ellipses
        if input.userEllipsesOnOff {
            insert_user_ell(
                &mut input,
                &mut accepted_poly,
                &mut int_pts,
                &mut pstats,
                &mut triple_points,
            );
        }

        // Insert all user ellipses by coordinates
        if input.userEllByCoord {
            insert_user_ell_by_coord(
                &mut input,
                &mut accepted_poly,
                &mut int_pts,
                &mut pstats,
                &mut triple_points,
            );
        }
    } else {
        // Insert all user ellipses first
        if input.userEllipsesOnOff {
            insert_user_ell(
                &mut input,
                &mut accepted_poly,
                &mut int_pts,
                &mut pstats,
                &mut triple_points,
            );
        }

        // Insert all user ellipses by coordinates
        if input.userEllByCoord {
            insert_user_ell_by_coord(
                &mut input,
                &mut accepted_poly,
                &mut int_pts,
                &mut pstats,
                &mut triple_points,
            );
        }

        // Insert user rects
        if input.userRectanglesOnOff {
            insert_user_rects(
                &mut input,
                &mut accepted_poly,
                &mut int_pts,
                &mut pstats,
                &mut triple_points,
            );
        }

        // Insert all user rectangles by coordinates
        if input.userRecByCoord {
            insert_user_rects_by_coord(
                &mut input,
                &mut accepted_poly,
                &mut int_pts,
                &mut pstats,
                &mut triple_points,
            );
        }
    }

    /*********  Probabilities (famProb) setup, CDF init  *****************/
    // 'CDF' size will shrink along when used with fracture intensity (P32) option
    let mut cdf = Vec::new();
    let mut cdf_size = total_families;

    if total_families > 0 {
        // Convert famProb to CDF
        cdf = create_cdf(&input.famProb);
    }

    let radii_folder = format!("{}/radii", cli.output_folder);
    let mut radii_all = None;

    if input.outputAllRadii {
        let file = radii_folder + "/radii_All.dat";
        radii_all = Some(File::open(file).unwrap());
        radii_all.as_mut().map(|f| f.write_all("Format: xRadius yRadius Family# (-1 = userRectangle, 0 = userEllipse, > 0 is family in order of famProb)\n".as_bytes()));
    }

    // Initialize uniform distribution on [0,1]
    let uniform_dist = Uniform::new(0., 1.);

    if total_families > 0 {
        // Holds index to current 'shapeFamily' being inserted
        let mut family_index;

        /******************************************************************************************/
        /**************************             MAIN LOOP            ******************************/

        // NOTE: p32Complete() works on global array 'p32Status'
        // p32Complete() only needs argument of the number of defined shape families
        // ********* Begin stochastic fracture insertion ***********
        while ((input.stopCondition == 0 && pstats.accepted_poly_count < input.nPoly)
            || (input.stopCondition == 1 && !p32_complete(&input, total_families)))
            && key != '~'
        {
            // cdfIdx holds the index to the CDF array for the current shape family being inserted
            let mut cdf_idx = 0;

            if input.stopCondition == 0 {
                // nPoly Option
                // Choose a family based purely on famProb probabilities
                family_index =
                    index_from_prob(&cdf, generator.clone().borrow_mut().sample(uniform_dist));
            }
            // Choose a family based on probabiliyis AND their target p32 completion status.
            // If a family has already met is fracture intinisty req. (p32) don't choose that family anymore
            else {
                // P32 Option
                family_index = index_from_prob_and_p32_status(
                    &mut input,
                    &cdf,
                    generator.clone().borrow_mut().sample(uniform_dist),
                    total_families,
                    &mut cdf_idx,
                );
            }

            let mut new_poly = generate_poly(
                &input,
                &mut shape_families[family_index],
                generator.clone(),
                distributions.clone(),
                family_index as isize,
                true,
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
                // Loop used to reinsert same poly with different translation
                // HOT KEY: check for keyboard input
                if let Ok(key_stroke) = stdout.read_char() {
                    key = key_stroke;
                }

                // Truncate poly if needed
                // 1 if poly is outside of domain or has less than 3 vertices
                if domain_truncation(&input, &mut new_poly, &input.domainSize) {
                    // Poly was completely outside domain, or was truncated to less than
                    // 3 vertices due to vertices being too close together
                    pstats.rejection_reasons.outside += 1;

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
                            &input,
                            &mut new_poly,
                            &shape_families[family_index],
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
                    &input,
                    &mut new_poly,
                    &mut accepted_poly,
                    &mut int_pts,
                    &mut pstats,
                    &mut triple_points,
                );

                #[cfg(feature = "testing")]
                if (reject_code != 0) {
                    return 1;
                }

                // IF POLY ACCEPTED:
                if reject_code == 0 {
                    // Intersections are ok
                    // Incriment counter of accepted polys
                    pstats.accepted_poly_count += 1;
                    pstats.accepted_from_fam[family_index] += 1;
                    // Make new rejection counter for next fracture attempt
                    pstats.rejects_per_attempt.push(0);

                    if new_poly.truncated {
                        pstats.truncated += 1;
                    }

                    // Calculate poly's area
                    new_poly.area = get_area(&new_poly);

                    // Update P32
                    if shape_families[family_index].layer == 0
                        && shape_families[family_index].region == 0
                    {
                        // Whole domain
                        shape_families[family_index].current_p32 += new_poly.area * 2. / dom_vol;
                    } else if shape_families[family_index].layer > 0
                        && shape_families[family_index].region == 0
                    {
                        // Layer
                        shape_families[family_index].current_p32 += new_poly.area * 2.
                            / input.layerVol[shape_families[family_index].layer - 1];
                    } else if shape_families[family_index].layer == 0
                        && shape_families[family_index].region > 0
                    {
                        // Region
                        shape_families[family_index].current_p32 += new_poly.area * 2.
                            / input.regionVol[shape_families[family_index].region - 1];
                    }

                    if input.stopCondition == 1 {
                        // If the last inserted pologon met the p32 reqirement, set that familiy to no longer
                        // insert any more fractures. ajust the CDF and familiy probabilites to account for this
                        if shape_families[family_index].current_p32
                            >= shape_families[family_index].p32_target
                        {
                            input.p32Status[family_index] = true; // Mark family as having its p32 requirement met
                            println!("P32 For Family {} Completed", family_index + 1);

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
                                    &mut input.famProb,
                                    &mut cdf_size,
                                    cdf_idx,
                                );
                            }
                        }
                    }

                    // Output to user: print running program status to user
                    if pstats.accepted_poly_count % 200 == 0 {
                        println!("Accepted {} fracturesn", pstats.accepted_poly_count);
                        println!("Rejected {} fracturesn", pstats.rejected_poly_count);
                        println!("Re-translated {} fractures", pstats.retranslated_poly_count);
                        println!("Current p32 values per family:n");

                        for (i, shape) in shape_families.iter().enumerate().take(total_families) {
                            if input.stopCondition == 0 {
                                println!(
                                    "{} family {} Current P32 = {:.8}",
                                    shape_type(shape),
                                    get_family_number(&input, i as isize, shape.shape_family),
                                    shape.current_p32
                                );
                            } else {
                                println!(
                                    "{} family {} target P32 = {:.8}, Current P32 = {}",
                                    shape_type(shape),
                                    get_family_number(&input, i as isize, shape.shape_family),
                                    shape.p32_target,
                                    shape.current_p32
                                );
                            }

                            if input.stopCondition == 1 && shape.p32_target <= shape.current_p32 {
                                println!("...Done");
                            } else {
                                println!();
                            }
                        }
                    }

                    // SAVING POLYGON (intersection and triple points saved witchin intersectionChecking())
                    accepted_poly.push(new_poly.clone()); // SAVE newPoly to accepted polys list
                } else {
                    // Poly rejected
                    // Inc reject counter for current poly
                    reject_counter += 1;
                    // Inc reject counter for current attempt
                    // (number of rejects until next fracture accepted)
                    pstats.rejects_per_attempt[pstats.accepted_poly_count] += 1;

                    if input.printRejectReasons {
                        print_reject_reason(reject_code, &new_poly);
                    }

                    if reject_counter >= input.rejectsPerFracture {
                        new_poly.vertices.clear(); // Delete manually, created with new[]
                        pstats.rejected_poly_count += 1;
                        pstats.rejected_from_fam[family_index] += 1;
                        // Stop retranslating polygon if its reached its reject limit
                        break; // Break will cause code to go to next poly
                    } else {
                        // Translate poly to new position
                        if input.printRejectReasons {
                            println!("Translating rejected fracture to new position");
                        }

                        pstats.retranslated_poly_count += 1;
                        re_translate_poly(
                            &input,
                            &mut new_poly,
                            &shape_families[family_index],
                            generator.clone(),
                        );
                    }
                } // End else poly rejected
            } // End loop while for re-translating polys option (reject == 1)
        } // !!!!  END MAIN LOOP !!!! end while loop for inserting polyons

        /************************** DFN GENERATION IS COMPLETE ***************************/

        // Remove last element off the rejects per attempt counter.
        // It will have one extra item due to how eachelement is initialized.
        if !pstats.rejects_per_attempt.is_empty() {
            let _ = pstats.rejects_per_attempt.pop();
        }
    } // End if totalFamilies != 0

    // The close to node check is inside of the close to edge check function
    // for optimization (only need do intersection close to node on one condition)
    // On close to node rejections, close to edge is counted as well.
    // To get the correct number we must subtract the close to node count
    // (they were counted in closeToEdge AND closeToNode)
    pstats.rejection_reasons.close_to_edge -= pstats.rejection_reasons.close_to_node;

    #[cfg(feature = "testing")]
    reset_terminal_mode();

    if input.outputAllRadii {
        // drop(radiiAll);
    }

    // // Assign apertures and permiability to accepted polygons
    // for (unsigned int i = 0; i < acceptedPoly.size(); i++) {
    //     assignAperture(acceptedPoly[i], generator);
    //     // NOTE: must assign aperture before permeability
    //     assignPermeability(acceptedPoly[i]);
    // }
    // Copy end of DFN generation stats to file, as well as print to screen
    let file_name = format!("{}/DFN_output.txt", cli.output_folder);
    let file_path = Path::new(&file_name);
    let _ = std::fs::create_dir_all(file_path.parent().unwrap());
    let mut file = File::create(file_name).unwrap();
    let now = SystemTime::now();
    log_msg(
        &mut file,
        "\n========================================================\n",
    );
    log_msg(&mut file, "            Network Generation Complete\n");
    log_msg(
        &mut file,
        "========================================================\n",
    );
    log_msg(&mut file, "Version of DFNGen: 2.2\n");
    log_msg(&mut file, &format!("Time Stamp: {:?}\n", &now));

    if input.stopCondition == 1 {
        println!("Final p32 values per family:");

        for (i, shape) in shape_families.iter().enumerate().take(total_families) {
            println!(
                "Family {} target P32 = {}, Final P32 = {}",
                i + 1,
                shape.p32_target,
                shape.current_p32
            );
            println!(
                "Family {} target P32 = {}, Final P32 = {}",
                i + 1,
                shape.p32_target,
                shape.current_p32,
            );
        }
    }

    println!("________________________________________________________");
    log_msg(
        &mut file,
        "\n________________________________________________________\n",
    );
    // Calculate total area, volume
    let mut user_defined_shapes_area = 0.;
    let mut family_area = Vec::new();

    if total_families > 0 {
        family_area = Vec::with_capacity(total_families); // Holds fracture area per family

        // Zero array
        for area in family_area.iter_mut().take(total_families) {
            *area = 0.;
        }
    }

    log_msg(
        &mut file,
        "\nStatistics Before Isolated Fractures Removed:\n\n",
    );
    log_msg(&mut file, &format!("Fractures: {}\n", accepted_poly.len()));
    log_msg(&mut file, &format!("Truncated: {}\n\n", pstats.truncated));

    // Calculate total fracture area, and area per family
    for i in 0..accepted_poly.len() {
        let area = accepted_poly[i].area;
        pstats.area_before_removal += area;

        if accepted_poly[i].family_num >= 0 {
            family_area[accepted_poly[i].family_num as usize] += area;
        } else {
            // User-defined polygon
            user_defined_shapes_area += area;
        }
    }

    log_msg(
        &mut file,
        &format!(
            "Total Surface Area:     {} m^2\n",
            pstats.area_before_removal * 2.
        ),
    );
    log_msg(
        &mut file,
        &format!(
            "Total Fracture Density   (P30): {}\n",
            accepted_poly.len() as f64 / dom_vol
        ),
    );
    log_msg(
        &mut file,
        &format!(
            "Total Fracture Intensity (P32): {}\n",
            (pstats.area_before_removal * 2.) / dom_vol
        ),
    );

    // Print family stats to user
    for i in 0..total_families {
        log_msg(&mut file, &format!("Family: {}\n", i + 1));
        log_msg(
            &mut file,
            &format!("    Accepted: {}\n", pstats.accepted_from_fam[i]),
        );
        log_msg(
            &mut file,
            &format!("    Rejected: {}\n", pstats.rejected_from_fam[i]),
        );

        if shape_families[i].layer > 0 {
            let idx = (shape_families[i].layer - 1) * 2;
            log_msg(
                &mut file,
                &format!("    Layer: {}\n", shape_families[i].layer),
            );
            log_msg(
                &mut file,
                &format!(
                    "    Layer {{-z, +z}}: {{{}, {}}}\n",
                    input.layers[idx],
                    input.layers[idx + 1]
                ),
            );
        } else {
            log_msg(&mut file, "    Layer: Whole Domain \n");
        }

        if shape_families[i].region > 0 {
            let idx = (shape_families[i].region - 1) * 6;
            log_msg(
                &mut file,
                &format!("    Region: {}\n", shape_families[i].region),
            );
            log_msg(
                &mut file,
                &format!(
                    "    {{-x,+x,-y,+y,-z,+z}}: {:?}",
                    &input.regions[idx..idx + 6]
                ),
            );
        } else {
            log_msg(&mut file, "    Region: Whole Domain");
        }

        log_msg(
            &mut file,
            &format!("    Surface Area: {} m^2\n", family_area[i] * 2.),
        );
        log_msg(
            &mut file,
            &format!(
                "    Fracture Intensity (P32): {}\n\n",
                shape_families[i].current_p32
            ),
        );
    }

    if user_defined_shapes_area > 0. {
        log_msg(&mut file, "User Defined:\n");
        log_msg(
            &mut file,
            &format!("    Surface Area: {} m^2\n", user_defined_shapes_area * 2.),
        );
        log_msg(
            &mut file,
            &format!(
                "    Fracture Intensity (P32): {}\n\n",
                user_defined_shapes_area * 2. / dom_vol
            ),
        );
    }

    if input.removeFracturesLessThan > 0. {
        log_msg(
            &mut file,
            &format!(
                "\nRemoving fractures with radius less than {} and rebuilding DFN\n",
                input.removeFracturesLessThan
            ),
        );
        let size = accepted_poly.len();
        remove_fractures(
            &input,
            input.removeFracturesLessThan,
            &mut accepted_poly,
            &mut int_pts,
            &mut triple_points,
            &mut pstats,
        );
        log_msg(
            &mut file,
            &format!(
                "Removed {} fractures with radius less than {}\n\n",
                size - accepted_poly.len(),
                input.removeFracturesLessThan
            ),
        );
    }

    if input.polygonBoundaryFlag {
        log_msg(
            &mut file,
            "\nExtracting fractures from a polygon boundary domain",
        );
        let size = accepted_poly.len();
        polygon_boundary(
            &input,
            &mut accepted_poly,
            &mut int_pts,
            &mut triple_points,
            &mut pstats,
        );
        log_msg(
            &mut file,
            &format!(
                "Removed {} fractures outside subdomain \n\n",
                size - accepted_poly.len()
            ),
        );
    }

    // Remove any isolated fractures and return
    // a list of polygon indices matching the users
    // boundaryFaces option. If input option
    // keepOnlyLargestCluster = 1, return largest
    // cluster matching users boundaryFaces option
    // If ignoreBoundaryFaces input option is on,
    // DFN will keep all fractures with intersections.
    let mut final_fractures = get_cluster(&input, &pstats);
    // Sort fracture indices to retain order by acceptance
    final_fractures.sort();
    // Error check for no boundary connection
    let mut print_connectivity_error = false;

    if final_fractures.is_empty() && !input.ignoreBoundaryFaces {
        print_connectivity_error = true;
        //if there is no fracture network connected users defined boundary faces
        //switch to ignore boundary faces option with notice to user that there is no connectivity
        final_fractures = get_cluster(&input, &pstats);
        //if still no fractures, there is no fracture network
    }

    if final_fractures.is_empty() {
        log_msg(&mut file, "\nError: DFN Generation has finished, however there are no intersecting fractures. Please adjust input parameters.\n");
        log_msg(
            &mut file,
            "Try increasing the fracture density, or shrinking the domain.\n",
        );
        std::process::exit(0)
    }

    // ************************* Print Statistics to User ***********************************
    log_msg(
        &mut file,
        "\n________________________________________________________\n\n",
    );
    log_msg(&mut file, "Statistics After Isolated Fractures Removed:\n");
    log_msg(
        &mut file,
        &format!("Final Number of Fractures: {}\n", final_fractures.len()),
    );
    log_msg(
        &mut file,
        &format!(
            "Isolated Fractures Removed: {}\n",
            accepted_poly.len() - final_fractures.len()
        ),
    );
    log_msg(
        &mut file,
        &format!(
            "Fractures before isolated fractures removed:: {}\n\n",
            accepted_poly.len()
        ),
    );
    // Reset totalArea to 0
    user_defined_shapes_area = 0.;

    if total_families > 0 {
        //zero out array.
        for area in family_area.iter_mut().take(total_families) {
            *area = 0.;
        }
    }

    // Calculate total fracture area, and area per family
    for i in 0..final_fractures.len() {
        let area = accepted_poly[final_fractures[i]].area;
        pstats.area_after_removal += area;

        if accepted_poly[final_fractures[i]].family_num >= 0 {
            family_area[accepted_poly[final_fractures[i]].family_num as usize] += area;
        } else {
            // User-defined polygon
            user_defined_shapes_area += area;
        }
    }

    // Re-count number of accepted fracture per family after isloated fractures were removed
    let mut accepted_from_fam_counters = Vec::new();

    if total_families > 0 {
        accepted_from_fam_counters = Vec::with_capacity(total_families);

        for counter in accepted_from_fam_counters.iter_mut().take(total_families) {
            // zero counters
            *counter = 0;
        }

        let size = final_fractures.len();

        for i in 0..size {
            if accepted_poly[final_fractures[i]].family_num >= 0 {
                accepted_from_fam_counters
                    [accepted_poly[final_fractures[i]].family_num as usize] += 1;
            }
        }
    }

    log_msg(
        &mut file,
        &format!(
            "Total Surface Area:     {} m^2\n",
            pstats.area_after_removal * 2.
        ),
    );
    log_msg(
        &mut file,
        &format!(
            "Total Fracture Density   (P30): {}\n",
            final_fractures.len() as f64 / dom_vol
        ),
    );
    log_msg(
        &mut file,
        &format!(
            "Total Fracture Intensity (P32): {}\n",
            (pstats.area_after_removal * 2.) / dom_vol
        ),
    );

    // Print family stats to user
    for i in 0..total_families {
        log_msg(&mut file, &format!("Family: {}\n", i + 1));
        log_msg(
            &mut file,
            &format!(
                "    Fractures After Isolated Fracture Removal: {}\n",
                accepted_from_fam_counters[i]
            ),
        );
        log_msg(
            &mut file,
            &format!(
                "    Isolated Fractures Removed: {}\n",
                pstats.accepted_from_fam[i] - accepted_from_fam_counters[i]
            ),
        );
        log_msg(
            &mut file,
            &format!("    Accepted: {}\n", pstats.accepted_from_fam[i]),
        );
        log_msg(
            &mut file,
            &format!("    Rejected: {}\n", pstats.rejected_from_fam[i]),
        );

        if shape_families[i].layer > 0 {
            let idx = (shape_families[i].layer - 1) * 2;
            log_msg(
                &mut file,
                &format!("    Layer: {}\n", shape_families[i].layer),
            );
            log_msg(
                &mut file,
                &format!(
                    "    Layer {{-z, +z}}: {{{}, {}}}\n",
                    input.layers[idx],
                    input.layers[idx + 1]
                ),
            );
        } else {
            log_msg(&mut file, "    Layer: Whole Domain\n");
        }

        if shape_families[i].region > 0 {
            let idx = (shape_families[i].region - 1) * 6;
            log_msg(
                &mut file,
                &format!("    Region: {}\n", shape_families[i].region),
            );
            log_msg(
                &mut file,
                &format!(
                    "    {{-x,+x,-y,+y,-z,+z}}: {:?}\n",
                    &input.regions[idx..idx + 6]
                ),
            );
        } else {
            log_msg(&mut file, "    Region: Whole Domain\n");
        }

        log_msg(
            &mut file,
            &format!("    Surface Area: {} m^2\n", family_area[i] * 2.),
        );
        log_msg(
            &mut file,
            &format!(
                "    Fracture Intensity (P32): {}\n\n",
                family_area[i] * 2. / dom_vol
            ),
        );
    }

    if user_defined_shapes_area > 0. {
        log_msg(&mut file, "User Defined Shapes: \n");
        log_msg(
            &mut file,
            &format!("    Surface Area: {} m^2\n", user_defined_shapes_area * 2.),
        );
        log_msg(
            &mut file,
            &format!(
                "    Fracture Intensity (P32): {}\n\n",
                user_defined_shapes_area * 2. / dom_vol
            ),
        );
    }

    accepted_from_fam_counters.clear();

    log_msg(
        &mut file,
        "\n________________________________________________________\n\n",
    );
    log_msg(
        &mut file,
        &format!(
            "\n{} Fractures Accepted (Before Isolated Fracture Removal)\n",
            accepted_poly.len()
        ),
    );
    log_msg(
        &mut file,
        &format!(
            "{} Final Fractures (After Isolated Fracture Removal)\n\n",
            final_fractures.len()
        ),
    );
    log_msg(
        &mut file,
        &format!("Total Fractures Rejected: {}\n", pstats.rejected_poly_count),
    );
    log_msg(
        &mut file,
        &format!(
            "Total Fractures Re-translated: {}\n",
            pstats.retranslated_poly_count
        ),
    );

    if print_connectivity_error {
        log_msg(
            &mut file,
            "\nERROR: DFN generation has finished but the formed\n",
        );
        log_msg(
            &mut file,
            "fracture network does not make a connection between\n",
        );
        log_msg(&mut file, "the user's specified boundary faces.\n");
        log_msg(
            &mut file,
            "Try increasing the fracture density, shrinking the domain\n",
        );
        log_msg(
            &mut file,
            "or consider using the 'ignoreBoundaryFaces' option.\n",
        );
        panic!();
    }

    //************ Intersection Stats ***************
    log_msg(
        &mut file,
        &format!(
            "\nNumber of Triple Intersection Points (Before Isolated Fracture Removal): {}\n",
            triple_points.len()
        ),
    );
    // Shrink intersection stats
    log_msg(&mut file, "\nIntersection Statistics:\n");
    log_msg(
        &mut file,
        &format!("    Number of Intersections: {} \n", int_pts.len()),
    );
    log_msg(
        &mut file,
        &format!(
            "    Intersections Shortened: {} \n",
            pstats.intersections_shortened
        ),
    );
    log_msg(
        &mut file,
        &format!(
            "    Original Intersection (Before Intersection Shrinking) Length: {} m\n",
            pstats.original_length
        ),
    );
    log_msg(
        &mut file,
        &format!(
            "    Intersection Length Discarded: {} m\n",
            pstats.discarded_length
        ),
    );
    log_msg(
        &mut file,
        &format!(
            "    Final Intersection Length: {} m\n",
            pstats.original_length - pstats.discarded_length
        ),
    );
    // *********** Rejection Stats *******************
    log_msg(&mut file, "\nRejection Statistics: \n");
    log_msg(
        &mut file,
        &format!(
            "    {} Short Intersections \n",
            pstats.rejection_reasons.short_intersection
        ),
    );
    log_msg(
        &mut file,
        &format!(
            "    {} Close to Node\n",
            pstats.rejection_reasons.close_to_node
        ),
    );
    log_msg(
        &mut file,
        &format!(
            "    {} Close to Edge\n",
            pstats.rejection_reasons.close_to_edge
        ),
    );
    log_msg(
        &mut file,
        &format!(
            "    {} Vertex Close to Edge\n",
            pstats.rejection_reasons.close_point_to_edge
        ),
    );
    log_msg(
        &mut file,
        &format!(
            "    {} Outside of Domain\n",
            pstats.rejection_reasons.outside
        ),
    );
    log_msg(
        &mut file,
        &format!(
            "    {} Triple intersection Rejections\n",
            pstats.rejection_reasons.triple
        ),
    );
    log_msg(
        &mut file,
        &format!(
            "    {} Intersections Close to Other Intersections\n\n",
            pstats.rejection_reasons.inter_close_to_inter
        ),
    );
    log_msg(
        &mut file,
        "\n________________________________________________________\n\n",
    );

    if total_families > 0 {
        log_msg(&mut file, "Fracture Estimation statistics:\n");
        log_msg(&mut file, "NOTE: If estimation and actual are very different, \nexpected family distributions might ");
        log_msg(&mut file, "not be accurate. \nIf this is the case, try increasing or decreasing \nthe 'radiiListIncrease' option ");
        log_msg(&mut file, "in the input file.\n\n");

        // Compare expected radii/poly size and actual
        for (i, shape) in shape_families.iter().enumerate().take(total_families) {
            if shape.distribution_type == 4 {
                // Constant
                log_msg(
                    &mut file,
                    &format!(
                        "{} Family {}\nUsing constant size\n\n",
                        shape_type(shape),
                        get_family_number(&input, i as isize, shape.shape_family)
                    ),
                );
            } else {
                log_msg(
                    &mut file,
                    &format!(
                        "{} Family {}\nEstimated: {}\n",
                        shape_type(shape),
                        get_family_number(&input, i as isize, shape.shape_family),
                        pstats.expected_from_fam[i]
                    ),
                );

                log_msg(
                    &mut file,
                    &format!(
                        "Actual:    {}\n\n",
                        pstats.accepted_from_fam[i] + pstats.rejected_from_fam[i]
                    ),
                );
            }
        }

        log_msg(
            &mut file,
            "\n________________________________________________________\n\n",
        );
    }

    log_msg(&mut file, &format!("Seed: {}\n", input.seed));

    // Write all output files
    write_output(
        &input,
        &cli.output_folder,
        &mut accepted_poly,
        &mut int_pts,
        &mut triple_points,
        &mut pstats,
        &final_fractures,
        &shape_families,
    );

    // Duplicate node counters are set in writeOutput(). Write output must happen before
    // duplicate node prints
    // Print number of duplicate nodes (pstats.intersectionsNodeCount is set in writeOutpu() )
    log_msg(
        &mut file,
        &format!(
            "\nLagrit Should Remove {} Nodes ({}/2 - {})\n",
            pstats.intersection_node_count / 2 - pstats.triple_node_count,
            pstats.intersection_node_count,
            pstats.triple_node_count
        ),
    );

    println!("DFNGen - Complete");
}

fn log_msg(file: &mut File, message: &str) {
    print!("{}", message);
    file.write_all(message.as_bytes()).unwrap();
}
