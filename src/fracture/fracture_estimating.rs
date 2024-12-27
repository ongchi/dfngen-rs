use std::{cell::RefCell, rc::Rc};

use rand::{distributions::Uniform, Rng};
use rand_mt::Mt19937GenRand64;

use super::domain::domain_truncation;
use crate::{
    distribution::{TruncExp, TruncLogNormal, TruncPowerLaw},
    fracture::insert_shape::{
        generate_poly, generate_poly_with_radius, get_family_number, get_largest_fracture_radius,
        poly_boundary, re_translate_poly, shape_type,
    },
    io::input::Input,
    math_functions::{
        adjust_cdf_and_fam_prob, create_cdf, get_area, index_from_prob_and_p32_status,
    },
    structures::Shape,
};

// **********************************************************************
// ****************  Sort Families Radii Lists  *************************
// Uses std::sort to sort each family's radii list from largest to smallest.
// This will allow the DFN gereration to start from largest to smallest
// fractures.
// Arg 1: vector<Shape> array of stochastic fracture families */
pub fn sort_radii(shape_fam: &mut Vec<Shape>) {
    for shape in shape_fam {
        shape.radii_list.sort_by(|a, b| b.partial_cmp(a).unwrap())
    }
}

// **********************************************************************
// **********************************************************************
// ***  Create Radii Lists for Shape Families When Using NPoly Option ***
//     Estimates the number of fractures needed for each family and
//     creates radii lists for each family based on their distribution.
//     Arg 1: vector<Shape> array of stochastic fracture families
//     Arg 2: Family probablity array ('famProb' in input file)
//     Arg 3: Random number generator, see std <random> library
//     Arg 4: Reference to Distributions class (used for exponential distribution)
pub fn generate_radii_lists_n_poly_option(
    force_large_fractures: bool,
    n_poly: usize,
    h: f64,
    n_fam_ell: usize,
    shape_families: &mut [Shape],
    fam_prob: &[f64],
    generator: Rc<RefCell<Mt19937GenRand64>>,
) {
    println!("Building radii lists for nPoly option...");

    if force_large_fractures {
        for shape in shape_families.iter_mut() {
            let radius = get_largest_fracture_radius(shape);
            shape.radii_list.push(radius);
        }
    }

    for i in 0..shape_families.len() {
        let amount_to_add = if force_large_fractures {
            (fam_prob[i] * (n_poly - shape_families.len()) as f64).ceil() as usize
        } else {
            (fam_prob[i] * n_poly as f64).ceil() as usize
        };

        add_radii(
            h,
            n_fam_ell,
            amount_to_add,
            i as isize,
            &mut shape_families[i],
            generator.clone(),
        );
    }

    println!("Building radii lists for nPoly option Complete");
}

// **********************************************************************
// ********************  Print Warining to User  ************************
//     This function prints a warning to the user when the random generation
//     of fracture radii lenths is continuously smaller than the minimum defined
//     radii allowed (defined by user in the input file)
//     Arg 1: Index to the family in vecotr<Shape> array the warning is refering to
//     Arg 2: Shape structure the warning is refering to
pub fn print_generating_fractures_less_than_hwarning(
    n_fam_ell: usize,
    fam_index: isize,
    shape_fam: &Shape,
) {
    println!(
        "WARNING: {} Family {} is attepting to populate fracture radii lists, however many fractures are being generated with radii less than 3*h (Minimum radius). Consider adjusting distribution parameters.",
        shape_type(shape_fam),
        get_family_number(n_fam_ell, fam_index, shape_fam.shape_family)
    );
}

// **********************************************************************
// *********  Add Percentage More Radii To Radii Lists  *****************
//     Function adds a percentage more radii to the fracture families
//     radii lists based on each families distribution.
//     This helps account for fracture rejections
//     Arg 1: Percentage to increase the list by. eg .10 will add %10 more radii
//     Arg 2: vector<Shape> array of stochastic fracture families
//     Arg 3: Random number generator (see std <random> library)
//     Arg 4: Distributions class (currently only used for exponential dist)
pub fn add_radii_to_lists(
    input: &Input,
    percent: f64,
    shape_families: &mut [Shape],
    generator: Rc<RefCell<Mt19937GenRand64>>,
) {
    for (i, shape) in shape_families.iter_mut().enumerate() {
        let amount_to_add = (shape.radii_list.len() as f64 * percent).ceil() as usize;
        add_radii(
            input.h,
            input.nFamEll,
            amount_to_add,
            i as isize,
            shape,
            generator.clone(),
        );
    }
}

// **********************************************************************
// ************  Add Radii To Shape Families Radii List  ****************
//     Adds 'amountToAdd' more radii to 'shapeFam's radii list
//     Arg 1: Number of fractures to add to the list
//     Arg 2: Family index to the global Shape structure array ('shapeFamilies' in main())
//            which the radii are being added to
//     Arg 3: The 'Shape' structure which the radii are being added to
//     Arg 4: Random number generator (see std <random> library)
//     Arg 5: Distributions class (currently only used for exponential dist)
pub fn add_radii(
    h: f64,
    n_fam_ell: usize,
    amount_to_add: usize,
    fam_idx: isize,
    shape_fam: &mut Shape,
    generator: Rc<RefCell<Mt19937GenRand64>>,
) {
    let min_radius = 3. * h;

    match shape_fam.distribution_type {
        // Lognormal
        1 => {
            let log_dist =
                TruncLogNormal::new(min_radius, f64::INFINITY, shape_fam.mean, shape_fam.sd)
                    .unwrap();

            for _ in 0..amount_to_add {
                match generator.borrow_mut().sample(&log_dist) {
                    Ok(radius) => shape_fam.radii_list.push(radius),
                    Err(_) => {
                        print_generating_fractures_less_than_hwarning(
                            n_fam_ell, fam_idx, shape_fam,
                        );
                        panic!()
                    }
                };
            }
        }

        // Truncated power-law
        2 => {
            let power_law = TruncPowerLaw::new(
                f64::max(shape_fam.min, min_radius),
                shape_fam.max,
                shape_fam.alpha,
            );

            for _ in 0..amount_to_add {
                let radius = generator.borrow_mut().sample(&power_law);
                shape_fam.radii_list.push(radius);
            }
        }

        // Exponential
        3 => {
            let exp_dist = TruncExp::new(
                min_radius,
                f64::INFINITY,
                shape_fam.exp_lambda,
                shape_fam.min_dist_input,
                shape_fam.max_dist_input,
            )
            .unwrap();

            for _ in 0..amount_to_add {
                match generator.clone().borrow_mut().sample(&exp_dist) {
                    Ok(radius) => shape_fam.radii_list.push(radius),
                    Err(_) => {
                        print_generating_fractures_less_than_hwarning(
                            n_fam_ell, fam_idx, shape_fam,
                        );
                        panic!()
                    }
                };
            }
        }

        _ => unreachable!(),
    }
}

// **********************************************************************
// ********  Estimate Number of Fractures When P32 Option is Used  ******
//     Inserts fractures into domain with FRAM disabled
//     Simply inserts and truncates fractures on the domain
//     until reqired P32 is met.
//     Used to estimate and generate radii lists for each fracture
//     family.
//     Arg 1: vector<Shape> array of stochastic fracture families
//     Arg 2: The probabilities for each families insertion into domain
//            (famProb) in input file
//     Arg 3: Random number generator (see std <random> library)
//     Arg 4: Distributions class (currently only used for exponential dist)
pub fn dry_run(
    input: &mut Input,
    shape_families: &mut [Shape],
    generator: Rc<RefCell<Mt19937GenRand64>>,
) {
    println!("Estimating number of fractures per family for defined fracture intensities (P32)...");
    let dom_vol = input.domainSize[0] * input.domainSize[1] * input.domainSize[2];
    let total_families = shape_families.len();
    let mut cdf_size = total_families; // This variable shrinks along with CDF when used with fracture intensity (P32) option
                                       // Create a copy of the family probablity
                                       // Algoithms used in this function modify this array,
                                       // we need to keep the original in its original state
    let mut fam_probability = input.famProb.clone();
    // std::copy(shapeProb, shapeProb + totalFamilies, famProbability);
    // Init uniform dist on [0,1)
    let uniform_dist = Uniform::new(0., 1.);
    /******  Convert famProb to CDF  *****/
    let mut cdf = create_cdf(&fam_probability);
    let mut family_index; // Holds index to shape family of fracture being generated
    let mut force_large_fract_count = 0;

    while input.p32Status.iter().any(|p| !p) {
        // Index to CDF array of current family being inserted
        let mut cdf_idx = 0;
        let mut reject_counter = 0;
        let mut new_poly =
            if (force_large_fract_count < shape_families.len()) && input.forceLargeFractures {
                let radius = get_largest_fracture_radius(&shape_families[force_large_fract_count]);
                family_index = force_large_fract_count;

                // Choose CDF randomly by family
                cdf_idx = (0..force_large_fract_count)
                    .filter(|i| !input.p32Status[*i])
                    .count()
                    - 1;

                force_large_fract_count += 1;
                let boundary = poly_boundary(
                    &input.domainSize,
                    &input.domainSizeIncrease,
                    &input.layers,
                    &input.regions,
                    &shape_families[force_large_fract_count - 1],
                );
                generate_poly_with_radius(
                    input.orientationOption,
                    input.eps,
                    radius,
                    &shape_families[force_large_fract_count - 1],
                    boundary,
                    generator.clone(),
                    family_index as isize,
                )
            } else {
                // Choose a family based on probabiliyis AND their target p32 completion status
                // if a family has already met is fracture intinisty req. (p32) dont choose that family anymore
                // Choose a family based on probabiliyis AND their target p32 completion status
                // if a family has already met is fracture intinisty reqirement (p32) dont choose that family anymore
                family_index = index_from_prob_and_p32_status(
                    &mut input.p32Status,
                    &cdf,
                    generator.clone().borrow_mut().sample(uniform_dist),
                    total_families,
                    &mut cdf_idx,
                );
                generate_poly(
                    input.h,
                    input.nFamEll,
                    &input.domainSize,
                    &input.domainSizeIncrease,
                    &input.layers,
                    &input.regions,
                    input.orientationOption,
                    input.eps,
                    &mut shape_families[family_index],
                    generator.clone(),
                    family_index as isize,
                    false,
                )
            };

        // Truncate poly if needed
        // Returns 1 if poly is outside of domain or truncated to less than 3 vertices
        // Vector for storing intersection boundaries
        let mut reject = false;

        while domain_truncation(input.h, input.eps, &mut new_poly, &input.domainSize) {
            // Poly is completely outside domain, or was truncated to
            // less than 3 vertices due to vertices being too close together
            reject_counter += 1; // Counter for re-trying a new translation

            // Test if newPoly has reached its limit of insertion attempts
            if reject_counter >= input.rejectsPerFracture {
                reject = true;
                break; // Reject poly, generate new polygon
            } else {
                // Retranslate poly and try again, preserving normal, size, and shape
                re_translate_poly(
                    input,
                    &mut new_poly,
                    &shape_families[family_index],
                    generator.clone(),
                );
            }
        }

        if reject {
            // Restart while loop
            // Generate new fracture
            continue;
        }

        // Calculate poly's area
        new_poly.area = get_area(&new_poly);

        // Update P32
        if shape_families[family_index].layer == 0 && shape_families[family_index].region == 0 {
            // Whole domain
            shape_families[family_index].current_p32 += new_poly.area * 2. / dom_vol;
        } else if shape_families[family_index].layer > 0 && shape_families[family_index].region == 0
        {
            // Layer
            shape_families[family_index].current_p32 +=
                new_poly.area * 2. / input.layerVol[shape_families[family_index].layer - 1];
        } else if shape_families[family_index].layer == 0 && shape_families[family_index].region > 0
        {
            // Region
            shape_families[family_index].current_p32 +=
                new_poly.area * 2. / input.regionVol[shape_families[family_index].region - 1];
        }

        // Save radius for real DFN generation
        shape_families[family_index]
            .radii_list
            .push(new_poly.xradius);

        // If the last inserted polygon met the p32 requirement, set that family to no longer
        // insert any more fractures. adjust the CDF and family probabilities
        if shape_families[family_index].current_p32 >= shape_families[family_index].p32_target {
            input.p32Status[family_index] = true; //mark family as having its p32 requirement met

            // Adjust CDF, PDF, and reduce their size by 1. Keep probabilities proportional.
            // Remove the completed families element in the CDF and famProb[]
            // Distribute the removed family probability evenly among the others and rebuild the CDF
            // familyIndex = index of family's probability
            // cdfIdx = index of the family's correspongding cdf, (index to elmt to remove)
            if cdf_size > 1 {
                // If there are still more families to insert ( cdfSize = 0 means no more families to insert)
                adjust_cdf_and_fam_prob(&mut cdf, &mut fam_probability, &mut cdf_size, cdf_idx);
            }
        }

        // No need to save any polygons,
        // We are just simulating dfn with no rejections
        // to get an idea of how many fractures we will
        // need for each family
        new_poly.vertices.clear();
    } // End while loop for inserting polyons

    // Reset p32 to 0
    for (i, shape) in shape_families.iter_mut().enumerate().take(total_families) {
        input.p32Status[i] = false;
        shape.current_p32 = 0.;
    }
}
