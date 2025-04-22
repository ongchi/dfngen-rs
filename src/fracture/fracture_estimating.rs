use std::{cell::RefCell, rc::Rc};

use rand::{distr::Uniform, Rng};
use rand_mt::Mt64;
use tracing::info;

use super::domain::domain_truncation;
use crate::{
    fracture::insert_shape::{
        generate_poly, generate_poly_with_radius, poly_boundary, re_translate_poly,
    },
    io::input::Input,
    math_functions::{adjust_cdf_and_fam_prob, cumsum, get_area, index_from_prob_and_p32_status},
    structures::FractureFamilyOption,
};

/// Estimate Number of Fractures When P32 Option is Used
///
/// Inserts fractures into domain with FRAM disabled
/// Simply inserts and truncates fractures on the domain
/// until reqired P32 is met.
/// Used to estimate and generate radii lists for each fracture
/// family.
///
/// # Arguments
///
/// * `input` - Input structure
/// * `frac_families` - Array of stochastic fracture families
/// * `generator` - Random number generator
pub fn dry_run(
    input: &mut Input,
    frac_families: &mut FractureFamilyOption,
    generator: Rc<RefCell<Mt64>>,
) {
    info!("Estimating number of fractures per family for defined fracture intensities (P32)...");
    let dom_vol = input.domainSize[0] * input.domainSize[1] * input.domainSize[2];
    let total_families = frac_families.families.len();
    let mut cdf_size = total_families; // This variable shrinks along with CDF when used with fracture intensity (P32) option
                                       // Create a copy of the family probablity
                                       // Algoithms used in this function modify this array,
                                       // we need to keep the original in its original state
                                       // std::copy(shapeProb, shapeProb + totalFamilies, famProbability);
                                       // Init uniform dist on [0,1)
    let uniform_dist = Uniform::new(0., 1.).unwrap();
    /******  Convert famProb to CDF  *****/
    let mut cdf = cumsum(&frac_families.probabilities);
    let mut family_index; // Holds index to shape family of fracture being generated
    let mut force_large_fract_count = 0;
    let mut p32_status = vec![false; frac_families.families.len()];

    while p32_status.iter().any(|p| !p) {
        // Index to CDF array of current family being inserted
        let mut cdf_idx = 0;
        let mut reject_counter = 0;
        let mut new_poly = if (force_large_fract_count < frac_families.families.len())
            && input.forceLargeFractures
        {
            let radius = frac_families.families[force_large_fract_count].radius.max;
            family_index = force_large_fract_count;

            // Choose CDF randomly by family
            cdf_idx = (0..force_large_fract_count)
                .filter(|i| !p32_status[*i])
                .count()
                - 1;

            force_large_fract_count += 1;
            let boundary = poly_boundary(
                &input.domainSize,
                &input.domainSizeIncrease,
                &input.layers,
                &input.regions,
                &frac_families.families[force_large_fract_count - 1],
            );
            generate_poly_with_radius(
                input.eps,
                radius,
                &frac_families.families[force_large_fract_count - 1],
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
                &mut p32_status,
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
                input.eps,
                &mut frac_families.families[family_index],
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
                    input.eps,
                    &input.domainSize,
                    &input.domainSizeIncrease,
                    &input.layers,
                    &input.regions,
                    &mut new_poly,
                    &frac_families.families[family_index],
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
        if frac_families.families[family_index].layer == 0
            && frac_families.families[family_index].region == 0
        {
            // Whole domain
            frac_families.families[family_index].current_p32 += new_poly.area * 2. / dom_vol;
        } else if frac_families.families[family_index].layer > 0
            && frac_families.families[family_index].region == 0
        {
            // Layer
            frac_families.families[family_index].current_p32 +=
                new_poly.area * 2. / input.layerVol[frac_families.families[family_index].layer - 1];
        } else if frac_families.families[family_index].layer == 0
            && frac_families.families[family_index].region > 0
        {
            // Region
            frac_families.families[family_index].current_p32 += new_poly.area * 2.
                / input.regionVol[frac_families.families[family_index].region - 1];
        }

        // Save radius for real DFN generation
        frac_families.families[family_index]
            .radii_list
            .push(new_poly.xradius);

        // If the last inserted polygon met the p32 requirement, set that family to no longer
        // insert any more fractures. adjust the CDF and family probabilities
        if frac_families.families[family_index].current_p32
            >= frac_families.families[family_index].p32_target
        {
            p32_status[family_index] = true; //mark family as having its p32 requirement met

            // Adjust CDF, PDF, and reduce their size by 1. Keep probabilities proportional.
            // Remove the completed families element in the CDF and famProb[]
            // Distribute the removed family probability evenly among the others and rebuild the CDF
            // familyIndex = index of family's probability
            // cdfIdx = index of the family's correspongding cdf, (index to elmt to remove)
            if cdf_size > 1 {
                // If there are still more families to insert ( cdfSize = 0 means no more families to insert)
                adjust_cdf_and_fam_prob(
                    &mut cdf,
                    &mut frac_families.probabilities,
                    &mut cdf_size,
                    cdf_idx,
                );
            }
        }

        // No need to save any polygons,
        // We are just simulating dfn with no rejections
        // to get an idea of how many fractures we will
        // need for each family
        new_poly.vertices.clear();
    } // End while loop for inserting polyons

    // Reset p32 to 0
    for shape in frac_families.families.iter_mut() {
        shape.current_p32 = 0.;
    }
}
