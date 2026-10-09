use std::{cell::RefCell, rc::Rc};

use rand::{distr::Uniform, Rng};
use rand_mt::Mt64;
use tracing::info;

use crate::fracture::fracture_family::{FractureFamilyCollection, RadiusOption};
use crate::fracture::insert_shape::re_translate_poly;
use crate::io::input::Input;
use crate::math_functions::{cumsum, index_from_prob_and_p32_status, remove_family_probability};

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
    input: &Input,
    frac_families: &mut FractureFamilyCollection,
    generator: Rc<RefCell<Mt64>>,
) {
    info!("Estimating number of fractures per family for defined fracture intensities (P32)...");
    let dom_vol = input.domain_size[0] * input.domain_size[1] * input.domain_size[2];
    // Probabilities of the families which have not met their P32 target.
    // This is a copy: the family probabilities must stay in their original state.
    let mut probabilities = frac_families.probabilities.clone();
    let uniform_dist = Uniform::new(0., 1.).unwrap();
    let mut cdf = cumsum(&probabilities);
    let mut family_index; // Holds index to shape family of fracture being generated
    let mut force_large_fract_count = 0;
    let mut p32_status = vec![false; frac_families.families.len()];

    while p32_status.iter().any(|p| !p) {
        // Index to CDF array of current family being inserted
        let cdf_idx;
        let mut reject_counter = 0;
        let mut new_poly = if (force_large_fract_count < frac_families.families.len())
            && input.force_large_fractures
        {
            // The CDF has no elements for completed families
            cdf_idx = (0..force_large_fract_count)
                .filter(|i| !p32_status[*i])
                .count();

            family_index = force_large_fract_count;
            force_large_fract_count += 1;

            frac_families.families[family_index].create_poly(
                input.eps,
                family_index,
                RadiusOption::MaxRadius,
                generator.clone(),
            )
        } else {
            // Choose a family based on probabiliyis AND their target p32 completion status
            // if a family has already met is fracture intinisty reqirement (p32) dont choose that family anymore
            (family_index, cdf_idx) = index_from_prob_and_p32_status(
                &p32_status,
                &cdf,
                generator.clone().borrow_mut().sample(uniform_dist),
            );

            frac_families.families[family_index].create_poly(
                input.eps,
                family_index,
                RadiusOption::FromRng,
                generator.clone(),
            )
        };

        // Truncate poly if needed
        // Returns 1 if poly is outside of domain or truncated to less than 3 vertices
        // Vector for storing intersection boundaries
        let mut reject = false;

        while new_poly.domain_truncation(input.h, input.eps, &input.domain_size) {
            // Poly is completely outside domain, or was truncated to
            // less than 3 vertices due to vertices being too close together
            reject_counter += 1; // Counter for re-trying a new translation

            // Test if newPoly has reached its limit of insertion attempts
            if reject_counter >= input.rejects_per_fracture {
                reject = true;
                break; // Reject poly, generate new polygon
            } else {
                // Retranslate poly and try again, preserving normal, size, and shape
                re_translate_poly(
                    input.eps,
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
        new_poly.assign_area();

        frac_families.families[family_index].add_area_to_p32(
            new_poly.area,
            dom_vol,
            &input.layer_volumes,
            &input.region_volumes,
        );

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

            // Remove the completed family from the probabilities and rebuild the CDF,
            // if there are still more families to insert
            if probabilities.len() > 1 {
                cdf = remove_family_probability(&mut probabilities, cdf_idx);
            }
        }

        // No need to save any polygons,
        // We are just simulating dfn with no rejections
        // to get an idea of how many fractures we will
        // need for each family
    } // End while loop for inserting polyons

    // Reset p32 to 0
    for shape in frac_families.families.iter_mut() {
        shape.current_p32 = 0.;
    }
}
