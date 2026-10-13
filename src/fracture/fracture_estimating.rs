use rand::{distr::Uniform, Rng};
use rand_mt::Mt64;
use tracing::info;

use crate::error::DfngenError;
use crate::fracture::family_probability::{
    cumulative_sum, index_from_probability_and_p32_status, remove_family_probability,
};
use crate::fracture::fracture_family::{FractureFamilyCollection, RadiusOption};
use crate::io::input::Input;

/// Estimates the number of fractures of each family for the P32 stop condition.
///
/// The function inserts fractures with no FRAM check until each family has its
/// target P32. The radius of each fracture goes into the radii list of its
/// family, for the real generation.
///
/// # Arguments
///
/// * `input` - Input options
/// * `family_collection` - Stochastic fracture families
/// * `generator` - Random number generator
pub fn dry_run(
    input: &Input,
    family_collection: &mut FractureFamilyCollection,
    generator: &mut Mt64,
) -> Result<(), DfngenError> {
    info!("Estimating number of fractures per family for defined fracture intensities (P32)...");
    let domain_volume = input.domain_volume;
    // Probabilities of the families which have not met their P32 target.
    // This is a copy: the family probabilities must stay in their original state.
    let mut probabilities = family_collection.probabilities.clone();
    let uniform_dist = Uniform::new(0., 1.)?;
    let mut cdf = cumulative_sum(&probabilities);
    let mut family_index; // Index of the family of the fracture being generated
    let mut forced_large_fracture_count = 0;
    let mut p32_status = vec![false; family_collection.families.len()];

    'insertion: while p32_status.iter().any(|p| !p) {
        // Index to CDF array of current family being inserted
        let cdf_index;
        let mut reject_counter = 0;
        let mut new_poly = if (forced_large_fracture_count < family_collection.families.len())
            && input.force_large_fractures
        {
            // The CDF has no elements for completed families
            cdf_index = (0..forced_large_fracture_count)
                .filter(|i| !p32_status[*i])
                .count();

            family_index = forced_large_fracture_count;
            forced_large_fracture_count += 1;

            family_collection.families[family_index].create_poly(
                input.eps,
                family_index,
                RadiusOption::MaxRadius,
                generator,
            )?
        } else {
            // Choose a family based on probabiliyis AND their target p32 completion status
            // if a family has already met is fracture intinisty reqirement (p32) dont choose that family anymore
            (family_index, cdf_index) = index_from_probability_and_p32_status(
                &p32_status,
                &cdf,
                generator.sample(uniform_dist),
            );

            family_collection.families[family_index].create_poly(
                input.eps,
                family_index,
                RadiusOption::FromRng,
                generator,
            )?
        };

        while new_poly.is_outside_box_domain(input.min_feature_size, input.eps, &input.domain_size)
        {
            // Poly is completely outside domain
            reject_counter += 1; // Counter for re-trying a new translation

            // Test if the new fracture has reached its limit of insertion attempts
            if reject_counter >= input.rejects_per_fracture {
                // Reject poly, generate new polygon
                continue 'insertion;
            }

            // Retranslate poly and try again, preserving normal, size, and family
            family_collection.families[family_index].retranslate(&mut new_poly, generator);
        }

        // Calculate poly's area
        new_poly.assign_area();

        family_collection.families[family_index].add_area_to_p32(
            new_poly.area,
            domain_volume,
            &input.region_volumes,
        );

        // Save radius for real DFN generation
        family_collection.families[family_index]
            .radii_list
            .push(new_poly.x_radius);

        // If the last inserted polygon met the p32 requirement, set that family to no longer
        // insert any more fractures. adjust the CDF and family probabilities
        if family_collection.families[family_index].current_p32
            >= family_collection.families[family_index].p32_target
        {
            p32_status[family_index] = true; //mark family as having its p32 requirement met

            // Remove the completed family from the probabilities and rebuild the CDF,
            // if there are still more families to insert
            if probabilities.len() > 1 {
                cdf = remove_family_probability(&mut probabilities, cdf_index);
            }
        }

        // No need to save any polygons,
        // We are just simulating dfn with no rejections
        // to get an idea of how many fractures we will
        // need for each family
    } // End while loop for inserting polyons

    // Reset p32 to 0
    for family in &mut family_collection.families {
        family.current_p32 = 0.;
    }

    Ok(())
}
