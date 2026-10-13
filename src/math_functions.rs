use tracing::{error, warn};

/// There are several spots in the code which used to use standard
/// deviation to know how to sort a list of points.
/// To increase performance, the standard deviation is now calculated without
/// the square root and division by N elements.
/// Return: The sum-deviation of array 'vec'
pub fn sum_deviation(vec: &[f64]) -> f64 {
    let mean = vec.iter().sum::<f64>() / (vec.len() as f64);
    vec.iter().map(|&v| (v - mean).powi(2)).sum()
}

/// Get Max Element's Index
///
/// Used to find index of element with max value from array4
/// Used in intersectionChecking()
///
/// # Arguments
///
/// * `vec` - Array of 4 elements
///
/// # Returns
///
/// Max element's array index
///
/// If elements happen to be equal, returns the first one
pub fn max_elmt_idx(vec: &[f64]) -> usize {
    vec.iter()
        .enumerate()
        .reduce(|max, el| if el.1 < max.1 { max } else { el })
        .map(|(i, _)| i)
        .unwrap_or(0)
}

/// Sort Array Indices
///
/// Similar to mathematica's Ordering[] funct.
///
/// # Arguments
///
/// * `vec` - Array of f64
///
/// # Returns
///
/// Array of sorted indicies to elements in 'v' sorted smallest to largest
pub fn sorted_index(vec: &[f64]) -> Vec<usize> {
    let mut v = vec.iter().enumerate().collect::<Vec<(usize, &f64)>>();
    v.sort_by(|a, b| (a.1).partial_cmp(b.1).unwrap());
    v.into_iter().map(|(i, _)| i).collect()
}

/// Index from Probability
///
/// CDF is 1 to 1 and algined with the stochastic family
/// shapes array (std vector). This chooses the family index based
/// on a random roll (random number between 0 and 1) and 'famProb'
///
/// # Arguments
///
/// * `cdf` - CDF of shape families based on famProb array in input file
/// * `roll` - Random number between 0 and 1
///
/// # Returns
///
/// Index to family based on the probability famProb
pub fn index_from_prob(cdf: &[f64], roll: f64) -> usize {
    for (i, v) in cdf.iter().enumerate() {
        if &roll <= v {
            return i;
        }
    }

    cdf.len() - 1
}

/// Choose Family Randomly Based On P32 and CDF
///
/// Use with fracture intensity (p32) stopping option
///
/// # Arguments
///
/// * `p32_status` - One element per family, true if the family has met its P32 target
/// * `cdf` - CDF of the families which have not met their P32 target, in family order
/// * `roll` - Random number between 0 and 1
///
/// # Returns
///
/// The index of the chosen family, and its index in `cdf`.
/// A family which has already met its fracture intensity requirement is never chosen.
pub fn index_from_prob_and_p32_status(
    p32_status: &[bool],
    cdf: &[f64],
    roll: f64,
) -> (usize, usize) {
    let cdf_idx = index_from_prob(cdf, roll);

    // The CDF has no elements for completed families, so cdf_idx counts
    // only the families which are not complete
    let family_index = p32_status
        .iter()
        .enumerate()
        .filter(|(_, is_complete)| !**is_complete)
        .nth(cdf_idx)
        .map(|(family_index, _)| family_index);

    match family_index {
        Some(family_index) => (family_index, cdf_idx),
        None => {
            error!("see indexFromProb_and_P32Status(), funct did not return anything");
            (p32_status.len() - 1, cdf_idx)
        }
    }
}

/// Cumulative sum of the elements along a given axis
pub fn cumsum(fam_prob: &[f64]) -> Vec<f64> {
    let cdf: Vec<f64> = fam_prob
        .iter()
        .scan(0., |acc, &x| {
            *acc += x;
            Some(*acc)
        })
        .collect();

    if let Some(last) = cdf.last() {
        if !(0.999..=1.001).contains(last) {
            warn!(
                "Family probabilities (probability of each family in input file) do not sum to 1"
            );
            warn!("sum = {:.17}", last);
            warn!("Please check input file.");
        }
    };

    cdf
}

/// Removes a family from the family probabilities. Used with P32 stopping
/// condition (see input file) when the P32 requirement of a family is met.
///
/// The probability of the removed family is divided equally among the
/// remaining families.
///
/// # Arguments
///
/// * `fam_probability` - Probability of each family which is not complete
/// * `idx2remove` - Index in `fam_probability` of the family to remove
///
/// # Returns
///
/// The CDF of the remaining families
pub fn remove_family_probability(fam_probability: &mut Vec<f64>, idx2remove: usize) -> Vec<f64> {
    let removed_probability = fam_probability.remove(idx2remove);
    let add_to_remaining_elmts = removed_probability / fam_probability.len() as f64;

    for probability in fam_probability.iter_mut() {
        *probability += add_to_remaining_elmts;
    }

    cumsum(fam_probability)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_p32_selection_skips_completed_families() {
        // Families 0 and 2 are complete. The CDF has families 1 and 3.
        let p32_status = [true, false, true, false];
        let cdf = [0.5, 1.0];

        assert_eq!(
            index_from_prob_and_p32_status(&p32_status, &cdf, 0.2),
            (1, 0)
        );
        assert_eq!(
            index_from_prob_and_p32_status(&p32_status, &cdf, 0.9),
            (3, 1)
        );
    }

    #[test]
    fn test_p32_selection_with_no_completed_family_follows_cdf() {
        let p32_status = [false, false, false];
        let cdf = [0.2, 0.5, 1.0];

        assert_eq!(
            index_from_prob_and_p32_status(&p32_status, &cdf, 0.4),
            (1, 1)
        );
    }

    #[test]
    fn test_removed_family_probability_is_shared_equally() {
        let mut probabilities = vec![0.5, 0.3, 0.2];

        let cdf = remove_family_probability(&mut probabilities, 1);

        assert_eq!(probabilities.len(), 2);
        assert!((probabilities[0] - 0.65).abs() < 1e-12);
        assert!((probabilities[1] - 0.35).abs() < 1e-12);
        assert!((cdf[1] - 1.0).abs() < 1e-12);
    }
}
