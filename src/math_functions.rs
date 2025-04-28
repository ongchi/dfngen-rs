use tracing::{error, warn};

/// Used for ORing arrays of bool for boundary face codes
/// Expectes array size is 6.
/// ORs dest with src, then saves to dest.
/// e.g destArray = destArray ^ srcArray
///
/// # Arguments
///
/// * `dest` - Destination array
/// * `src` - Source array
pub fn or(dest: &mut [bool; 6], src: &[bool; 6]) {
    for (d, s) in dest.iter_mut().zip(src.iter()) {
        *d |= *s
    }
}

/// There are several spots in the code which used to use standard
/// deviation to know how to sort a list of points.
/// To increase performance, the standard deviation is now calculated without
/// the square root and division by N elements.
/// Return: The sum-deviation of array 'vec'
fn sum_deviation(vec: &[f64]) -> f64 {
    let mean = vec.iter().cloned().sum::<f64>() / (vec.len() as f64);
    vec.iter().cloned().map(|v| v - mean).map(|v| v * v).sum()
}

/// Sum deviaiton array X 3
///
/// Used in intersectionChecking()
///
/// How function is used:
/// v is a pointer to an array of 4 points, each point
/// contains x, y, and z coordinates
/// The 4 sumDeviation is computed on all points' x, y, and z coord.
/// An array of three elements is returned containing the sum deviation
/// of all x's, all y's, and all z's
///
/// # Arguments
///
/// * `v` - Array of 12 elements: 4 points, {x1, y1, z1, ... , x4, y4, z4}
///
/// # Returns
///
/// Array of 3 elements x,y,z with each stdDev, respectively
pub fn sum_dev_ary3(v: &[f64; 12]) -> [f64; 3] {
    let x = [v[0], v[3], v[6], v[9]];
    let y = [v[1], v[4], v[7], v[10]];
    let z = [v[2], v[5], v[8], v[11]];
    [sum_deviation(&x), sum_deviation(&y), sum_deviation(&z)]
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
/// * `p32_status` - Array of bools, 1 to 1 with number of families
/// * `cdf` - CDF of shape families based on famProb array in input file
/// * `roll` - Random number between 0 and 1
/// * `fam_size` - Number of stochastic families
/// * `cdf_idx` - Index of CDF element which was chosen randomly
///
/// # Returns
///
/// Index of chosen family based on its family probability, and p32Status
/// (avoids inserting a fracture from a family which has already met it's fracture intensity requrement)
pub fn index_from_prob_and_p32_status(
    p32_status: &mut [bool],
    cdf: &[f64],
    roll: f64,
    fam_size: usize,
    cdf_idx: &mut usize,
) -> usize {
    // The p32Status bool array stays 1 to 1 with the number of total families.
    // The p32Status array with element = 1 means that family has reached its p32 requirement.
    // The CDF contains an amount of elements equal to the number of families which has not met its intensity requirement.
    // To get the correct familyShape[] index from cdf, me must ignore any families their p32 requirement alrady met
    // Index we hit based on random roll,
    *cdf_idx = index_from_prob(cdf, roll);
    // If cdfIndex = 2 (thrid element including 0)
    // we need to find the families with p32Status still 0, and choose the thrid one
    // this family will be used to create the next polygon
    let mut count = 0; // Count of families encountered with p32Status = 0

    for p32 in p32_status.iter_mut() {
        *p32 = false;
    }

    for (i, p32) in p32_status.iter().enumerate() {
        if cdf_idx == &count {
            return i; // Returns family index we need to build poly with.
        } else if !p32 {
            count += 1; // Count number of 0's ( number of families not having met their p32 req. )
        }
    }

    error!("see indexFromProb_and_P32Status(), funct did not return anything");

    fam_size - 1
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
            warn!("Familiy probabilities (famProb in input file) do not sum to 1");
            warn!("sum = {:.17}", cdf.last().unwrap());
            warn!("Please check input file.");
        }
    };

    cdf
}

/// FIXME: vector size changed, return new vector instead of change inplace.
///
/// Adjusts the CDF and the famProb array. Used with P32 stopping contdition (see input file)
/// When a family's P32 requrement is met, the CDF is adjusted to remove that family from
/// the array. The famProb is also adjusted in the same way.
///
/// # Arguments
///
/// * `cdf` - CDF array
/// * `fam_probability` - Probability of each family
/// * `cdf_size` - Number of elements in CDF array
/// * `idx2remove` - Index to the element in the famProb array which is being removed
pub fn adjust_cdf_and_fam_prob(
    cdf: &mut [f64],
    fam_probability: &mut [f64],
    cdf_size: &mut usize,
    idx2remove: usize,
) {
    *cdf_size -= 1;
    let mut new_probs = Vec::with_capacity(cdf.len());
    // Adjust probabilities, remove element while keeping the rest of probabilities in proportion
    // Take probability of the familiy being removed, and divide it equally among remaining probabilities
    let add_to_remaining_elmts = fam_probability[idx2remove] / (*cdf_size as f64); // Distribute removed probability among leftore famillies probabilities
                                                                                   // cdfIdx is index to that families cdf index AND familily probability index
    let mut idx = 0;

    for (i, prob) in fam_probability.iter().enumerate().take(*cdf_size + 1) {
        // cdfSize+1 becuase of cdfSize-- at the beginning of function. This
        // makes the loop cover the all probabilities in the ary before one is removed
        // and distributed to the rest
        if i != idx2remove {
            new_probs[idx] = prob + add_to_remaining_elmts;
            idx += 1;
        }
    }

    for (t, s) in fam_probability.iter_mut().zip(new_probs) {
        *t = s;
    }
    for (t, s) in cdf.iter_mut().zip(cumsum(fam_probability)) {
        *t = s;
    }
}
