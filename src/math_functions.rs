use crate::read_input::Input;
use crate::structures::Poly;
use crate::vector_functions::{cross_product, magnitude};

// Used for ORing arrays of bool for boundary face codes
// Expectes array size is 6.
// ORs dest with src, then saves to dest.
// e.g destArray = destArray ^ srcArray
//
// Arg 1: dest array
// Arg 2: souce array
pub fn or(dest: &mut [bool; 6], src: &[bool; 6]) {
    for (d, s) in dest.iter_mut().zip(src.iter()) {
        *d |= *s
    }
}

// There are several spots in the code which used to use standard
// deviation to know how to sort a list of points.
// To increase performance, the standard deviation is now calculated without
// the square root and division by N elements.
// Return: The sum-deviation of array 'vec'
fn sum_deviation(vec: &[f64]) -> f64 {
    let mean = vec.iter().cloned().sum::<f64>() / (vec.len() as f64);
    vec.iter().cloned().map(|v| v - mean).map(|v| v * v).sum()
}

// *************  Sum deviaiton array X 3  **************************
// Used in intersectionChecking()
// How function is used:
// v is a pointer to an array of 4 points, each point
// contains x, y, and z coordinates
// The 4 sumDeviation is computed on all points' x, y, and z coord.
// An array of three elements is returned containing the sum deviation
// of all x's, all y's, and all z's
//
// Arg 1: Array of 12 elements: 4 points, {x1, y1, z1, ... , x4, y4, z4}
// Return: Array of 3 elements x,y,z with each stdDev, respectively */
pub fn sum_dev_ary3(v: &[f64; 12]) -> [f64; 3] {
    let x = [v[0], v[3], v[6], v[9]];
    let y = [v[1], v[4], v[7], v[10]];
    let z = [v[2], v[5], v[8], v[11]];
    [sum_deviation(&x), sum_deviation(&y), sum_deviation(&z)]
}

// *******************  Get Max Element's Index  ********************
// Used to find index of element with max value from array4
// Used in intersectionChecking()
// Arg 1: Array of doubles
// Arg 2: Size of array
// Return: Max element's array index
//
// If elements happen to be equal, returns the first one
pub fn max_elmt_idx(vec: &[f64]) -> usize {
    vec.iter()
        .enumerate()
        .reduce(|max, el| if el.1 < max.1 { max } else { el })
        .map(|(i, _)| i)
        .unwrap_or(0)
}

// ********************  Sort Array Indices  ************************
// Similar to mathematica's Ordering[] funct.
// Arg 1: Pointer to array of doubles
// Arg 2: Size of array, number of elements
// Return: Array of sorted indicies to elements in 'v'
//         sorted smallest to largest  */
pub fn sorted_index(vec: &[f64]) -> Vec<usize> {
    let mut v = vec.iter().enumerate().collect::<Vec<(usize, &f64)>>();
    v.sort_by(|a, b| (a.1).partial_cmp(b.1).unwrap());
    v.into_iter().map(|(i, _)| i).collect()
}

// ********************  Get Poly's Area  ***************************
// Calculate exact area of polygon (after truncation)
// Summary of algorithm:
// 1: Creates a point on the inside of the polygon (insidePt)
// 2: Uses 'insidePt' to create vectors to all outside vertices, breaking the polygon into triangles
// 3: Uses .5 * (magnitude of cross product) for each triangle for area calculation
// 4: Sums the areas for each trianlge for total area of polygon
// Arg 1: Polygon
// Return: Area of polygon */
pub fn get_area(poly: &Poly) -> f64 {
    if poly.number_of_nodes == 3 {
        //area = 1/2 mag of xProd
        let v1 = [
            poly.vertices[3] - poly.vertices[0],
            poly.vertices[4] - poly.vertices[1],
            poly.vertices[5] - poly.vertices[2],
        ];
        let v2 = [
            poly.vertices[6] - poly.vertices[0],
            poly.vertices[7] - poly.vertices[1],
            poly.vertices[8] - poly.vertices[2],
        ];
        let x_prod = cross_product(&v1, &v2);

        0.5 * magnitude(x_prod[0], x_prod[1], x_prod[2])
    } else {
        // More than 3 vertices
        let mut poly_area = 0.; // For summing area over trianlges of polygon
                                // Get coordinate within polygon
        let idx_across = (poly.number_of_nodes / 2 * 3) as usize;
        let inside_pt = [
            poly.vertices[0] + (0.5 * (poly.vertices[idx_across] - poly.vertices[0])),
            poly.vertices[1] + (0.5 * (poly.vertices[idx_across + 1] - poly.vertices[1])),
            poly.vertices[2] + (0.5 * (poly.vertices[idx_across + 2] - poly.vertices[2])),
        ];

        for i in 0..poly.number_of_nodes {
            let idx = (i * 3) as usize;
            let v1 = [
                poly.vertices[idx] - inside_pt[0],
                poly.vertices[idx + 1] - inside_pt[1],
                poly.vertices[idx + 2] - inside_pt[2],
            ];
            let v2 = [
                poly.vertices[idx + 3] - inside_pt[0],
                poly.vertices[idx + 4] - inside_pt[1],
                poly.vertices[idx + 5] - inside_pt[2],
            ];
            let x_prod = cross_product(&v1, &v2);
            let area = 0.5 * magnitude(x_prod[0], x_prod[1], x_prod[2]);
            poly_area += area; // Accumulate area
        }

        // Last portion of polygon, insidePt to first vertice and insidePt to last vertice
        let last = (3 * (poly.number_of_nodes - 1)) as usize;
        let v1 = [
            poly.vertices[0] - inside_pt[0],
            poly.vertices[1] - inside_pt[1],
            poly.vertices[2] - inside_pt[2],
        ];
        let v2 = [
            poly.vertices[last] - inside_pt[0],
            poly.vertices[last + 1] - inside_pt[1],
            poly.vertices[last + 2] - inside_pt[2],
        ];
        let x_prod = cross_product(&v1, &v2);
        let area = 0.5 * magnitude(x_prod[0], x_prod[1], x_prod[2]);
        poly_area += area; // Accumulate area

        poly_area
    }
}

// *************  Index from Probability  ***************************
// CDF is 1 to 1 and algined with the stochastic family
// shapes array (std vector). This chooses the family index based
// on a random roll (random number between 0 and 1) and 'famProb'
// Arg 1: CDF of shape families based on famProb array in input file
// Arg 2: Random number between 0 and 1
// Return: Index to family based on the probability famProb  */
pub fn index_from_prob(cdf: &[f64], roll: f64) -> usize {
    for (i, v) in cdf.iter().enumerate() {
        if &roll <= v {
            return i;
        }
    }

    cdf.len() - 1
}

// ******  Choose Family Randomly Based On P32 and CDF  ***********************
// Use with fracture intensity (p32) stopping option
// Arg 1: Pointer to CDF array
// Arg 2: Random variable between 0 and 1
// Arg 3: Number of stochastic families
// Arg 4: Number of elements in CDF array
// Arg 5: OUTPUT, index of CDF element which was chosen randomly
// Return: Index of chosen family based on its family probability, and p32Status
//         (avoids inserting a fracture from a family which has already met it's fracture intensity requrement) */
pub fn index_from_prob_and_p32_status(
    input: &mut Input,
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

    for i in 0..fam_size {
        input.p32Status[i] = false;
    }

    for i in 0..fam_size {
        if cdf_idx == &count && !input.p32Status[i] {
            return i; // Returns family index we need to build poly with.
        } else if !input.p32Status[i] {
            count += 1; // Count number of 0's ( number of families not having met their p32 req. )
        }
    }

    println!("ERROR: see indexFromProb_and_P32Status(), funct did not return anything");

    fam_size - 1
}

// ******  Choose Family Randomly Based On P32 and CDF  ***********************
// Given a family number, return its corresponding CDF array index
// Arg 1: CDF array
// Arg 2: p32Status array
// Arg 3: Family index of family whos CDF index to return
// Return: Index of CDF which belongs to shapeFamily[famIdx] */
pub fn cdf_idx_from_fam_num(p32_status: &[bool], fam_idx: usize) -> usize {
    let mut idx = 0;
    // Check how many 0's (p32's not complete) before the given family index
    // This gives the index in the CDF array for the given family
    for i in (0..fam_idx).rev() {
        if !p32_status[i] {
            idx += 1;
        }
    }

    idx - 1
}

// ************************  Create CDF *****************************
// Creates CDF from famProb[]
// Arg 1: Pointer to famProb array (see input file, and readInput())
pub fn create_cdf(fam_prob: &[f64]) -> Vec<f64> {
    // Convert famProb to CDF
    let mut cdf = vec![0.; fam_prob.len()];
    cdf[0] = fam_prob[0];

    for i in 1..fam_prob.len() {
        cdf[i] = cdf[i - 1] + fam_prob[i];
    }

    if (cdf.last().unwrap() < &0.999) || (cdf.last().unwrap() > &1.001) {
        println!("WARNING: Familiy probabilities (famProb in input file) do not sum to 1");
        println!("sum = {:.17}", cdf.last().unwrap());
        println!("Please check input file.");
    }

    cdf
}

// FIXME: vector size changed, return new vector instead of change inplace.
// Adjusts the CDF and the famProb array. Used with P32 stopping contdition (see input file)
// When a family's P32 requrement is met, the CDF is adjusted to remove that family from
// the array. The famProb is also adjusted in the same way.
// Arg 1: Pointer to CDF array
// Arg 2: Pointer to famProb array
// Arg 3: Number of elements in CDF array
// Arg 4: Index to the element in the famProb array which is being removed */
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
    for (t, s) in cdf.iter_mut().zip(create_cdf(fam_probability)) {
        *t = s;
    }
}
