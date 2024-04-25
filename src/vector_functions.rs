use core::panic;

use crate::read_input::Input;

// Calculates crossproduct of v1 and v2
// Arg 1: Pointer to array of three elements
// Arg 2: Pointer to array of three elements
// Return: Pointer to cross product, array of three elements
pub fn cross_product(v1: &[f64; 3], v2: &[f64; 3]) -> [f64; 3] {
    [
        v1[1] * v2[2] - v1[2] * v2[1],
        v1[2] * v2[0] - v1[0] * v2[2],
        v1[0] * v2[1] - v1[1] * v2[0],
    ]
}

// Normalizes vector passed into fucntion.
// Arg 1: Vector (3 element array) to be normalized.
pub fn normalize(vec: &mut [f64; 3]) {
    let inv_mag = 1.0 / (vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]).sqrt();

    if inv_mag.is_infinite() {
        panic!("ERROR: Attempted to normalize a vector with magnitude = 0")
    }

    vec[0] *= inv_mag;
    vec[1] *= inv_mag;
    vec[2] *= inv_mag;
}

// Calculates the dot product of vector A with B
// Arg 1: Pointer to array of three elements
// Arg 2: Pointer to array of three elements
// Return: Pointer to dot product, array of three elements
pub fn dot_product(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

// // Prints vertices with {x, y, z} format.
// // Arg 1: Pointer to vertice array. Expects array
// //        length to be a multple of three
// // Arg 2: Number of vertices in array */
// pub fn printVertices(vert: &[&[f64; 3]]) {
//     println!("Vertices:");
//     for v in vert {
//         println!("\t{{{}, {}, {}}}", v[0], v[1], v[2])
//     }
// }

// Calculates magnitude of a vector
// Arg 1: x
// Arg 2: y
// Arg 3: z
pub fn magnitude(x: f64, y: f64, z: f64) -> f64 {
    ((x * x) + (y * y) + (z * z)).sqrt()
}

// Calculates the square magnitude of a vector
// Arg 1: x
// Arg 2: y
// Arg 3: z
// Return: Square magnitude of {x,y,z}
pub fn sqr_magnitude(x: f64, y: f64, z: f64) -> f64 {
    (x * x) + (y * y) + (z * z)
}

// Check if two vectors are parallel
// Arg 1: Pointer to vector 1, array of three doubles
// Arg 2: Pointer to vector 2, array of three doubles
// Output: True if vectors are parallel
//         False otherwise
pub fn parallel(input: &Input, v1: &mut [f64; 3], v2: &mut [f64; 3]) -> bool {
    normalize(v1);
    normalize(v2);
    let dot_prod = dot_product(v1, v2);

    1. - input.eps < dot_prod && dot_prod < 1. + input.eps
}

// // Projection funcion without sqrt()
// // Projects v1 onto v2
// // Arg 1: Pointer to vector 1 containing {x,y,z}, all doubles
// // Arg 2: Pointer to vector 2 containing {x,y,z}, all doubles
// // Return: Resulting vector.
// pub fn projection(v1: &[f64; 3], v2: &[f64; 3]) -> [f64; 3] {
//     let v2_ls = v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2];
//     let mut result = [0., 0., 0.];
//
//     if v2_ls >= eps {
//         let temp = (v2[0] * v1[0] + v2[1] * v1[1] + v2[2] * v1[2]) / v2_ls;
//         result[0] = v2[0] * temp;
//         result[1] = v2[1] * temp;
//         result[2] = v2[2] * temp;
//     }
//
//     result
// }

// Calculates the distance between two points in 3D space.
// Arg 1: Pointer to 3 element double array {x,y,z}
// Arg 2: Pointer to 3 element double array {x,y,z}
// Return: Distance between point A and point B
pub fn euclidean_distance(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    let temp1 = a[0] - b[0];
    let temp2 = a[1] - b[1];
    let temp3 = a[2] - b[2];
    (temp1 * temp1 + temp2 * temp2 + temp3 * temp3).sqrt()
}

// // Claculates the angle between two vectors
// // Arg 1: Pointer to 3 element double array, {x,y,z} vector
// // Arg 2: Pointer to 3 element double array, {x,y,z} vector
// // Return: Angle between both vectors in radians
// pub fn angleBeteenVectors(vector1: &[f64; 3], vector2: &[f64; 3]) -> f64 {
//     // cos(x) = u.v / (mag(u) * mag(v))
//     // x = arccos(x) = (mag(u) * mag(v)) / u.v
//     let dot = dotProduct(vector1, vector2);
//     let v1_mag = magnitude(vector1[0], vector1[1], vector1[2]);
//     let v2_mag = magnitude(vector2[0], vector2[1], vector2[2]);
//
//     (dot / (v1_mag * v2_mag)).acos()
// }
