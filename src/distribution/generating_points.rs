use std::{cell::RefCell, rc::Rc};

use crate::io::input::Input;
use parry3d::na::{Point3, Vector3};
use rand::distributions::Uniform;
use rand::Rng;
use rand_mt::Mt19937GenRand64;

// ********************  Discretize Intersection  ***************************
// Discretizes intersetion
// Arg 1: End point 1, array of three doubles {x, y, z}
// Arg 2: End point 2, array of three doubles {x, y, z}
// Return: List of 3D points of the discretized nodes, including end points */
pub fn discretize_line_of_intersection(
    input: &Input,
    pt1: &Point3<f64>,
    pt2: &Point3<f64>,
    dist: f64,
) -> Vec<Point3<f64>> {
    // If reduced mesh, just save endpoints
    if input.visualizationMode {
        let mut points_list = Vec::new();
        let pt = Point3::new(pt1[0], pt1[1], pt1[2]);
        points_list.push(pt);
        let pt = Point3::new(pt2[0], pt2[1], pt2[2]);
        points_list.push(pt);
        return points_list;
    }

    let v = pt2 - pt1;
    let p = *pt1;

    let nprime = (2. * dist / input.h).ceil();
    let hprime = 1. / nprime;
    let nprime = nprime as usize;
    let mut xx = Vec::with_capacity(nprime + 1);
    let mut temp = 0.;

    for _ in 0..nprime {
        // Array from 0 - 1 with step = hprime
        xx.push(temp);
        temp += hprime;
    }

    if let Some(last) = xx.last_mut() {
        *last = 1.;
    }
    let mut points_list = Vec::with_capacity(nprime + 1);

    for x in xx.iter().cloned().take(points_list.capacity()) {
        points_list.push(line_function_3d(&v, &p, x));
    }

    points_list
}

// *********************** Parametric Line Function *************************
// Returns a point on the line/vector v at point t,  0 <= t <= 1
// Arg 1: Array of three doubles {x, y, z}, vector of line segment
// Arg 2: Array of three doubles {x, y, z}, end point on line
// Arg 3: t, 0<= t <= 1
// Return: Point on line
pub fn line_function_3d(v: &Vector3<f64>, point: &Point3<f64>, t: f64) -> Point3<f64> {
    Point3::new(
        point[0] + v[0] * t,
        point[1] + v[1] * t,
        point[2] + v[2] * t,
    )
}

// ****** Fisher Distributions for Generating polygons Normal Vectors *******
// Creates and returns an x,y,z array of doubles using Fisher distribution.
//
// NOTE: Uses new[] to return an array. NEED TO USE delete[] TO FREE THE MEMORY AFTER USE
//
// Arg 1: theta, the angle the normal vector makes with the z-axis
// Arg 2: phi, the angle the projection of the normal onto the x-y plane makes with the x-axis
// Arg 3: kappa, parameter for the Fisher distribnShaprutions
// Arg 4: Random generator, see std c++ <random> library
// Return: A Fisher distribution array {x, y, z}. Used for random generation of
// polygon normal vectors.
pub fn fisher_distribution(
    input: &Input,
    angle1: f64,
    angle2: f64,
    kappa: f64,
    rng: Rc<RefCell<Mt19937GenRand64>>,
) -> Vector3<f64> {
    let ck = (kappa.exp() - (-kappa).exp()) / kappa;

    let v1 = if input.orientationOption == 0 {
        // Spherical Coordinates
        // angleOne = Theta
        // angleTwo = Phi
        Vector3::new(
            angle1.sin() * angle2.cos(),
            angle1.sin() * angle2.sin(),
            angle1.cos(),
        )
    } else if input.orientationOption == 1 {
        // Trend and Plunge
        // angleOne = Trend
        // angleTwo = Plunge
        Vector3::new(
            angle1.cos() * angle2.cos(),
            angle1.sin() * angle2.cos(),
            angle2.cos(),
        )
    } else if input.orientationOption == 2 {
        // Dip and Strike
        // angleOne = Dip
        // angleTwo = Strike
        Vector3::new(
            angle1.sin() * angle2.sin(),
            -angle1.sin() * angle2.cos(),
            angle1.cos(),
        )
    } else {
        unreachable!()
    };

    let u = Vector3::new(0., 0., 1.);
    let x_prod = u.cross(&v1);

    // Get rotation matrix if normal vectors are not the same (if xProd is not zero vector)
    let r = if !(x_prod[0].abs() <= input.eps
        && x_prod[1].abs() <= input.eps
        && x_prod[2].abs() <= input.eps)
    {
        // Since vectors are normalized, sin = magnitude(AxB) and cos = A . B
        let sin = (x_prod[0] * x_prod[0] + x_prod[1] * x_prod[1] + x_prod[2] * x_prod[2]).sqrt();
        let cos = u.dot(&v1);
        let v = [
            0., -x_prod[2], x_prod[1], x_prod[2], 0., -x_prod[0], -x_prod[1], x_prod[0], 0.,
        ];
        let scalar = (1.0 - cos) / (sin * sin);
        let v_squared = [
            (v[0] * v[0] + v[1] * v[3] + v[2] * v[6]) * scalar,
            (v[0] * v[1] + v[1] * v[4] + v[2] * v[7]) * scalar,
            (v[0] * v[2] + v[1] * v[5] + v[2] * v[8]) * scalar,
            (v[3] * v[0] + v[4] * v[3] + v[5] * v[6]) * scalar,
            (v[3] * v[1] + v[4] * v[4] + v[5] * v[7]) * scalar,
            (v[3] * v[2] + v[4] * v[5] + v[5] * v[8]) * scalar,
            (v[6] * v[0] + v[7] * v[3] + v[8] * v[6]) * scalar,
            (v[6] * v[1] + v[7] * v[4] + v[8] * v[7]) * scalar,
            (v[6] * v[2] + v[7] * v[5] + v[8] * v[8]) * scalar,
        ];

        [
            1. + v[0] + v_squared[0],
            0. + v[1] + v_squared[1],
            0. + v[2] + v_squared[2],
            0. + v[3] + v_squared[3],
            1. + v[4] + v_squared[4],
            0. + v[5] + v_squared[5],
            0. + v[6] + v_squared[6],
            0. + v[7] + v_squared[7],
            1. + v[8] + v_squared[8],
        ]
    } else {
        // Identity Matrix
        [1., 0., 0., 0., 1., 0., 0., 0., 1.]
    };

    let theta_dist = Uniform::new(0., 2. * std::f64::consts::PI);
    let distribution = Uniform::new(0., 1.);
    let theta_random = rng.borrow_mut().sample(theta_dist);
    let y = rng.borrow_mut().sample(distribution);

    let w = 1. / kappa * ((-kappa).exp() + kappa * ck * y).ln();
    let temp = (1. - (w * w)).sqrt();
    let v = [temp * theta_random.cos(), temp * theta_random.sin()];

    // Matrix multiply with R
    Vector3::new(
        v[0] * r[0] + v[1] * r[1] + w * r[2],
        v[0] * r[3] + v[1] * r[4] + w * r[5],
        v[0] * r[6] + v[1] * r[7] + w * r[8],
    )
}

// ******************* Returns random TRANSLATION ***************************
// Uses new[] to pass a vector/array. NEED TO USE delete[] TO FREE THE MEMORY AFTER USE
// Arg 1: Random generator, see std c++ <random> library
// Arg 2: minimum x for random x
// Arg 3: maximum x for random x
// Arg 4: maximum y for random y
// Arg 5: minimum y for random y
// Arg 6: maximum z for random z
// Arg 7: minimum z for random z
// Return: Pointer to random ranslation, array of three doubles {x, y, z}
pub fn random_translation(
    rng: Rc<RefCell<Mt19937GenRand64>>,
    x_min: f64,
    x_max: f64,
    y_min: f64,
    y_max: f64,
    z_min: f64,
    z_max: f64,
) -> [f64; 3] {
    let distribution_x = Uniform::new(x_min, x_max);
    let distribution_y = Uniform::new(y_min, y_max);
    let distribution_z = Uniform::new(z_min, z_max);

    let x = rng.borrow_mut().sample(distribution_x);
    let y = rng.borrow_mut().sample(distribution_y);
    let z = rng.borrow_mut().sample(distribution_z);

    [x, y, z]
}

// *******************  Truncated Power-Law  ********************************
// Distrubution function for truncated power-law
// randomNum must be between 0 and 1
// This distribution should be sampled uniformly between 0 and 1 to
// produce a truncated power-law distribution
// Arg 1: Random variable between 0 and 1
// Arg 2: Minimum number which can be returned from the distribution
// Arg 3: Maximum number which can be returned from the distribution
// Arg 4: Power-law's alpha
// Return: Random float adhering to the truncated power law distribution */
pub fn truncated_power_law(random_num: f64, min: f64, max: f64, alpha: f64) -> f64 {
    let temp = 1. - random_num + (random_num * (min / max).powf(alpha));
    min * temp.powf(-1. / alpha)
}

// **********  Generates Theta Array for Generating Ellipses  ***************
// Integrate diff eq for theta as function of arc length using RK2
// Used once for each ell family, saves theta array to shape structures
// Arg 1: OUTPUT, Theta array used for ellipse generation
// Arg 2: Aspect ratio of ellipse family
// Arg 3: Number of points being used for ellipse family */
pub fn generate_theta(theta_array: &mut Vec<f64>, aspect_ratio: f64, n_points: usize) {
    let a = 1.;
    let b = aspect_ratio;
    let mut temp1 = (a - b) / (a + b);
    temp1 = temp1 * temp1;
    let c = std::f64::consts::PI * (a + b) * (1. + (3. * temp1) / (10. + (4. - 3. * temp1).sqrt()));
    let del = c / (n_points as f64);

    theta_array.push(0.);
    for i in 1..n_points {
        let mut tmp =
            (b * theta_array[i - 1].cos()).powf(2.) + (a * theta_array[i - 1].sin()).powf(2.);
        let f_tmp = del / tmp.sqrt();
        tmp = theta_array[i - 1] + f_tmp;
        theta_array.push(
            theta_array[i - 1]
                + 0.5
                    * del
                    * ((b * tmp.cos()).powf(2.) + (a * tmp.sin()).powf(2.))
                        .sqrt()
                        .powf(-1.)
                + 0.5 * f_tmp,
        );
    }
}
