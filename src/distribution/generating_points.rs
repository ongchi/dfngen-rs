use std::{cell::RefCell, rc::Rc};

use parry3d_f64::na::{Point3, Vector3};
use rand::distr::Uniform;
use rand::Rng;
use rand_mt::Mt64;

/// Discretize Intersection
///
/// # Arguments
///
/// * `h` - Minimum feature size
/// * `visualization_mode` - If false, creates a fine mesh, according to h parameter.
///     If true, produce only first round of triangulations. In this case no modeling of flow and transport is possible.
/// * `pt1` - End point 1
/// * `pt2` - End point 2
///
/// # Returns
///
/// List of 3D points of the discretized nodes, including end points
pub fn discretize_line_of_intersection(
    h: f64,
    visualization_mode: bool,
    pt1: &Point3<f64>,
    pt2: &Point3<f64>,
    dist: f64,
) -> Vec<Point3<f64>> {
    // If reduced mesh, just save endpoints
    if visualization_mode {
        let mut points_list = Vec::new();
        let pt = Point3::new(pt1[0], pt1[1], pt1[2]);
        points_list.push(pt);
        let pt = Point3::new(pt2[0], pt2[1], pt2[2]);
        points_list.push(pt);
        return points_list;
    }

    let v = pt2 - pt1;
    let p = *pt1;

    let nprime = (2. * dist / h).ceil();
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
        points_list.push(p + v * x);
    }

    points_list
}

/// Returns random TRANSLATION
///
/// # Arguments
///
/// * `rng` - Random generator, see std c++ <random> library
/// * `x_min` - Minimum x for random x
/// * `x_max` - Maximum x for random x
/// * `y_min` - Minimum y for random y
/// * `y_max` - Maximum y for random y
/// * `z_min` - Minimum z for random z
/// * `z_max` - Maximum z for random z
pub fn random_translation(
    rng: Rc<RefCell<Mt64>>,
    x_min: f64,
    x_max: f64,
    y_min: f64,
    y_max: f64,
    z_min: f64,
    z_max: f64,
) -> Vector3<f64> {
    let distribution_x = Uniform::new(x_min, x_max).unwrap();
    let distribution_y = Uniform::new(y_min, y_max).unwrap();
    let distribution_z = Uniform::new(z_min, z_max).unwrap();

    let x = rng.borrow_mut().sample(distribution_x);
    let y = rng.borrow_mut().sample(distribution_y);
    let z = rng.borrow_mut().sample(distribution_z);

    Vector3::new(x, y, z)
}

/// Generates Theta Array for Generating Ellipses
///
/// Integrate diff eq for theta as function of arc length using RK2
/// Used once for each ell family, saves theta array to shape structures
///
/// # Arguments
///
/// * `theta_array` - Theta array used for ellipse generation
/// * `aspect_ratio` - Aspect ratio of ellipse family
/// * `n_points` - Number of points being used for ellipse family
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
