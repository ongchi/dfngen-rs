use std::{cell::RefCell, rc::Rc};

use parry3d_f64::bounding_volume::Aabb;
use parry3d_f64::na::Vector3;
use rand::distr::Uniform;
use rand::Rng;
use rand_mt::Mt64;

/// Returns random position within a bounding box
pub fn random_position(boundary: Aabb, rng: Rc<RefCell<Mt64>>) -> Vector3<f64> {
    let distribution_x = Uniform::new(boundary.mins.x, boundary.maxs.x).unwrap();
    let distribution_y = Uniform::new(boundary.mins.y, boundary.maxs.y).unwrap();
    let distribution_z = Uniform::new(boundary.mins.z, boundary.maxs.z).unwrap();

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
pub fn generate_theta(aspect_ratio: f64, n_points: usize) -> Vec<f64> {
    let a = 1.;
    let b = aspect_ratio;
    let mut temp1 = (a - b) / (a + b);
    temp1 = temp1 * temp1;
    let c = std::f64::consts::PI * (a + b) * (1. + (3. * temp1) / (10. + (4. - 3. * temp1).sqrt()));
    let del = c / (n_points as f64);

    let mut theta_array: Vec<f64> = Vec::new();
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

    theta_array
}
