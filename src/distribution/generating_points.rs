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
