use parry3d_f64::na::{Rotation3, Vector3};
use rand::distr::{Distribution, Uniform};

use crate::cg::{normal_from_dip_strike, normal_from_trend_plunge, rotation_between};

/// The fisher distribution of polygon normal vector generation
#[derive(Debug)]
pub struct Fisher {
    orientation: Vector3<f64>,
    kappa: f64,
    eps: f64,
    uniform: Uniform<f64>,
    theta_uniform: Uniform<f64>,
}

impl Fisher {
    fn new_with_orientation(orientation: Vector3<f64>, kappa: f64, eps: f64) -> Self {
        Self {
            orientation,
            kappa,
            eps,
            uniform: Uniform::new(0., 1.).unwrap(),
            theta_uniform: Uniform::new(0., 2. * std::f64::consts::PI).unwrap(),
        }
    }

    pub fn new_with_theta_phi(theta: f64, phi: f64, kappa: f64, eps: f64) -> Self {
        Self::new_with_orientation(
            Vector3::new(
                theta.sin() * phi.cos(),
                theta.sin() * phi.sin(),
                theta.cos(),
            ),
            kappa,
            eps,
        )
    }

    pub fn new_with_trend_plunge(trend: f64, plunge: f64, kappa: f64, eps: f64) -> Self {
        Self::new_with_orientation(normal_from_trend_plunge(trend, plunge), kappa, eps)
    }

    pub fn new_with_dip_strike(dip: f64, strike: f64, kappa: f64, eps: f64) -> Self {
        Self::new_with_orientation(normal_from_dip_strike(dip, strike), kappa, eps)
    }
}

impl Distribution<Vector3<f64>> for Fisher {
    fn sample<R: rand::Rng + ?Sized>(&self, rng: &mut R) -> Vector3<f64> {
        let ck = (self.kappa.exp() - (-self.kappa).exp()) / self.kappa;

        let u = Vector3::new(0., 0., 1.);
        let x_prod = u.cross(&self.orientation);

        // Rotation is necessary only if the normal vectors are not the same (xProd is not zero vector)
        let rotation = if x_prod.magnitude() > self.eps {
            rotation_between(&u, &self.orientation)
        } else {
            Rotation3::identity()
        };

        let theta_random = rng.sample(self.theta_uniform);
        let y = rng.sample(self.uniform);

        let w = 1. / self.kappa * ((-self.kappa).exp() + self.kappa * ck * y).ln();
        let temp = (1. - (w * w)).sqrt();
        let v = [temp * theta_random.cos(), temp * theta_random.sin()];

        rotation * Vector3::new(v[0], v[1], w)
    }
}
