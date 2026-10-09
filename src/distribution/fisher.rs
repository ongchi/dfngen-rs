use std::fmt::{Display, Formatter};

use parry3d_f64::na::Vector3;
use rand::distr::{Distribution, Uniform};

use crate::cg::{rotate, rotation_matrix};

#[derive(Debug)]
enum AngleOption {
    ThetaPhi(f64, f64),
    TrendPlunge(f64, f64),
    DipStrike(f64, f64),
}

/// The fisher distribution of polygon normal vector generation
#[derive(Debug)]
pub struct Fisher {
    angle_option: AngleOption,
    orientation: Vector3<f64>,
    kappa: f64,
    eps: f64,
    uniform: Uniform<f64>,
    theta_uniform: Uniform<f64>,
}

impl Fisher {
    fn new_with_orientation(
        orientation: Vector3<f64>,
        kappa: f64,
        eps: f64,
        angle_option: AngleOption,
    ) -> Self {
        Self {
            orientation,
            kappa,
            eps,
            uniform: Uniform::new(0., 1.).unwrap(),
            theta_uniform: Uniform::new(0., 2. * std::f64::consts::PI).unwrap(),
            angle_option,
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
            AngleOption::ThetaPhi(theta, phi),
        )
    }

    pub fn new_with_trend_plunge(trend: f64, plunge: f64, kappa: f64, eps: f64) -> Self {
        Self::new_with_orientation(
            Vector3::new(
                trend.cos() * plunge.cos(),
                trend.sin() * plunge.cos(),
                plunge.cos(),
            ),
            kappa,
            eps,
            AngleOption::TrendPlunge(trend, plunge),
        )
    }

    pub fn new_with_dip_strike(dip: f64, strike: f64, kappa: f64, eps: f64) -> Self {
        Self::new_with_orientation(
            Vector3::new(
                dip.sin() * strike.sin(),
                -dip.sin() * strike.cos(),
                dip.cos(),
            ),
            kappa,
            eps,
            AngleOption::DipStrike(dip, strike),
        )
    }
}

impl Display for Fisher {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let rad_to_deg = 180. / std::f64::consts::PI;

        match self.angle_option {
            AngleOption::ThetaPhi(theta, phi) => {
                writeln!(f, "Theta-rad: {}", theta)?;
                writeln!(f, "Theta-deg: {}", theta * rad_to_deg)?;
                writeln!(f, "Phi-rad: {}", phi)?;
                writeln!(f, "Phi-deg: {}", phi * rad_to_deg)?;
            }
            AngleOption::TrendPlunge(trend, plunge) => {
                writeln!(f, "Trend-rad: {}", trend)?;
                writeln!(f, "Trend-deg: {}", trend * rad_to_deg)?;
                writeln!(f, "Plunge-rad: {}", plunge)?;
                writeln!(f, "Plunge-deg: {}", plunge * rad_to_deg)?;
            }
            AngleOption::DipStrike(dip, strike) => {
                writeln!(f, "Dip-rad: {}", dip)?;
                writeln!(f, "Dip-deg: {}", dip * rad_to_deg)?;
                writeln!(f, "Strike-rad: {}", strike)?;
                writeln!(f, "Strike-deg: {}", strike * rad_to_deg)?;
            }
        }

        writeln!(f, "Kappa: {}", self.kappa)
    }
}

impl Distribution<Vector3<f64>> for Fisher {
    fn sample<R: rand::Rng + ?Sized>(&self, rng: &mut R) -> Vector3<f64> {
        let ck = (self.kappa.exp() - (-self.kappa).exp()) / self.kappa;

        let u = Vector3::new(0., 0., 1.);
        let x_prod = u.cross(&self.orientation);

        // Get rotation matrix if normal vectors are not the same (if xProd is not zero vector)
        let r = if x_prod.magnitude() > self.eps {
            rotation_matrix(&u, &self.orientation)
        } else {
            // Identity Matrix
            [1., 0., 0., 0., 1., 0., 0., 0., 1.]
        };

        let theta_random = rng.sample(self.theta_uniform);
        let y = rng.sample(self.uniform);

        let w = 1. / self.kappa * ((-self.kappa).exp() + self.kappa * ck * y).ln();
        let temp = (1. - (w * w)).sqrt();
        let v = [temp * theta_random.cos(), temp * theta_random.sin()];

        // Matrix multiply with R
        Vector3::from(rotate(&r, v[0], v[1], w))
    }
}
