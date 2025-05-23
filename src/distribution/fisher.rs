use std::fmt::{Display, Formatter};

use parry3d_f64::na::Vector3;
use rand::distr::{Distribution, Uniform};

#[derive(Debug)]
enum AngleOption {
    ThetaPhi(f64, f64),
    TrendPlunge(f64, f64),
    DipStrike(f64, f64),
}

/// The fisher distribution of polygon normal vector generation
#[derive(Debug)]
pub struct Fisher {
    angle_option: Option<AngleOption>,
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
            angle_option: Some(angle_option),
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
            Some(AngleOption::ThetaPhi(theta, phi)) => {
                writeln!(f, "Theta-rad: {}", theta)?;
                writeln!(f, "Theta-deg: {}", theta * rad_to_deg)?;
                writeln!(f, "Phi-rad: {}", phi)?;
                writeln!(f, "Phi-deg: {}", phi * rad_to_deg)?;
            }
            Some(AngleOption::TrendPlunge(trend, plunge)) => {
                writeln!(f, "Trend-rad: {}", trend)?;
                writeln!(f, "Trend-deg: {}", trend * rad_to_deg)?;
                writeln!(f, "Plunge-rad: {}", plunge)?;
                writeln!(f, "Plunge-deg: {}", plunge * rad_to_deg)?;
            }
            Some(AngleOption::DipStrike(dip, strike)) => {
                writeln!(f, "Dip-rad: {}", dip)?;
                writeln!(f, "Dip-deg: {}", dip * rad_to_deg)?;
                writeln!(f, "Strike-rad: {}", strike)?;
                writeln!(f, "Strike-deg: {}", strike * rad_to_deg)?;
            }
            None => unreachable!(),
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
        let r = if !(x_prod[0].abs() <= self.eps
            && x_prod[1].abs() <= self.eps
            && x_prod[2].abs() <= self.eps)
        {
            // Since vectors are normalized, sin = magnitude(AxB) and cos = A . B
            let sin =
                (x_prod[0] * x_prod[0] + x_prod[1] * x_prod[1] + x_prod[2] * x_prod[2]).sqrt();
            let cos = u.dot(&self.orientation);
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

        let theta_random = rng.sample(self.theta_uniform);
        let y = rng.sample(self.uniform);

        let w = 1. / self.kappa * ((-self.kappa).exp() + self.kappa * ck * y).ln();
        let temp = (1. - (w * w)).sqrt();
        let v = [temp * theta_random.cos(), temp * theta_random.sin()];

        // Matrix multiply with R
        Vector3::new(
            v[0] * r[0] + v[1] * r[1] + w * r[2],
            v[0] * r[3] + v[1] * r[4] + w * r[5],
            v[0] * r[6] + v[1] * r[7] + w * r[8],
        )
    }
}
