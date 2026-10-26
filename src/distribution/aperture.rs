use rand::Rng;
use rand_distr::StandardNormal;
use rand_mt::Mt64;
use serde::Deserialize;

use crate::error::InvalidAperture;

/// Distribution of the hydraulic aperture of a fracture (m), as the
/// transmissivity options of dfnWorks, with the aperture in place of the
/// transmissivity. r is the radius of the fracture.
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(tag = "distribution", rename_all = "snake_case", deny_unknown_fields)]
pub enum ApertureDistribution {
    /// The same aperture for each fracture
    Constant { value: f64 },
    /// ln(aperture) has a normal distribution
    LogNormal { log_mean: f64, sd: f64 },
    /// aperture = alpha × r^beta
    Correlated { alpha: f64, beta: f64 },
    /// ln(aperture) = ln(alpha × r^beta) + sd × N(0, 1)
    SemiCorrelated { alpha: f64, beta: f64, sd: f64 },
}

impl ApertureDistribution {
    /// Checks the parameters. The aperture must be finite and above 0.
    pub fn validate(&self) -> Result<(), InvalidAperture> {
        let is_positive = |value: f64| value.is_finite() && value > 0.;
        let is_valid_sd = |sd: f64| sd.is_finite() && sd >= 0.;

        let is_valid = match *self {
            Self::Constant { value } => is_positive(value),
            Self::LogNormal { log_mean, sd } => log_mean.is_finite() && is_valid_sd(sd),
            Self::Correlated { alpha, beta } => is_positive(alpha) && beta.is_finite(),
            Self::SemiCorrelated { alpha, beta, sd } => {
                is_positive(alpha) && beta.is_finite() && is_valid_sd(sd)
            }
        };

        if is_valid {
            Ok(())
        } else {
            Err(InvalidAperture {
                aperture: self.clone(),
            })
        }
    }

    /// Aperture of a fracture with the radius. Only `log_normal` and
    /// `semi_correlated` use the generator.
    pub fn sample(&self, radius: f64, generator: &mut Mt64) -> f64 {
        match *self {
            Self::Constant { value } => value,
            Self::LogNormal { log_mean, sd } => {
                let normal: f64 = generator.sample(StandardNormal);
                (log_mean + sd * normal).exp()
            }
            Self::Correlated { alpha, beta } => alpha * radius.powf(beta),
            Self::SemiCorrelated { alpha, beta, sd } => {
                let normal: f64 = generator.sample(StandardNormal);
                alpha * radius.powf(beta) * (sd * normal).exp()
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_correlated_aperture_is_alpha_times_radius_to_beta() {
        let distribution = ApertureDistribution::Correlated {
            alpha: 1e-4,
            beta: 0.5,
        };

        let aperture = distribution.sample(4., &mut Mt64::new(1));

        assert!((aperture - 2e-4).abs() < 1e-18);
    }

    #[test]
    fn test_log_normal_aperture_without_spread_is_exp_of_log_mean() {
        let distribution = ApertureDistribution::LogNormal {
            log_mean: -7.,
            sd: 0.,
        };

        let aperture = distribution.sample(1., &mut Mt64::new(1));

        assert!((aperture - (-7f64).exp()).abs() < 1e-18);
    }

    #[test]
    fn test_zero_constant_aperture_is_not_valid() {
        let distribution = ApertureDistribution::Constant { value: 0. };

        assert!(distribution.validate().is_err());
    }
}
