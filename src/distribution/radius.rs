use rand::Rng;
use rand_mt::Mt64;

use crate::distribution::{LogNormalSampling, TruncExp, TruncLogNormal, TruncPowerLaw};
use crate::error::DfngenError;

pub enum RadiusFunction {
    LogNormal {
        mu: f64,
        sigma: f64,
        sampling: LogNormalSampling,
    },
    TruncatedPowerLaw {
        alpha: f64,
    },
    Exponential {
        lambda: f64,
    },
    Constant(f64),
}

pub struct RadiusDistribution {
    pub min: f64,
    pub max: f64,
    pub function: RadiusFunction,
}

impl RadiusDistribution {
    pub fn new_truncated_power_law(alpha: f64, min: f64, max: f64) -> Self {
        Self {
            function: RadiusFunction::TruncatedPowerLaw { alpha },
            min,
            max,
        }
    }

    pub fn new_log_normal(
        mu: f64,
        sigma: f64,
        min: f64,
        max: f64,
        sampling: LogNormalSampling,
    ) -> Self {
        Self {
            function: RadiusFunction::LogNormal {
                mu,
                sigma,
                sampling,
            },
            min,
            max,
        }
    }

    pub fn new_exponential(lambda: f64, min: f64, max: f64) -> Self {
        Self {
            function: RadiusFunction::Exponential { lambda },
            min,
            max,
        }
    }

    pub fn new_constant(value: f64) -> Self {
        Self {
            function: RadiusFunction::Constant(value),
            min: value,
            max: value,
        }
    }

    /// Sampling from distribution function.
    pub fn sample(&self, generator: &mut Mt64) -> Result<f64, DfngenError> {
        let radius = match self.function {
            RadiusFunction::LogNormal {
                mu,
                sigma,
                sampling,
            } => {
                let distr = TruncLogNormal::new(self.min, self.max, mu, sigma, sampling)?;
                generator.sample(&distr)?
            }
            RadiusFunction::TruncatedPowerLaw { alpha } => {
                let distr = TruncPowerLaw::new(self.min, self.max, alpha);
                generator.sample(&distr)
            }
            RadiusFunction::Exponential { lambda } => {
                let distr = TruncExp::new(self.min, self.max, lambda)?;
                generator.sample(&distr)?
            }
            RadiusFunction::Constant(c) => c,
        };

        Ok(radius)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_samples_are_within_limits(distribution: &RadiusDistribution) {
        let mut generator = Mt64::new(1);

        for _ in 0..10_000 {
            let radius = distribution.sample(&mut generator).unwrap();
            assert!(
                (distribution.min..=distribution.max).contains(&radius),
                "radius {radius} is not in [{}, {}]",
                distribution.min,
                distribution.max
            );
        }
    }

    #[test]
    fn test_log_normal_radius_is_within_min_and_max() {
        assert_samples_are_within_limits(&RadiusDistribution::new_log_normal(
            1.,
            1.,
            1.,
            4.,
            LogNormalSampling::Rejection,
        ));
    }

    #[test]
    fn test_exponential_radius_is_within_min_and_max() {
        assert_samples_are_within_limits(&RadiusDistribution::new_exponential(0.4, 1., 4.));
    }

    #[test]
    fn test_truncated_power_law_radius_is_within_min_and_max() {
        assert_samples_are_within_limits(&RadiusDistribution::new_truncated_power_law(2.5, 1., 4.));
    }
}
