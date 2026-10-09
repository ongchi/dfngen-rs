use rand::distr::Distribution;
use rand_distr::LogNormal;
use serde::Deserialize;
use statrs::distribution::{ContinuousCDF, Normal};

/// The method which gives the samples of the truncated log-normal distribution
#[derive(Deserialize, Default, Debug, Clone, Copy, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum LogNormalSampling {
    /// Sample the full distribution and reject each value outside the range.
    /// Fast when the range contains most of the distribution.
    #[default]
    Rejection,
    /// One uniform number for each sample. The time does not depend on the range.
    InverseCdf,
}

#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error("The standard deviation of the log-normal distribution must be greater than 0")]
    SigmaTooSmall,
    #[error("The range [{min}, {max}] is too far from the mean of the log-normal distribution. Move the range nearer to the mean.")]
    RangeTooFarFromMean { min: f64, max: f64 },
    #[error("Too many samples of the log-normal distribution are outside the range [{min}, {max}]. Set `sampling = \"inverse_cdf\"` in the `radius` table of the family.")]
    InefficientSampling { min: f64, max: f64 },
}

/// The truncated log-normal distribution
pub struct TruncLogNormal {
    min: f64,
    max: f64,
    method: Method,
}

enum Method {
    Rejection(LogNormal<f64>),
    InverseCdf {
        mu: f64,
        sigma: f64,
        /// For a range in the upper tail, the CDF values are near 1 and have few
        /// significant digits. Then the range is mirrored to the lower tail.
        is_mirrored: bool,
        cdf_low: f64,
        cdf_high: f64,
    },
}

impl TruncLogNormal {
    pub fn new(
        min: f64,
        max: f64,
        mu: f64,
        sigma: f64,
        sampling: LogNormalSampling,
    ) -> Result<Self, Error> {
        if sigma.is_nan() || sigma <= 0. {
            return Err(Error::SigmaTooSmall);
        }

        let method = match sampling {
            LogNormalSampling::Rejection => {
                Method::Rejection(LogNormal::new(mu, sigma).map_err(|_| Error::SigmaTooSmall)?)
            }
            LogNormalSampling::InverseCdf => inverse_cdf_method(min, max, mu, sigma)?,
        };

        Ok(Self { min, max, method })
    }

    fn sample_with_rejection<R: rand::Rng + ?Sized>(
        &self,
        log_normal: LogNormal<f64>,
        rng: &mut R,
    ) -> Result<f64, Error> {
        for _ in 0..1000 {
            let value = rng.sample(log_normal);

            if value >= self.min && value <= self.max {
                return Ok(value);
            }
        }

        Err(Error::InefficientSampling {
            min: self.min,
            max: self.max,
        })
    }
}

fn inverse_cdf_method(min: f64, max: f64, mu: f64, sigma: f64) -> Result<Method, Error> {
    let z_min = (min.ln() - mu) / sigma;
    let z_max = (max.ln() - mu) / sigma;
    let is_mirrored = z_min > 0.;
    let (z_low, z_high) = if is_mirrored {
        (-z_max, -z_min)
    } else {
        (z_min, z_max)
    };

    let cdf_low = Normal::standard().cdf(z_low);
    let cdf_high = Normal::standard().cdf(z_high);
    if min < max && cdf_low >= cdf_high {
        return Err(Error::RangeTooFarFromMean { min, max });
    }

    Ok(Method::InverseCdf {
        mu,
        sigma,
        is_mirrored,
        cdf_low,
        cdf_high,
    })
}

impl Distribution<Result<f64, Error>> for TruncLogNormal {
    fn sample<R: rand::Rng + ?Sized>(&self, rng: &mut R) -> Result<f64, Error> {
        match self.method {
            Method::Rejection(log_normal) => self.sample_with_rejection(log_normal, rng),
            Method::InverseCdf {
                mu,
                sigma,
                is_mirrored,
                cdf_low,
                cdf_high,
            } => {
                let cdf_value = cdf_low + (cdf_high - cdf_low) * rng.random::<f64>();
                let z_low_tail = Normal::standard().inverse_cdf(cdf_value);
                let z = if is_mirrored { -z_low_tail } else { z_low_tail };

                // Rounding can put the value a little outside the limits.
                Ok(f64::exp(mu + sigma * z).clamp(self.min, self.max))
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use rand::Rng;
    use rand_mt::Mt64;

    use super::LogNormalSampling::{InverseCdf, Rejection};
    use super::*;

    fn sample_values(distr: &TruncLogNormal, count: usize) -> Vec<f64> {
        let mut generator = Mt64::new(1);
        (0..count)
            .map(|_| generator.sample(distr).unwrap())
            .collect()
    }

    fn fraction_at_or_below(values: &[f64], limit: f64) -> f64 {
        values.iter().filter(|&&value| value <= limit).count() as f64 / values.len() as f64
    }

    fn assert_samples_are_inside_the_limits(sampling: LogNormalSampling) {
        let (min, max) = (1.4, 5.6);
        let distr = TruncLogNormal::new(min, max, 0.85, 0.4, sampling).unwrap();

        let values = sample_values(&distr, 10_000);

        assert!(values.iter().all(|&value| value >= min && value <= max));
    }

    /// Z = ln(X) is standard normal, truncated to [-1, 2].
    fn assert_central_samples_follow_the_truncated_cdf(sampling: LogNormalSampling) {
        let distr = TruncLogNormal::new(f64::exp(-1.), f64::exp(2.), 0., 1., sampling).unwrap();

        let values = sample_values(&distr, 100_000);

        // (Phi(0) - Phi(-1)) / (Phi(2) - Phi(-1))
        let expected = 0.416989;
        assert!((fraction_at_or_below(&values, 1.) - expected).abs() < 0.01);
    }

    #[test]
    fn test_rejection_samples_are_inside_the_limits() {
        assert_samples_are_inside_the_limits(Rejection);
    }

    #[test]
    fn test_inverse_cdf_samples_are_inside_the_limits() {
        assert_samples_are_inside_the_limits(InverseCdf);
    }

    #[test]
    fn test_rejection_samples_follow_the_truncated_cdf() {
        assert_central_samples_follow_the_truncated_cdf(Rejection);
    }

    #[test]
    fn test_inverse_cdf_samples_follow_the_truncated_cdf() {
        assert_central_samples_follow_the_truncated_cdf(InverseCdf);
    }

    #[test]
    fn test_inverse_cdf_upper_tail_samples_follow_the_truncated_cdf() {
        // Z = ln(X) is standard normal, truncated to [4, 5].
        let distr = TruncLogNormal::new(f64::exp(4.), f64::exp(5.), 0., 1., InverseCdf).unwrap();

        let values = sample_values(&distr, 100_000);

        // (Q(4) - Q(4.2)) / (Q(4) - Q(5)), with Q = 1 - Phi
        let expected = 0.58390;
        assert!((fraction_at_or_below(&values, f64::exp(4.2)) - expected).abs() < 0.01);
    }

    #[test]
    fn test_inverse_cdf_lower_tail_samples_follow_the_truncated_cdf() {
        // Z = ln(X) is standard normal, truncated to [-5, -4].
        let distr = TruncLogNormal::new(f64::exp(-5.), f64::exp(-4.), 0., 1., InverseCdf).unwrap();

        let values = sample_values(&distr, 100_000);

        // (Phi(-4.2) - Phi(-5)) / (Phi(-4) - Phi(-5))
        let expected = 0.41610;
        assert!((fraction_at_or_below(&values, f64::exp(-4.2)) - expected).abs() < 0.01);
    }

    #[test]
    fn test_inverse_cdf_with_equal_limits_gives_the_limit() {
        let distr = TruncLogNormal::new(2., 2., 0., 1., InverseCdf).unwrap();

        let values = sample_values(&distr, 10);

        assert!(values.iter().all(|&value| value == 2.));
    }

    #[test]
    fn test_rejection_in_a_tail_gives_an_error_which_names_inverse_cdf() {
        // Z = ln(X) is standard normal, truncated to [6, 7].
        let distr = TruncLogNormal::new(f64::exp(6.), f64::exp(7.), 0., 1., Rejection).unwrap();

        let error = Mt64::new(1).sample(&distr).unwrap_err();

        assert!(error.to_string().contains("sampling = \"inverse_cdf\""));
    }
}
