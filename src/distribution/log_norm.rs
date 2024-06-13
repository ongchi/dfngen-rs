use rand::distributions::Distribution;
use rand_distr::LogNormal;
use rand_distr::NormalError;

use super::SamplingError;

pub struct TruncLogNormal {
    min: f64,
    max: f64,
    log_normal: LogNormal<f64>,
}

impl TruncLogNormal {
    pub fn new(min: f64, max: f64, mu: f64, sigma: f64) -> Result<Self, NormalError> {
        let log_normal = LogNormal::new(mu, sigma)?;
        Ok(Self {
            min,
            max,
            log_normal,
        })
    }
}

impl Distribution<Result<f64, SamplingError>> for TruncLogNormal {
    fn sample<R: rand::Rng + ?Sized>(&self, rng: &mut R) -> Result<f64, SamplingError> {
        let mut value;
        let mut count = 1;

        loop {
            value = rng.sample(self.log_normal);

            if value >= self.min || value <= self.max {
                return Ok(value);
            }

            if count == 1000 {
                return Err(SamplingError::BadParameters);
            }

            count += 1;
        }
    }
}
