use rand::distr::{Distribution, Uniform};
use std::sync::OnceLock;

use super::SamplingError;

static MAX_VALUE: OnceLock<f64> = OnceLock::new();

fn max_val_below_one() -> f64 {
    let rep = (f64::MANTISSA_DIGITS as f64 * std::f64::consts::LOG10_2).floor() as usize + 2;
    format!("0.{}", "9".repeat(rep)).parse().unwrap()
}

/// The exponential distribution
pub struct TruncExp {
    min: f64,
    max: f64,
    lambda_inverse: f64,
    uniform: Uniform<f64>,
}

#[derive(Debug)]
pub enum Error {
    LambdaTooSmall,
    InvalidNormMinValue,
    InvalidNormMaxValue,
}

impl TruncExp {
    pub fn new(
        min: f64,
        max: f64,
        lambda: f64,
        norm_min: f64,
        norm_max: f64,
    ) -> Result<TruncExp, Error> {
        if lambda < 0. {
            return Err(Error::LambdaTooSmall);
        }
        if !(0. ..=1.).contains(&norm_min) {
            return Err(Error::InvalidNormMinValue);
        }
        if !(0. ..=1.).contains(&norm_max) {
            return Err(Error::InvalidNormMaxValue);
        }

        Ok(TruncExp {
            min,
            max,
            lambda_inverse: 1. / lambda,
            uniform: Uniform::new(norm_min, norm_max).unwrap(),
        })
    }
}

impl Distribution<Result<f64, SamplingError>> for TruncExp {
    fn sample<R: rand::prelude::Rng + ?Sized>(&self, rng: &mut R) -> Result<f64, SamplingError> {
        let mut value;
        let mut count = 1;

        loop {
            let rand_val = rng.sample(self.uniform);
            let rand_val = match rand_val {
                1. => *MAX_VALUE.get_or_init(max_val_below_one),
                v => v,
            };
            value = -(1. - rand_val).ln() * self.lambda_inverse;

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
