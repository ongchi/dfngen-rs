use rand::distributions::{Distribution, Uniform};
use std::sync::OnceLock;

static MAX_VALUE: OnceLock<f64> = OnceLock::new();

fn max_val_below_one() -> f64 {
    let rep = (f64::MANTISSA_DIGITS as f64 * std::f64::consts::LOG10_2).floor() as usize + 2;
    format!("0.{}", "9".repeat(rep)).parse().unwrap()
}

/// The exponential distribution
pub struct Exp {
    lambda_inverse: f64,
    uniform: Uniform<f64>,
}

#[derive(Debug)]
pub enum Error {
    LambdaTooSmall,
    InvalidMinValue,
    InvalidMaxValue,
}

impl Exp {
    pub fn new(lambda: f64, min: f64, max: f64) -> Result<Exp, Error> {
        if lambda < 0. {
            return Err(Error::LambdaTooSmall);
        }
        if !(0. ..=1.).contains(&min) {
            return Err(Error::InvalidMinValue);
        }
        if !(0. ..=1.).contains(&max) {
            return Err(Error::InvalidMaxValue);
        }

        Ok(Exp {
            lambda_inverse: 1. / lambda,
            uniform: Uniform::new(min, max),
        })
    }
}

impl Distribution<f64> for Exp {
    fn sample<R: rand::prelude::Rng + ?Sized>(&self, rng: &mut R) -> f64 {
        let rand_val = rng.sample(self.uniform);
        let rand_val = match rand_val {
            1. => *MAX_VALUE.get_or_init(max_val_below_one),
            val => val,
        };

        -(1. - rand_val).ln() * self.lambda_inverse
    }
}
