use rand::distributions::{Distribution, Uniform};

/// The truncated power-law distribution
pub struct TruncPowerLaw {
    min: f64,
    max: f64,
    alpha: f64,
    uniform: Uniform<f64>,
}

impl TruncPowerLaw {
    pub fn new(min: f64, max: f64, alpha: f64) -> Self {
        Self {
            min,
            max,
            alpha,
            uniform: Uniform::new(0., 1.),
        }
    }
}

impl Distribution<f64> for TruncPowerLaw {
    fn sample<R: rand::Rng + ?Sized>(&self, rng: &mut R) -> f64 {
        let rand_val = rng.sample(self.uniform);
        let tmp = 1. - rand_val + (rand_val * (self.min / self.max).powf(self.alpha));
        self.min * tmp.powf(-1. / self.alpha)
    }
}
