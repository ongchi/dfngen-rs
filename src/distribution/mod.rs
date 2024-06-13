mod exp;
mod fisher;
pub mod generating_points;
mod log_norm;
mod power_law;

pub use exp::TruncExp;
pub use fisher::Fisher;
pub use log_norm::TruncLogNormal;
pub use power_law::TruncPowerLaw;

#[derive(Debug)]
pub enum SamplingError {
    BadParameters,
}
