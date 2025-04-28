#[derive(thiserror::Error, Debug)]
pub enum DfngenError {
    #[error(transparent)]
    SystemTime(#[from] std::time::SystemTimeError),
    #[error(transparent)]
    Io(#[from] std::io::Error),
    #[error(transparent)]
    DistrUniform(#[from] rand::distr::uniform::Error),
    #[error(transparent)]
    DistrNormal(#[from] rand_distr::NormalError),
    #[error(transparent)]
    ExponentialSampling(#[from] crate::distribution::exp::Error),
    #[error("Too many fractures are being generated with radii less than minimum radius.")]
    InefficientSampling,
    #[error("{0}")]
    ValueError(String),
}
