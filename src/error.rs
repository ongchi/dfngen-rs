use crate::structures::Shape;

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
    Sampling(#[from] crate::distribution::SamplingError),
    #[error(transparent)]
    ExponentialSampling(#[from] crate::distribution::exp::Error),
    #[error("{shape} Family {id} is attepting to populate fracture radii lists, however many fractures are being generated with radii less than minimum radius.")]
    TooManySmallFractures { shape: Shape, id: usize },
}
