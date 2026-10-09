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
    #[error("Invalid input file: {0}")]
    InvalidInput(String),
    #[error("There are no intersecting fractures. Increase the fracture density, or decrease the domain size.")]
    NoIntersectingFractures,
    #[error("The fracture network does not connect the specified boundary faces. Increase the fracture density, decrease the domain size, or use the clusters.ignore_boundary_faces option.")]
    NoBoundaryConnection,
}
