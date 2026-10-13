#[derive(thiserror::Error, Debug)]
pub enum DfngenError {
    #[error(transparent)]
    SystemTime(#[from] std::time::SystemTimeError),
    #[error(transparent)]
    Io(#[from] std::io::Error),
    #[error(transparent)]
    DistrUniform(#[from] rand::distr::uniform::Error),
    #[error(transparent)]
    LogNormalSampling(#[from] crate::distribution::log_norm::Error),
    #[error(transparent)]
    ExponentialSampling(#[from] crate::distribution::exp::Error),
    #[error("Family {family}: cannot sample a fracture radius: {source}")]
    FamilyRadius {
        family: usize,
        source: Box<DfngenError>,
    },
    #[error("Too many fractures are being generated with radii less than minimum radius.")]
    InefficientSampling,
    #[error("Invalid input file: {0}")]
    InvalidInput(String),
    #[error("No fracture is defined. Add a fracture family or a user-defined fracture in the domain to the input file.")]
    NoFractureDefined,
    #[error("There are no intersecting fractures. Increase the fracture density, or decrease the domain size.")]
    NoIntersectingFractures,
}
