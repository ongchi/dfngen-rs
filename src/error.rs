use crate::distribution::aperture::ApertureDistribution;

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
    #[error("Invalid input file: cannot read {file}: {source}")]
    ReadInput {
        file: String,
        source: std::io::Error,
    },
    /// The TOML text: syntax, an unknown or a missing key, a wrong type
    #[error("Invalid input file: {0}")]
    InvalidToml(String),
    /// The content of the domain mesh (`.vtu`) or of the user-defined fracture file (`.vtp`)
    #[error("Invalid input file: {file}: {message}")]
    InvalidFile { file: String, message: String },
    #[error("Invalid input file: the text of the user-defined fracture file {file} is necessary")]
    UserFractureFileMissing { file: String },
    #[error("Invalid input file: min_feature_size must be greater than 0")]
    MinFeatureSize,
    /// `fracture` starts at 1, in the order of the user-defined fracture file
    #[error("Invalid input file: user-defined fracture {fracture}: {source}")]
    UserFracture {
        fracture: usize,
        source: UserFractureError,
    },
    /// The polygons of one fracture (same `fracture_id`) are one fracture of the mesh
    #[error(
        "Invalid input file: the pieces of the user-defined fracture with fracture_id {fracture_id} have different apertures ({first} and {other}). Give them one aperture."
    )]
    UserFracturePieceApertures {
        fracture_id: u64,
        first: f64,
        other: f64,
    },
    /// `family` starts at 1, in the order of the input file
    #[error("Invalid input file: family {family}: {source}")]
    Family { family: usize, source: FamilyError },
    #[error(
        "Invalid input file: [mesh] needs [region_boundaries] with enable = true: its boundary lines give the boundaries of the mesh"
    )]
    MeshNeedsRegionBoundaries,
    #[error(
        "Invalid input file: clusters.keep_only_connected_to needs [region_boundaries] with enable = true: its boundary lines give the boundary surfaces which a fracture touches"
    )]
    ClustersNeedRegionBoundaries,
    #[error(
        "Invalid input file: clusters.keep_only_connected_to needs a domain mesh with boundary surfaces (cell data arrays boundary__<name>)"
    )]
    ClustersNeedBoundarySurfaces,
    /// `group` starts at 1, in the order of the input file
    #[error("Invalid input file: group {group} of clusters.keep_only_connected_to has no boundary surface. Add one.")]
    EmptyClusterSurfaceGroup { group: usize },
    #[error(
        "Invalid input file: clusters.keep_only_connected_to: the domain mesh has no boundary surface {name} (cell data array boundary__{name})"
    )]
    UnknownClusterSurface { name: String },
    #[error("Invalid input file: mesh.max_size_factor must be 1 or more")]
    MeshMaxSizeFactor,
    #[error("Invalid input file: mesh.slope must be greater than 0")]
    MeshSlope,
    #[error("No fracture is defined. Add a fracture family or a user-defined fracture in the domain to the input file.")]
    NoFractureDefined,
    #[error("Mesh: {0}")]
    Mesh(String),
    #[error("There are no intersecting fractures. Increase the fracture density, or decrease the domain size.")]
    NoIntersectingFractures,
    /// The bytes of a checkpoint
    #[error("Invalid checkpoint: {0}")]
    CheckpointFormat(String),
    #[error("The checkpoint is of another input: the input files changed after the checkpoint. Use the input files of the run which wrote it.")]
    CheckpointInputChanged,
    /// `fracture` starts at 1, in the order of the accepted stochastic fractures
    #[error("The restart from the checkpoint failed: FRAM rejected accepted fracture {fracture} ({reason}). The checkpoint and this version of dfngen do not agree.")]
    CheckpointRebuildFailed {
        fracture: usize,
        reason: &'static str,
    },
}

/// A region name of a family or of a user-defined fracture which the domain mesh cannot give
#[derive(thiserror::Error, Debug)]
pub enum RegionError {
    #[error("the domain mesh has no region {region} (no cell data array region__{region})")]
    Unknown { region: String },
    #[error("no cell of the domain mesh is in the region {region}")]
    Empty { region: String },
}

/// A value of a user-defined fracture which is not valid
#[derive(thiserror::Error, Debug)]
pub enum UserFractureError {
    #[error(transparent)]
    Region(#[from] RegionError),
    #[error("the aperture {aperture} of the aperture array must be above 0")]
    Aperture { aperture: f64 },
}

/// A value of a fracture family which is not valid
#[derive(thiserror::Error, Debug)]
pub enum FamilyError {
    #[error(transparent)]
    Region(#[from] RegionError),
    #[error("number_of_points must be 3 or more")]
    NumberOfPoints,
    #[error("aspect_ratio must be greater than 0")]
    AspectRatio,
    #[error("orientation must have one pair: theta and phi, trend and plunge, or dip and strike")]
    OrientationPair,
    #[error("p32_target is necessary when generation.n_poly is not set")]
    P32TargetMissing,
    #[error("p32_target is not allowed when generation.n_poly is set")]
    P32TargetNotAllowed,
    #[error(transparent)]
    Aperture(#[from] InvalidAperture),
    /// A value which the input code did not give to the family builder (a defect of the code)
    #[error("{0} is required")]
    MissingValue(&'static str),
}

/// An aperture distribution with a value out of its range
#[derive(thiserror::Error, Debug)]
#[error("aperture {aperture:?}: value and alpha must be above 0, sd must not be negative")]
pub struct InvalidAperture {
    pub aperture: ApertureDistribution,
}
