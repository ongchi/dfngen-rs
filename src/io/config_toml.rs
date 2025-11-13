/// TOML Configuration Format Support
///
/// This module provides serde-compatible structures for parsing TOML configuration files.
/// The TOML format provides a more structured, readable, and maintainable alternative
/// to the original custom text-based format.

use serde::{Deserialize, Serialize};
use parry3d_f64::na::Vector3;

/// Top-level TOML configuration structure
#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct TomlConfig {
    pub domain: DomainConfig,
    pub stopping_condition: StoppingConditionConfig,
    pub fram: FramConfig,
    pub output: OutputConfig,
    pub fractures: FracturesConfig,
    pub boundaries: BoundariesConfig,

    #[serde(default)]
    pub ellipse_families: Vec<EllipseFamilyConfig>,

    #[serde(default)]
    pub rectangle_families: Vec<RectangleFamilyConfig>,

    #[serde(default)]
    pub user_fractures: Option<UserFracturesConfig>,

    #[serde(default)]
    pub optional: Option<OptionalConfig>,
}

/// Domain configuration
#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct DomainConfig {
    /// Domain size (x, y, z) in meters
    pub size: DomainSize,

    /// Minimum feature size for FRAM algorithm (h parameter)
    #[serde(default = "default_h")]
    pub h: f64,

    /// Size increase for fractures extending beyond domain
    #[serde(default)]
    pub size_increase: DomainSize,
}

fn default_h() -> f64 {
    1.0
}

/// Domain size with x, y, z components
#[derive(Debug, Serialize, Deserialize, Clone)]
#[serde(untagged)]
pub enum DomainSize {
    /// Format: {x = 100, y = 100, z = 100}
    Named { x: f64, y: f64, z: f64 },
    /// Format: [100, 100, 100]
    Array([f64; 3]),
}

impl DomainSize {
    pub fn to_vector3(&self) -> Vector3<f64> {
        match self {
            DomainSize::Named { x, y, z } => Vector3::new(*x, *y, *z),
            DomainSize::Array([x, y, z]) => Vector3::new(*x, *y, *z),
        }
    }
}

impl Default for DomainSize {
    fn default() -> Self {
        DomainSize::Named {
            x: 100.0,
            y: 100.0,
            z: 100.0,
        }
    }
}

/// Stopping condition configuration
#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct StoppingConditionConfig {
    /// Stopping condition: "nPoly" (by count) or "P32" (by fracture intensity)
    pub mode: String,

    /// Number of polygons to generate (used when mode = "nPoly")
    #[serde(default)]
    pub n_poly: usize,

    /// P32 target (used when mode = "P32")
    #[serde(default)]
    pub p32_target: Option<f64>,

    /// Percentage increase to pre-generated radii lists (e.g., 0.2 = 20%)
    #[serde(default = "default_radii_increase")]
    pub radii_list_increase: f64,
}

fn default_radii_increase() -> f64 {
    0.2
}

/// FRAM algorithm configuration
#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct FramConfig {
    /// Enable FRAM (Feature Rejection Algorithm for Meshing)
    #[serde(default = "default_true")]
    pub enable: bool,

    /// Use relaxed FRAM (less strict validation)
    #[serde(default)]
    pub relaxed: bool,

    /// Check and enforce triple intersection constraints
    #[serde(default = "default_true")]
    pub triple_intersections: bool,

    /// Enable visualization mode (first triangulation round only)
    #[serde(default)]
    pub visualization_mode: bool,
}

fn default_true() -> bool {
    true
}

/// Output configuration
#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct OutputConfig {
    /// Print fracture rejection reasons to console
    #[serde(default)]
    pub print_reject_reasons: bool,

    /// Output all generated radii (accepted + rejected)
    #[serde(default)]
    pub output_all_radii: bool,

    /// Output accepted radii per family (before isolated fracture removal)
    #[serde(default)]
    pub output_accepted_radii_per_family: bool,

    /// Output final radii per family (after isolated fracture removal)
    #[serde(default)]
    pub output_final_radii_per_family: bool,
}

/// Fractures configuration
#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct FracturesConfig {
    /// Number of re-translation attempts per fracture
    #[serde(default = "default_rejects_per_fracture")]
    pub rejects_per_fracture: usize,

    /// Force insertion of largest fractures per family
    #[serde(default)]
    pub force_large_fractures: bool,

    /// Remove fractures smaller than this radius
    #[serde(default)]
    pub remove_smaller_than: f64,

    /// Orientation definition option: "spherical" (theta/phi), "trend_plunge", or "dip_strike"
    #[serde(default = "default_orientation_option")]
    pub orientation_option: String,
}

fn default_rejects_per_fracture() -> usize {
    1
}

fn default_orientation_option() -> String {
    "spherical".to_string()
}

/// Boundary faces configuration
#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct BoundariesConfig {
    /// Which boundary faces must be connected (order: +x, -x, +y, -y, +z, -z)
    #[serde(default = "default_boundary_faces")]
    pub faces: BoundaryFaces,

    /// Keep only largest cluster connecting boundary faces
    #[serde(default)]
    pub keep_only_largest_cluster: bool,

    /// Keep isolated fractures not connected to any boundary
    #[serde(default = "default_keep_isolated")]
    pub keep_isolated_fractures: bool,

    /// Ignore boundary face connection requirements
    #[serde(default)]
    pub ignore_boundary_faces: bool,
}

fn default_boundary_faces() -> BoundaryFaces {
    BoundaryFaces {
        faces: [true, true, true, true, true, true],
    }
}

fn default_keep_isolated() -> bool {
    true
}

/// Boundary faces array wrapper
#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct BoundaryFaces {
    /// Array of 6 booleans for boundary faces [+x, -x, +y, -y, +z, -z]
    pub faces: [bool; 6],
}

impl BoundaryFaces {
    pub fn to_array(&self) -> [bool; 6] {
        self.faces
    }
}

/// Ellipse fracture family configuration
#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct EllipseFamilyConfig {
    /// Family name (optional, for reference)
    #[serde(default)]
    pub name: Option<String>,

    /// Probability of selecting this family during generation
    #[serde(default)]
    pub probability: f64,

    /// Number of vertices to approximate the ellipse
    #[serde(default = "default_num_points")]
    pub num_points: u8,

    /// Target P32 (fracture intensity) for this family
    #[serde(default)]
    pub p32_target: Option<f64>,

    /// Aspect ratio (yradius / xradius)
    #[serde(default)]
    pub aspect_ratio: f64,

    /// Beta rotation angle (degrees)
    #[serde(default)]
    pub beta: f64,

    /// Beta has a distribution (true) or is constant (false)
    #[serde(default)]
    pub beta_distribution: bool,

    /// Layer ID (0 = whole domain, >0 = specific layer)
    #[serde(default)]
    pub layer: usize,

    /// Region ID (0 = whole domain, >0 = specific region)
    #[serde(default)]
    pub region: usize,

    /// Orientation parameters
    pub orientation: OrientationConfig,

    /// Radius distribution parameters
    pub radius: RadiusConfig,
}

fn default_num_points() -> u8 {
    8
}

/// Rectangle fracture family configuration
#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct RectangleFamilyConfig {
    /// Family name (optional, for reference)
    #[serde(default)]
    pub name: Option<String>,

    /// Probability of selecting this family during generation
    #[serde(default)]
    pub probability: f64,

    /// Target P32 (fracture intensity) for this family
    #[serde(default)]
    pub p32_target: Option<f64>,

    /// Aspect ratio (yradius / xradius)
    #[serde(default)]
    pub aspect_ratio: f64,

    /// Beta rotation angle (degrees)
    #[serde(default)]
    pub beta: f64,

    /// Beta has a distribution (true) or is constant (false)
    #[serde(default)]
    pub beta_distribution: bool,

    /// Layer ID (0 = whole domain, >0 = specific layer)
    #[serde(default)]
    pub layer: usize,

    /// Region ID (0 = whole domain, >0 = specific region)
    #[serde(default)]
    pub region: usize,

    /// Orientation parameters
    pub orientation: OrientationConfig,

    /// Radius distribution parameters
    pub radius: RadiusConfig,
}

/// Orientation configuration (angle distribution)
#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct OrientationConfig {
    /// Primary angle: theta (spherical), trend (trend/plunge), or dip (dip/strike)
    pub angle1: f64,

    /// Secondary angle: phi (spherical), plunge (trend/plunge), or strike (dip/strike)
    pub angle2: f64,

    /// Concentration parameter for Fisher distribution (kappa)
    #[serde(default)]
    pub kappa: f64,
}

/// Radius distribution configuration
#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct RadiusConfig {
    /// Distribution type: "lognormal", "power_law", "exponential", or "constant"
    pub distribution_type: String,

    // Lognormal parameters
    #[serde(default)]
    pub log_mean: Option<f64>,
    #[serde(default)]
    pub log_std: Option<f64>,
    #[serde(default)]
    pub log_min: Option<f64>,
    #[serde(default)]
    pub log_max: Option<f64>,

    // Power-law parameters
    #[serde(default)]
    pub alpha: Option<f64>,

    // Exponential parameters
    #[serde(default)]
    pub exp_mean: Option<f64>,
    #[serde(default)]
    pub exp_min: Option<f64>,
    #[serde(default)]
    pub exp_max: Option<f64>,

    // Constant value
    #[serde(default)]
    pub constant: Option<f64>,

    // General min/max bounds
    #[serde(default)]
    pub min: Option<f64>,
    #[serde(default)]
    pub max: Option<f64>,
}

/// User-defined fractures configuration
#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct UserFracturesConfig {
    /// Insert user-defined rectangles before ellipses
    #[serde(default)]
    pub insert_rectangles_first: bool,

    /// Path to file with user-defined ellipses
    #[serde(default)]
    pub ellipses_file: Option<String>,

    /// Path to file with user-defined rectangles
    #[serde(default)]
    pub rectangles_file: Option<String>,

    /// Path to file with user-defined polygons (by coordinates)
    #[serde(default)]
    pub polygons_by_coord_file: Option<String>,

    /// Path to file with user-defined ellipses (by coordinates)
    #[serde(default)]
    pub ellipses_by_coord_file: Option<String>,

    /// Path to file with user-defined rectangles (by coordinates)
    #[serde(default)]
    pub rectangles_by_coord_file: Option<String>,
}

/// Optional configuration sections
#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct OptionalConfig {
    /// Seed for random number generator (0 = use system time)
    #[serde(default)]
    pub seed: u64,

    /// Layer boundaries configuration
    #[serde(default)]
    pub layers: Option<LayersConfig>,

    /// Region boundaries configuration
    #[serde(default)]
    pub regions: Option<RegionsConfig>,

    /// Polygon boundary configuration
    #[serde(default)]
    pub polygon_boundary: Option<PolygonBoundaryConfig>,
}

/// Layer boundaries configuration
#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct LayersConfig {
    /// Layer boundary z-values
    pub boundaries: Vec<f64>,
}

/// Region boundaries configuration
#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct RegionsConfig {
    /// Region boundary coordinates [xmin, xmax, ymin, ymax, zmin, zmax, ...]
    pub boundaries: Vec<f64>,
}

/// Polygon boundary configuration (2D domain boundary)
#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct PolygonBoundaryConfig {
    /// Path to file with polygon boundary vertices
    pub vertices_file: String,
}
