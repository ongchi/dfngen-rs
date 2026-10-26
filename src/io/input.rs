use std::collections::HashMap;
use std::path::Path;
use std::rc::Rc;

use parry3d_f64::na::Vector3;
use serde::Deserialize;
use tracing::{debug, info, warn};

use crate::distribution::{ApertureDistribution, Fisher, LogNormalSampling, RadiusDistribution};
use crate::domain_mesh::{DomainMesh, MeshRegion};
use crate::error::{DfngenError, FamilyError, RegionError, UserFractureError};
use crate::fracture::fracture_family::{
    FractureFamily, FractureFamilyBuilder, FractureFamilyCollection,
};
use crate::fracture::user_defined_fractures::UserFracture;
use crate::fram::FramOptions;
use crate::io::vtp::{parse_user_fractures, read_user_fractures, FracturePolygon};
use crate::io::vtu::{parse_domain_mesh, read_domain_mesh};
use crate::mesh::MeshOptions;
use crate::region_boundaries::RegionBoundaryOptions;

/// The input has `generation.n_poly` for `NPoly`, and no `n_poly` for `P32`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StopCondition {
    /// Stop at this number of accepted fractures.
    NPoly(usize),

    /// Stop when each family has its target fracture intensity.
    P32,
}

#[derive(Debug)]
pub struct Input {
    /// DFN generation stop condition.
    pub stop_condition: StopCondition,

    /// Size of the bounding box of the domain mesh. The generation has the
    /// center of this box at the origin.
    pub domain_size: Vector3<f64>,

    /// Volume of the domain: the sum of the cell volumes of the domain mesh.
    pub domain_volume: f64,

    /// The cells which make the domain.
    pub domain_mesh: Rc<DomainMesh>,

    /// Minimum feature size, FRAM parameter.
    pub min_feature_size: f64,

    /// epsilon
    pub eps: f64,

    /// Options of FRAM. `None`: FRAM is off, there will be no fracture rejections.
    pub fram: Option<FramOptions>,

    /// Options of the boundary lines of the regions. `None`: no lines.
    pub region_boundaries: Option<RegionBoundaryOptions>,

    /// Options of the triangle mesh of the network. `None`: no mesh.
    pub mesh: Option<MeshOptions>,

    /// 0 - Keep each cluster of 2 or more fractures
    /// 1 - Keep only the largest cluster
    pub keep_only_largest_cluster: bool,

    /// 0 - remove isolated fractures and clusters
    /// 1 - Keep isolated fractures and clusters
    pub keep_isolated_fractures: bool,

    /// Inserts the largest possible fracture for each defined fracture family,
    /// defined by the user-defined maxium radius
    ///     0 - Off (Do not force insertion of larest fractures)
    ///     1 - On  (Force insertion of largest fractures)
    pub force_large_fractures: bool,

    /// Seed for random generator.
    pub seed: u64,

    /// 0 - Ignore this option, keep all fractures.
    ///
    /// (>0) - Size of minimum fracture radius. Fractures smaller than
    /// defined radius will be removed AFTER DFN generation.
    ///
    /// Minimum and maximum size options under fracture family
    /// distributions will still be used while generating the DFN.
    pub remove_fractures_less_than: f64,

    /// Volume of each region of the domain mesh, in the order of
    /// `DomainMesh::regions`.
    pub region_volumes: Vec<f64>,

    /// User-defined fractures, in the order of their file. The vertices have
    /// the coordinates of the generation.
    pub user_fractures: Vec<UserFracture>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct InputFile {
    seed: u64,
    min_feature_size: f64,
    /// Path of an ASCII `.vtp` file with the user-defined fractures (optional).
    /// A relative path starts at the folder of the input file.
    user_fractures: Option<String>,
    /// Path of an ASCII `.vtu` file with the mesh which gives the domain.
    /// A relative path starts at the folder of the input file.
    domain: String,
    generation: GenerationSection,
    #[serde(default)]
    fram: FramSection,
    region_boundaries: RegionBoundariesSection,
    mesh: MeshSection,
    clusters: ClustersSection,
    #[serde(default, rename = "family")]
    families: Vec<FamilySection>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct GenerationSection {
    n_poly: Option<usize>,
    force_large_fractures: bool,
    remove_fractures_less_than: f64,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
#[serde(default)]
struct FramSection {
    enable: bool,
    relaxed: bool,
    triple_intersections: bool,
    check_user_fractures: bool,
    write_rejected_fractures: bool,
    rejects_per_fracture: usize,
}

/// Without the section or a key: FRAM is off, each option is false, and
/// `rejects_per_fracture` is 10, as in dfnWorks.
impl Default for FramSection {
    fn default() -> Self {
        Self {
            enable: false,
            relaxed: false,
            triple_intersections: false,
            check_user_fractures: false,
            write_rejected_fractures: false,
            rejects_per_fracture: 10,
        }
    }
}

impl FramSection {
    /// The FRAM options, or `None` if FRAM is off. The other keys then have no effect.
    fn options(&self) -> Option<FramOptions> {
        self.enable.then_some(FramOptions {
            relaxed: self.relaxed,
            triple_intersections: self.triple_intersections,
            check_user_fractures: self.check_user_fractures,
            write_rejected_fractures: self.write_rejected_fractures,
            rejects_per_fracture: self.rejects_per_fracture,
        })
    }
}

/// Necessary section and key: `enable = false` gives no boundary lines
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct RegionBoundariesSection {
    enable: bool,
}

impl RegionBoundariesSection {
    /// The options of the boundary lines, or `None` if the option is off
    fn options(&self) -> Option<RegionBoundaryOptions> {
        self.enable.then_some(RegionBoundaryOptions)
    }
}

/// Necessary section and keys: `enable = false` gives no mesh, and the other
/// keys then have no effect
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct MeshSection {
    enable: bool,
    max_size_factor: f64,
    slope: f64,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct ClustersSection {
    keep_only_largest: bool,
    keep_isolated_fractures: bool,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct FamilySection {
    number_of_points: usize,
    probability: f64,
    aspect_ratio: f64,
    /// Name of a region of the domain mesh (optional). Without it, the
    /// family uses the whole mesh.
    region: Option<String>,
    beta: Option<f64>,
    p32_target: Option<f64>,
    /// Necessary; only the mesh uses it
    aperture: ApertureDistribution,
    radius: RadiusSection,
    orientation: OrientationSection,
}

#[derive(Deserialize)]
#[serde(tag = "distribution", rename_all = "snake_case", deny_unknown_fields)]
enum RadiusSection {
    LogNormal {
        log_mean: f64,
        sd: f64,
        min: f64,
        max: f64,
        sampling: Option<LogNormalSampling>,
    },
    TruncatedPowerLaw {
        alpha: f64,
        min: f64,
        max: f64,
    },
    Exponential {
        mean: f64,
        min: f64,
        max: f64,
    },
    Constant {
        value: f64,
    },
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct OrientationSection {
    theta: Option<f64>,
    phi: Option<f64>,
    trend: Option<f64>,
    plunge: Option<f64>,
    dip: Option<f64>,
    strike: Option<f64>,
    kappa: f64,
}

/// Reads the TOML input file.
///
/// Returns the input options and the stochastic fracture families.
pub fn read_input(input_file: &str) -> Result<(Input, FractureFamilyCollection), DfngenError> {
    info!("DFN Generator Input File: {}", input_file);

    let text = std::fs::read_to_string(input_file).map_err(|source| DfngenError::ReadInput {
        file: input_file.to_string(),
        source,
    })?;

    let file = parse_input_file(&text)?;
    let input_folder = Path::new(input_file).parent().unwrap_or(Path::new(""));
    let domain_mesh = read_domain_mesh(&input_folder.join(&file.domain))?;
    let user_fractures = match &file.user_fractures {
        Some(fracture_file) => read_user_fractures(&input_folder.join(fracture_file))?,
        None => Vec::new(),
    };

    build_input(&file, domain_mesh, &user_fractures)
}

/// Names of the files which a TOML input file names
pub struct InputFileNames {
    /// Value of the `domain` key: the domain mesh (`.vtu`)
    pub domain: String,
    /// Value of the `user_fractures` key: the user-defined fractures (`.vtp`)
    pub user_fractures: Option<String>,
}

/// Reads the names of the files which a TOML input file names. It reads no file.
///
/// # Arguments
///
/// * `toml_text` - Text of the TOML input file
pub fn read_input_file_names(toml_text: &str) -> Result<InputFileNames, DfngenError> {
    let file = parse_input_file(toml_text)?;

    Ok(InputFileNames {
        domain: file.domain,
        user_fractures: file.user_fractures,
    })
}

/// Builds the input from the text of the input files. It reads no file.
///
/// Returns the input options and the stochastic fracture families, as `read_input` does.
///
/// # Arguments
///
/// * `toml_text` - Text of the TOML input file
/// * `domain_mesh_text` - Text of the file which the `domain` key names
/// * `user_fractures_text` - Text of the file which the `user_fractures` key
///   names. It is ignored if the input file does not have the key.
pub fn build_input_from_text(
    toml_text: &str,
    domain_mesh_text: &str,
    user_fractures_text: Option<&str>,
) -> Result<(Input, FractureFamilyCollection), DfngenError> {
    let file = parse_input_file(toml_text)?;
    let domain_mesh =
        parse_domain_mesh(domain_mesh_text).map_err(|message| DfngenError::InvalidFile {
            file: file.domain.clone(),
            message,
        })?;
    let user_fractures = match &file.user_fractures {
        Some(fracture_file) => {
            let text = user_fractures_text.ok_or_else(|| DfngenError::UserFractureFileMissing {
                file: fracture_file.clone(),
            })?;

            parse_user_fractures(text).map_err(|message| DfngenError::InvalidFile {
                file: fracture_file.clone(),
                message,
            })?
        }
        None => Vec::new(),
    };

    build_input(&file, domain_mesh, &user_fractures)
}

fn parse_input_file(text: &str) -> Result<InputFile, DfngenError> {
    toml::from_str(text).map_err(|e| DfngenError::InvalidToml(e.to_string()))
}

/// Position in the mesh regions (starts at 1) of the region name of a family
/// or of a user-defined fracture. No name gives 0, the whole mesh.
fn find_mesh_region_id(
    region_name: Option<&str>,
    mesh_regions: &[MeshRegion],
) -> Result<usize, RegionError> {
    let Some(region_name) = region_name else {
        return Ok(0);
    };

    match mesh_regions
        .iter()
        .position(|region| region.name == region_name)
    {
        Some(position) if mesh_regions[position].volume > 0. => Ok(position + 1),
        Some(_) => Err(RegionError::Empty {
            region: region_name.to_string(),
        }),
        None => Err(RegionError::Unknown {
            region: region_name.to_string(),
        }),
    }
}

/// Builds the input options and the stochastic fracture families.
///
/// # Arguments
///
/// * `file` - Content of the TOML input file
/// * `domain_mesh` - Domain mesh of `domain`
/// * `fracture_polygons` - User-defined fractures, as their file gives them
fn build_input(
    file: &InputFile,
    domain_mesh: DomainMesh,
    fracture_polygons: &[FracturePolygon],
) -> Result<(Input, FractureFamilyCollection), DfngenError> {
    // Each radius is at least 3 times min_feature_size, thus never 0
    if file.min_feature_size.is_nan() || file.min_feature_size <= 0. {
        return Err(DfngenError::MinFeatureSize);
    }

    let eps = file.min_feature_size * 1e-8;
    let domain_mesh = Rc::new(domain_mesh);

    let fram = file.fram.options();
    match &fram {
        None => info!("FRAM IS DISABLED"),
        Some(fram_options) if fram_options.relaxed => {
            info!("Running with relaxed FRAM. Mesh may not be fully conforming");
        }
        Some(_) => {}
    }

    if fram
        .as_ref()
        .is_some_and(|fram_options| !fram_options.check_user_fractures)
        && !fracture_polygons.is_empty()
    {
        info!("FRAM does not check two user-defined fractures which intersect");
    }

    let region_volumes: Vec<f64> = domain_mesh
        .regions()
        .iter()
        .map(|region| region.volume)
        .collect();

    debug!(
        "Input: min_feature_size = {}, region volumes {:?} m^3",
        file.min_feature_size, region_volumes
    );

    let stop_condition = file
        .generation
        .n_poly
        .map_or(StopCondition::P32, StopCondition::NPoly);

    let user_fractures = fracture_polygons
        .iter()
        .enumerate()
        .map(|(i, polygon)| {
            build_user_fracture(polygon, &domain_mesh).map_err(|source| DfngenError::UserFracture {
                fracture: i + 1,
                source,
            })
        })
        .collect::<Result<Vec<_>, _>>()?;
    check_piece_apertures(&user_fractures)?;

    let family_context = FamilyContext {
        stop_condition,
        min_feature_size: file.min_feature_size,
        eps,
        domain_mesh: &domain_mesh,
    };

    let mut families = Vec::with_capacity(file.families.len());
    let mut probabilities = Vec::with_capacity(file.families.len());
    for (i, section) in file.families.iter().enumerate() {
        let family = build_family(section, i + 1, &family_context).map_err(|source| {
            DfngenError::Family {
                family: i + 1,
                source,
            }
        })?;
        families.push(family);
        probabilities.push(section.probability);
    }

    let region_boundaries = file.region_boundaries.options();
    let mesh = build_mesh_options(file, region_boundaries.is_some())?;

    if mesh.is_some() && fram.is_none() {
        warn!("[mesh] with FRAM off: features smaller than min_feature_size give bad triangles");
    }

    let input = Input {
        stop_condition,
        domain_size: domain_mesh.size(),
        domain_volume: domain_mesh.volume(),
        domain_mesh,
        min_feature_size: file.min_feature_size,
        eps,
        fram,
        region_boundaries,
        mesh,
        keep_only_largest_cluster: file.clusters.keep_only_largest,
        keep_isolated_fractures: file.clusters.keep_isolated_fractures,
        force_large_fractures: file.generation.force_large_fractures,
        seed: file.seed,
        remove_fractures_less_than: file.generation.remove_fractures_less_than,
        region_volumes,
        user_fractures,
    };

    Ok((
        input,
        FractureFamilyCollection {
            families,
            probabilities,
        },
    ))
}

/// The mesh options, or `None` if `[mesh]` is off. The apertures then
/// have no effect. `build_family` has checked the family apertures.
///
/// # Arguments
///
/// * `file` - Content of the TOML input file
/// * `has_boundary_lines` - True if `[region_boundaries]` is on
fn build_mesh_options(
    file: &InputFile,
    has_boundary_lines: bool,
) -> Result<Option<MeshOptions>, DfngenError> {
    let section = &file.mesh;

    if !section.enable {
        return Ok(None);
    }

    if !has_boundary_lines {
        return Err(DfngenError::MeshNeedsRegionBoundaries);
    }
    if !(section.max_size_factor.is_finite() && section.max_size_factor >= 1.) {
        return Err(DfngenError::MeshMaxSizeFactor);
    }
    if !(section.slope.is_finite() && section.slope > 0.) {
        return Err(DfngenError::MeshSlope);
    }

    Ok(Some(MeshOptions {
        max_size_factor: section.max_size_factor,
        slope: section.slope,
        family_apertures: file
            .families
            .iter()
            .map(|family| family.aperture.clone())
            .collect(),
    }))
}

/// Moves a polygon of the user-defined fracture file to the coordinates of
/// the generation (the center of the mesh at the origin), finds its regions
/// in the mesh, and checks its aperture.
fn build_user_fracture(
    polygon: &FracturePolygon,
    domain_mesh: &DomainMesh,
) -> Result<UserFracture, UserFractureError> {
    if !(polygon.aperture.is_finite() && polygon.aperture > 0.) {
        return Err(UserFractureError::Aperture {
            aperture: polygon.aperture,
        });
    }

    let domain_center = domain_mesh.center();

    Ok(UserFracture {
        vertices: polygon
            .vertices
            .iter()
            .map(|&vertex| (Vector3::from(vertex) - domain_center).into())
            .collect(),
        region_ids: polygon
            .regions
            .iter()
            .map(|region| find_mesh_region_id(Some(region), domain_mesh.regions()))
            .collect::<Result<_, _>>()?,
        fracture_id: polygon.fracture_id,
        aperture: polygon.aperture,
    })
}

/// The pieces of a user-defined fracture (same `fracture_id`) must have the
/// same aperture: the mesh gives the fracture one aperture.
fn check_piece_apertures(user_fractures: &[UserFracture]) -> Result<(), DfngenError> {
    let mut first_apertures: HashMap<u64, f64> = HashMap::new();

    for fracture in user_fractures {
        let Some(fracture_id) = fracture.fracture_id else {
            continue;
        };
        let first = *first_apertures
            .entry(fracture_id)
            .or_insert(fracture.aperture);

        if fracture.aperture != first {
            return Err(DfngenError::UserFracturePieceApertures {
                fracture_id,
                first,
                other: fracture.aperture,
            });
        }
    }

    Ok(())
}

/// Values of the input file which each family needs.
struct FamilyContext<'a> {
    stop_condition: StopCondition,
    min_feature_size: f64,
    eps: f64,
    domain_mesh: &'a Rc<DomainMesh>,
}

fn build_family(
    section: &FamilySection,
    id: usize,
    context: &FamilyContext,
) -> Result<FractureFamily, FamilyError> {
    let region_id = find_mesh_region_id(section.region.as_deref(), context.domain_mesh.regions())?;

    let mut builder = FractureFamilyBuilder::default();

    builder
        .id(id)
        .number_of_points(valid_number_of_points(section.number_of_points)?)
        .radius(build_radius_distribution(
            &section.radius,
            3.0 * context.min_feature_size,
        ))
        .aspect_ratio(valid_aspect_ratio(section.aspect_ratio)?)
        .orientation(build_orientation(&section.orientation, context.eps)?)
        .region_id(region_id)
        .domain_mesh(context.domain_mesh.clone());

    section.aperture.validate()?;

    if let Some(beta) = section.beta {
        builder.beta(beta.to_radians());
    }

    match (context.stop_condition, section.p32_target) {
        (StopCondition::P32, Some(p32_target)) => {
            builder.p32_target(p32_target);
        }
        (StopCondition::P32, None) => {
            return Err(FamilyError::P32TargetMissing);
        }
        (StopCondition::NPoly(_), Some(_)) => {
            return Err(FamilyError::P32TargetNotAllowed);
        }
        (StopCondition::NPoly(_), None) => {}
    }

    builder.build()
}

/// Returns `radius`, or `min_radius` with a warning if `radius` is too small for the mesh.
fn limit_to_min_radius(radius: f64, min_radius: f64, description: &str) -> f64 {
    if radius >= min_radius {
        return radius;
    }

    warn!(
        "{} ({}) is less than 3 times min_feature_size, reset to {}",
        description, radius, min_radius
    );
    min_radius
}

fn build_radius_distribution(section: &RadiusSection, min_radius: f64) -> RadiusDistribution {
    match *section {
        RadiusSection::LogNormal {
            log_mean,
            sd,
            min,
            max,
            sampling,
        } => RadiusDistribution::new_log_normal(
            log_mean,
            sd,
            limit_to_min_radius(min, min_radius, "Minimum radius of log-normal distribution"),
            max,
            sampling.unwrap_or_default(),
        ),
        RadiusSection::TruncatedPowerLaw { alpha, min, max } => {
            RadiusDistribution::new_truncated_power_law(
                alpha,
                limit_to_min_radius(
                    min,
                    min_radius,
                    "Minimum radius of truncated power-law distribution",
                ),
                max,
            )
        }
        RadiusSection::Exponential { mean, min, max } => RadiusDistribution::new_exponential(
            1. / mean,
            limit_to_min_radius(
                min,
                min_radius,
                "Minimum radius of exponential distribution",
            ),
            max,
        ),
        RadiusSection::Constant { value } => RadiusDistribution::new_constant(limit_to_min_radius(
            value,
            min_radius,
            "Constant radius",
        )),
    }
}

/// The angles of the section are in degrees.
fn build_orientation(section: &OrientationSection, eps: f64) -> Result<Fisher, FamilyError> {
    let kappa = section.kappa;

    match (
        (section.theta, section.phi),
        (section.trend, section.plunge),
        (section.dip, section.strike),
    ) {
        ((Some(theta), Some(phi)), (None, None), (None, None)) => Ok(Fisher::new_with_theta_phi(
            theta.to_radians(),
            phi.to_radians(),
            kappa,
            eps,
        )),
        ((None, None), (Some(trend), Some(plunge)), (None, None)) => Ok(
            Fisher::new_with_trend_plunge(trend.to_radians(), plunge.to_radians(), kappa, eps),
        ),
        ((None, None), (None, None), (Some(dip), Some(strike))) => Ok(Fisher::new_with_dip_strike(
            dip.to_radians(),
            strike.to_radians(),
            kappa,
            eps,
        )),
        _ => Err(FamilyError::OrientationPair),
    }
}

fn valid_number_of_points(number_of_points: usize) -> Result<usize, FamilyError> {
    if number_of_points < 3 {
        return Err(FamilyError::NumberOfPoints);
    }
    Ok(number_of_points)
}

/// The fracture area must be greater than 0. A NaN aspect ratio is rejected too.
fn valid_aspect_ratio(aspect_ratio: f64) -> Result<f64, FamilyError> {
    if aspect_ratio.is_nan() || aspect_ratio <= 0. {
        return Err(FamilyError::AspectRatio);
    }
    Ok(aspect_ratio)
}

#[cfg(test)]
pub mod test_inputs {
    use super::*;
    use crate::domain_mesh::test_meshes::l_shaped_mesh;

    pub const GENERAL_OPTIONS: &str = r#"
        seed = 12345
        min_feature_size = 0.1
        domain = "domain.vtu"

        [generation]
        n_poly = 60
        force_large_fractures = false
        remove_fractures_less_than = 0.0

        [fram]
        enable = true
        relaxed = false
        triple_intersections = true
        write_rejected_fractures = false
        rejects_per_fracture = 10

        [region_boundaries]
        enable = false

        [mesh]
        enable = false
        max_size_factor = 10.0
        slope = 0.5

        [clusters]
        keep_only_largest = false
        keep_isolated_fractures = false
    "#;

    /// Reads the sections with the L-shaped mesh as the domain. The mesh has
    /// the regions `south` and `north`, and its center is (15, 25, 35).
    pub fn parse(sections: &[&str]) -> Result<(Input, FractureFamilyCollection), DfngenError> {
        build_input(&parse_input_file(&sections.concat())?, l_shaped_mesh(), &[])
    }
}

#[cfg(test)]
mod tests {
    use super::test_inputs::{parse, GENERAL_OPTIONS};
    use super::*;
    use crate::distribution::RadiusFunction;
    use crate::domain_mesh::test_meshes::{box_tetrahedra, l_shaped_mesh, region_names};

    const RECTANGLE_FAMILY: &str = r#"
        [[family]]
        number_of_points = 4
        probability = 0.6
        aspect_ratio = 1.5
        aperture = { distribution = "constant", value = 1e-4 }
        radius = { distribution = "log_normal", log_mean = 0.5, sd = 0.4, min = 1.0, max = 4.0 }
        orientation = { trend = 90.0, plunge = 30.0, kappa = 10.0 }
    "#;

    const ELLIPSE_FAMILY: &str = r#"
        [[family]]
        number_of_points = 8
        probability = 0.4
        aspect_ratio = 1.0
        beta = 90.0
        aperture = { distribution = "constant", value = 1e-4 }
        radius = { distribution = "truncated_power_law", alpha = 2.5, min = 1.0, max = 5.0 }
        orientation = { theta = 0.0, phi = 0.0, kappa = 8.0 }
    "#;

    fn error_message(result: Result<(Input, FractureFamilyCollection), DfngenError>) -> String {
        input_error(result).to_string()
    }

    fn input_error(result: Result<(Input, FractureFamilyCollection), DfngenError>) -> DfngenError {
        result.err().expect("the input must be rejected")
    }

    /// Path of a file of the test case `mesh_domain`
    fn mesh_domain_case_path(file_name: &str) -> String {
        format!(
            "{}/tests/data/mesh_domain/{}",
            env!("CARGO_MANIFEST_DIR"),
            file_name
        )
    }

    #[test]
    fn test_input_from_text_is_the_input_from_the_files() {
        let read_text =
            |file_name| std::fs::read_to_string(mesh_domain_case_path(file_name)).unwrap();
        let toml_text = read_text("mesh_domain.toml");
        let names = read_input_file_names(&toml_text).unwrap();
        let fracture_file = names.user_fractures.unwrap();

        let (from_text, families_from_text) = build_input_from_text(
            &toml_text,
            &read_text(&names.domain),
            Some(&read_text(&fracture_file)),
        )
        .unwrap();
        let (from_files, families_from_files) =
            read_input(&mesh_domain_case_path("mesh_domain.toml")).unwrap();

        let summary = |input: &Input, families: &FractureFamilyCollection| {
            (
                input.seed,
                input.domain_size,
                input.domain_volume,
                input.region_volumes.clone(),
                input.user_fractures.len(),
                families.families.len(),
            )
        };
        assert_eq!(
            summary(&from_text, &families_from_text),
            summary(&from_files, &families_from_files)
        );
    }

    #[test]
    fn test_input_from_text_without_the_fracture_text_is_rejected() {
        let read_text =
            |file_name| std::fs::read_to_string(mesh_domain_case_path(file_name)).unwrap();

        let error = input_error(build_input_from_text(
            &read_text("mesh_domain.toml"),
            &read_text("mesh_domain.vtu"),
            None,
        ));

        assert!(
            matches!(&error, DfngenError::UserFractureFileMissing { file } if file == "mesh_domain_user_fractures.vtp"),
            "{error}"
        );
    }

    #[test]
    fn test_general_options_are_read() {
        let (input, _) = parse(&[GENERAL_OPTIONS]).unwrap();

        assert_eq!(
            (input.seed, input.stop_condition, input.domain_size),
            (12345, StopCondition::NPoly(60), Vector3::new(10., 10., 10.))
        );
    }

    #[test]
    fn test_input_without_domain_is_rejected() {
        let options = GENERAL_OPTIONS.replace("domain = \"domain.vtu\"", "");

        let message = error_message(parse(&[&options]));

        assert!(message.contains("missing field `domain`"), "{message}");
    }

    #[test]
    fn test_input_with_domain_size_is_rejected() {
        let options = GENERAL_OPTIONS.replace(
            "domain = \"domain.vtu\"",
            "domain = \"domain.vtu\"\nsize = [10.0, 10.0, 10.0]",
        );

        let message = error_message(parse(&[&options]));

        assert!(message.contains("unknown field `size`"), "{message}");
    }

    #[test]
    fn test_mesh_gives_the_domain_size_and_volume() {
        let (input, _) = parse(&[GENERAL_OPTIONS]).unwrap();

        assert_eq!(
            (input.domain_size, input.domain_volume.round()),
            (Vector3::new(10., 10., 10.), 750.)
        );
    }

    /// `RECTANGLE_FAMILY` with the given value of the `region` key
    fn rectangle_family_with_region(region: &str) -> String {
        RECTANGLE_FAMILY.replace(
            "aspect_ratio = 1.5",
            &format!("aspect_ratio = 1.5\nregion = {region}"),
        )
    }

    #[test]
    fn test_family_region_selects_the_mesh_region_with_that_name() {
        let family = rectangle_family_with_region(r#""north""#);
        let (_, families) = parse(&[GENERAL_OPTIONS, &family]).unwrap();

        assert_eq!(families.families[0].region_id, 2);
    }

    #[test]
    fn test_family_without_region_uses_the_whole_mesh() {
        let (_, families) = parse(&[GENERAL_OPTIONS, RECTANGLE_FAMILY]).unwrap();

        assert_eq!(families.families[0].region_id, 0);
    }

    #[test]
    fn test_family_with_region_which_the_mesh_does_not_have_is_rejected() {
        let family = rectangle_family_with_region(r#""east""#);
        let error = input_error(parse(&[GENERAL_OPTIONS, &family]));

        assert!(
            matches!(
                &error,
                DfngenError::Family { family: 1, source: FamilyError::Region(RegionError::Unknown { region }) }
                    if region == "east"
            ),
            "{error}"
        );
    }

    #[test]
    fn test_family_with_a_number_as_region_is_rejected() {
        let family = rectangle_family_with_region("0");
        let message = error_message(parse(&[GENERAL_OPTIONS, &family]));

        assert!(message.contains("expected a string"), "{message}");
    }

    #[test]
    fn test_family_with_region_which_has_no_cell_is_rejected() {
        let tetrahedra = box_tetrahedra([0.; 3], [1.; 3], &[]);
        let mesh = DomainMesh::new(&tetrahedra, &region_names(&["empty"])).unwrap();
        let options = [GENERAL_OPTIONS, &rectangle_family_with_region(r#""empty""#)].concat();

        let error = input_error(build_input(&parse_input_file(&options).unwrap(), mesh, &[]));

        assert!(
            matches!(
                &error,
                DfngenError::Family { family: 1, source: FamilyError::Region(RegionError::Empty { region }) }
                    if region == "empty"
            ),
            "{error}"
        );
    }

    /// Builds the input with one triangle in the L-shaped mesh as the
    /// user-defined fracture. The first vertex is the center of the mesh.
    fn parse_with_triangle_of_regions(
        regions: &[&str],
    ) -> Result<(Input, FractureFamilyCollection), DfngenError> {
        parse_options_with_triangle(GENERAL_OPTIONS, regions)
    }

    /// As `parse_with_triangle_of_regions`, with other options
    fn parse_options_with_triangle(
        options: &str,
        regions: &[&str],
    ) -> Result<(Input, FractureFamilyCollection), DfngenError> {
        let triangle = FracturePolygon {
            regions: regions.iter().map(ToString::to_string).collect(),
            ..triangle(2e-4, Some(4))
        };

        build_input(
            &parse_input_file(options).unwrap(),
            l_shaped_mesh(),
            &[triangle],
        )
    }

    /// A triangle in the L-shaped mesh with no region; its first vertex is the center of the mesh
    fn triangle(aperture: f64, fracture_id: Option<u64>) -> FracturePolygon {
        FracturePolygon {
            vertices: vec![[15., 25., 35.], [16., 25., 35.], [16., 26., 35.]],
            aperture,
            regions: Vec::new(),
            fracture_id,
        }
    }

    /// Builds the input of `GENERAL_OPTIONS` with these user-defined fractures
    fn parse_with_triangles(
        triangles: &[FracturePolygon],
    ) -> Result<(Input, FractureFamilyCollection), DfngenError> {
        build_input(
            &parse_input_file(GENERAL_OPTIONS).unwrap(),
            l_shaped_mesh(),
            triangles,
        )
    }

    #[test]
    fn test_user_fracture_coordinates_are_moved_by_the_mesh_center() {
        let (input, _) = parse_with_triangle_of_regions(&[]).unwrap();

        assert_eq!(input.user_fractures[0].vertices[0], [0., 0., 0.]);
    }

    #[test]
    fn test_user_fracture_regions_select_the_mesh_regions_with_those_names() {
        let (input, _) = parse_with_triangle_of_regions(&["north", "south"]).unwrap();

        assert_eq!(input.user_fractures[0].region_ids, [2, 1]);
    }

    #[test]
    fn test_user_fracture_keeps_the_fracture_id_of_its_polygon() {
        let (input, _) = parse_with_triangle_of_regions(&[]).unwrap();

        assert_eq!(input.user_fractures[0].fracture_id, Some(4));
    }

    #[test]
    fn test_user_fracture_with_region_which_the_mesh_does_not_have_is_rejected() {
        let error = input_error(parse_with_triangle_of_regions(&["east"]));

        assert!(
            matches!(
                &error,
                DfngenError::UserFracture {
                    fracture: 1,
                    source: UserFractureError::Region(RegionError::Unknown { region })
                } if region == "east"
            ),
            "{error}"
        );
    }

    #[test]
    fn test_input_without_user_fractures_key_has_no_user_fractures() {
        let (input, _) = parse(&[GENERAL_OPTIONS]).unwrap();

        assert!(input.user_fractures.is_empty());
    }

    #[test]
    fn test_user_fracture_table_is_rejected() {
        let message = error_message(parse(&[
            GENERAL_OPTIONS,
            "[[user_fracture]]\nvertices = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]]",
        ]));

        assert!(
            message.contains("unknown field `user_fracture`"),
            "{message}"
        );
    }

    #[test]
    fn test_eps_is_derived_from_h() {
        let (input, _) = parse(&[GENERAL_OPTIONS]).unwrap();

        assert_eq!(input.eps, 0.1 * 1e-8);
    }

    #[test]
    fn test_family_ids_follow_the_file_order() {
        let (_, families) = parse(&[GENERAL_OPTIONS, RECTANGLE_FAMILY, ELLIPSE_FAMILY]).unwrap();

        let ids: Vec<usize> = families.families.iter().map(|family| family.id).collect();
        assert_eq!((ids, families.probabilities), (vec![1, 2], vec![0.6, 0.4]));
    }

    #[test]
    fn test_family_beta_is_converted_to_radians() {
        let (_, families) = parse(&[GENERAL_OPTIONS, ELLIPSE_FAMILY]).unwrap();

        assert_eq!(families.families[0].beta, Some(std::f64::consts::FRAC_PI_2));
    }

    #[test]
    fn test_family_without_beta_has_uniform_rotation() {
        let (_, families) = parse(&[GENERAL_OPTIONS, RECTANGLE_FAMILY]).unwrap();

        assert_eq!(families.families[0].beta, None);
    }

    #[test]
    fn test_minimum_radius_is_limited_to_three_times_h() {
        let family = ELLIPSE_FAMILY.replace("min = 1.0", "min = 0.1");

        let (_, families) = parse(&[GENERAL_OPTIONS, &family]).unwrap();

        assert_eq!(families.families[0].radius.min, 3.0 * 0.1);
    }

    #[test]
    fn test_missing_key_is_rejected() {
        let options = GENERAL_OPTIONS.replace("seed = 12345", "");

        let message = error_message(parse(&[&options]));

        assert!(message.contains("missing field `seed`"), "{}", message);
    }

    /// The options with this text in place of the section `[name]`
    fn replace_section(options: &str, name: &str, section: &str) -> String {
        let start = options.find(&format!("[{name}]")).unwrap();
        let end = start + 1 + options[start + 1..].find('[').unwrap();
        format!("{}{section}\n{}", &options[..start], &options[end..])
    }

    /// `GENERAL_OPTIONS` with this text in place of its `[fram]` section
    fn replace_fram_section(fram_section: &str) -> String {
        replace_section(GENERAL_OPTIONS, "fram", fram_section)
    }

    #[test]
    fn test_print_reject_reasons_key_is_rejected() {
        let options = replace_fram_section("[fram]\nenable = true\nprint_reject_reasons = true\n");

        let message = error_message(parse(&[&options]));

        assert!(
            message.contains("unknown field `print_reject_reasons`"),
            "{message}"
        );
    }

    #[test]
    fn test_fram_is_off_without_its_section() {
        let (input, _) = parse(&[&replace_fram_section("")]).unwrap();

        assert_eq!(input.fram, None);
    }

    #[test]
    fn test_fram_section_with_only_enable_has_the_defaults() {
        let options = replace_fram_section("[fram]\nenable = true\n");

        let (input, _) = parse(&[&options]).unwrap();

        assert_eq!(
            input.fram,
            Some(FramOptions {
                relaxed: false,
                triple_intersections: false,
                check_user_fractures: false,
                write_rejected_fractures: false,
                rejects_per_fracture: 10,
            })
        );
    }

    #[test]
    fn test_fram_is_off_with_enable_false_and_other_keys() {
        let options = GENERAL_OPTIONS.replace("enable = true", "enable = false");

        let (input, _) = parse(&[&options]).unwrap();

        assert_eq!(input.fram, None);
    }

    /// `GENERAL_OPTIONS` with this text in place of its `[region_boundaries]` section
    fn replace_region_boundaries_section(section: &str) -> String {
        replace_section(GENERAL_OPTIONS, "region_boundaries", section)
    }

    #[test]
    fn test_input_without_region_boundaries_section_is_rejected() {
        let options = replace_region_boundaries_section("");

        let message = error_message(parse(&[&options]));

        assert!(
            message.contains("missing field `region_boundaries`"),
            "{message}"
        );
    }

    #[test]
    fn test_region_boundaries_without_enable_is_rejected() {
        let options = replace_region_boundaries_section("[region_boundaries]\n");

        let message = error_message(parse(&[&options]));

        assert!(message.contains("missing field `enable`"), "{message}");
    }

    #[test]
    fn test_region_boundaries_are_on_with_enable_true() {
        let options = replace_region_boundaries_section("[region_boundaries]\nenable = true\n");

        let (input, _) = parse(&[&options]).unwrap();

        assert_eq!(input.region_boundaries, Some(RegionBoundaryOptions));
    }

    #[test]
    fn test_region_boundaries_are_off_with_enable_false() {
        let (input, _) = parse(&[GENERAL_OPTIONS]).unwrap();

        assert_eq!(input.region_boundaries, None);
    }

    #[test]
    fn test_unknown_region_boundaries_key_is_rejected() {
        let options =
            replace_region_boundaries_section("[region_boundaries]\nenable = true\ncut = true\n");

        let message = error_message(parse(&[&options]));

        assert!(message.contains("unknown field `cut`"), "{message}");
    }

    #[test]
    fn test_user_fracture_check_is_on_with_its_key_true() {
        let options = replace_fram_section("[fram]\nenable = true\ncheck_user_fractures = true\n");

        let (input, _) = parse(&[&options]).unwrap();

        assert!(input.fram.unwrap().check_user_fractures);
    }

    #[test]
    fn test_unknown_key_is_rejected() {
        let options = GENERAL_OPTIONS.replace("seed = 12345", "seed = 12345\nsede = 1");

        let message = error_message(parse(&[&options]));

        assert!(message.contains("unknown field `sede`"), "{}", message);
    }

    #[test]
    fn test_unknown_radius_distribution_is_rejected() {
        let family = ELLIPSE_FAMILY.replace("truncated_power_law", "power_law");

        let message = error_message(parse(&[GENERAL_OPTIONS, &family]));

        assert!(
            message.contains("unknown variant `power_law`"),
            "{}",
            message
        );
    }

    fn log_normal_sampling(family: &str) -> LogNormalSampling {
        let (_, families) = parse(&[GENERAL_OPTIONS, family]).unwrap();
        match families.families[0].radius.function {
            RadiusFunction::LogNormal { sampling, .. } => sampling,
            _ => panic!("the family must have a log-normal radius"),
        }
    }

    #[test]
    fn test_input_with_truncate_fractures_is_rejected() {
        let options = GENERAL_OPTIONS.replace(
            "domain = \"domain.vtu\"",
            "domain = \"domain.vtu\"\ntruncate_fractures = false",
        );

        let message = error_message(parse(&[&options]));

        assert!(
            message.contains("unknown field `truncate_fractures`"),
            "{message}"
        );
    }

    #[test]
    fn test_log_normal_sampling_is_rejection_without_the_key() {
        assert_eq!(
            log_normal_sampling(RECTANGLE_FAMILY),
            LogNormalSampling::Rejection
        );
    }

    #[test]
    fn test_log_normal_sampling_inverse_cdf_is_read() {
        let family = RECTANGLE_FAMILY.replace("max = 4.0", "max = 4.0, sampling = \"inverse_cdf\"");

        assert_eq!(log_normal_sampling(&family), LogNormalSampling::InverseCdf);
    }

    #[test]
    fn test_unknown_log_normal_sampling_is_rejected() {
        let family = RECTANGLE_FAMILY.replace("max = 4.0", "max = 4.0, sampling = \"ziggurat\"");

        let message = error_message(parse(&[GENERAL_OPTIONS, &family]));

        assert!(
            message.contains("unknown variant `ziggurat`"),
            "{}",
            message
        );
    }

    #[test]
    fn test_sampling_key_of_other_radius_distribution_is_rejected() {
        let family = ELLIPSE_FAMILY.replace("max = 5.0", "max = 5.0, sampling = \"inverse_cdf\"");

        let message = error_message(parse(&[GENERAL_OPTIONS, &family]));

        assert!(message.contains("unknown field `sampling`"), "{}", message);
    }

    #[test]
    fn test_parameter_of_other_radius_distribution_is_rejected() {
        let family = ELLIPSE_FAMILY.replace("alpha = 2.5", "alpha = 2.5, sd = 1.0");

        let message = error_message(parse(&[GENERAL_OPTIONS, &family]));

        assert!(message.contains("unknown field `sd`"), "{}", message);
    }

    #[test]
    fn test_orientation_with_two_angle_pairs_is_rejected() {
        let family = ELLIPSE_FAMILY.replace("theta = 0.0", "theta = 0.0, dip = 1.0, strike = 1.0");

        let error = input_error(parse(&[GENERAL_OPTIONS, &family]));

        assert!(
            matches!(
                error,
                DfngenError::Family {
                    family: 1,
                    source: FamilyError::OrientationPair
                }
            ),
            "{error}"
        );
    }

    #[test]
    fn test_input_without_n_poly_has_p32_stop_condition() {
        let options = GENERAL_OPTIONS.replace("n_poly = 60", "");
        let family = ELLIPSE_FAMILY.replace("beta = 90.0", "beta = 90.0\np32_target = 0.5");

        let (input, _) = parse(&[&options, &family]).unwrap();

        assert_eq!(input.stop_condition, StopCondition::P32);
    }

    #[test]
    fn test_p32_stop_condition_needs_family_target() {
        let options = GENERAL_OPTIONS.replace("n_poly = 60", "");

        let error = input_error(parse(&[&options, ELLIPSE_FAMILY]));

        assert!(
            matches!(
                error,
                DfngenError::Family {
                    family: 1,
                    source: FamilyError::P32TargetMissing
                }
            ),
            "{error}"
        );
    }

    #[test]
    fn test_family_target_with_n_poly_is_rejected() {
        let family = ELLIPSE_FAMILY.replace("beta = 90.0", "beta = 90.0\np32_target = 0.5");

        let error = input_error(parse(&[GENERAL_OPTIONS, &family]));

        assert!(
            matches!(
                error,
                DfngenError::Family {
                    family: 1,
                    source: FamilyError::P32TargetNotAllowed
                }
            ),
            "{error}"
        );
    }

    #[test]
    fn test_stop_condition_key_is_rejected() {
        let options =
            GENERAL_OPTIONS.replace("n_poly = 60", "n_poly = 60\nstop_condition = \"n_poly\"");

        let message = error_message(parse(&[&options]));

        assert!(
            message.contains("unknown field `stop_condition`"),
            "{}",
            message
        );
    }

    #[test]
    fn test_output_section_is_rejected() {
        let options = format!("{GENERAL_OPTIONS}\n[output]\nprint_reject_reasons = false\n");

        let message = error_message(parse(&[&options]));

        assert!(message.contains("unknown field `output`"), "{message}");
    }

    #[test]
    fn test_rejects_per_fracture_in_generation_section_is_rejected() {
        let options =
            GENERAL_OPTIONS.replace("n_poly = 60", "n_poly = 60\nrejects_per_fracture = 10");

        let message = error_message(parse(&[&options]));

        assert!(
            message.contains("unknown field `rejects_per_fracture`"),
            "{message}"
        );
    }

    #[test]
    fn test_radii_list_increase_key_is_rejected() {
        let options =
            GENERAL_OPTIONS.replace("n_poly = 60", "n_poly = 60\nradii_list_increase = 0.1");

        let message = error_message(parse(&[&options]));

        assert!(
            message.contains("unknown field `radii_list_increase`"),
            "{message}"
        );
    }

    #[test]
    fn test_angle_unit_is_rejected() {
        let options = GENERAL_OPTIONS.replace(
            "min_feature_size = 0.1",
            "min_feature_size = 0.1\nangle_unit = \"radian\"",
        );

        let message = error_message(parse(&[&options]));

        assert!(
            message.contains("unknown field `angle_unit`"),
            "{}",
            message
        );
    }

    #[test]
    fn test_family_with_two_points_is_rejected() {
        let family = ELLIPSE_FAMILY.replace("number_of_points = 8", "number_of_points = 2");

        let error = input_error(parse(&[GENERAL_OPTIONS, &family]));

        assert!(
            matches!(
                error,
                DfngenError::Family {
                    family: 1,
                    source: FamilyError::NumberOfPoints
                }
            ),
            "{error}"
        );
    }

    #[test]
    fn test_family_with_zero_aspect_ratio_is_rejected() {
        let family = ELLIPSE_FAMILY.replace("aspect_ratio = 1.0", "aspect_ratio = 0.0");

        let error = input_error(parse(&[GENERAL_OPTIONS, &family]));

        assert!(
            matches!(
                error,
                DfngenError::Family {
                    family: 1,
                    source: FamilyError::AspectRatio
                }
            ),
            "{error}"
        );
    }

    /// `GENERAL_OPTIONS` with the boundary lines on and this `[mesh]` section
    fn options_with_mesh_section(mesh_section: &str) -> String {
        let options = replace_section(
            GENERAL_OPTIONS,
            "region_boundaries",
            "[region_boundaries]\nenable = true\n",
        );
        replace_section(&options, "mesh", mesh_section)
    }

    #[test]
    fn test_mesh_is_off_with_enable_false() {
        let (input, _) = parse(&[GENERAL_OPTIONS, ELLIPSE_FAMILY]).unwrap();

        assert_eq!(input.mesh, None);
    }

    #[test]
    fn test_mesh_options_are_read_with_enable_true() {
        let options = options_with_mesh_section(
            "[mesh]\nenable = true\nmax_size_factor = 4.0\nslope = 0.3\n",
        );

        let (input, _) = parse(&[&options, ELLIPSE_FAMILY]).unwrap();

        assert_eq!(
            input.mesh,
            Some(MeshOptions {
                max_size_factor: 4.,
                slope: 0.3,
                family_apertures: vec![ApertureDistribution::Constant { value: 1e-4 }],
            })
        );
    }

    #[test]
    fn test_input_without_mesh_section_is_rejected() {
        let options = replace_section(GENERAL_OPTIONS, "mesh", "");

        let message = error_message(parse(&[&options]));

        assert!(message.contains("missing field `mesh`"), "{message}");
    }

    #[test]
    fn test_mesh_section_without_a_key_is_rejected() {
        for key in ["enable", "max_size_factor", "slope"] {
            let mesh_section: String =
                "[mesh]\nenable = false\nmax_size_factor = 10.0\nslope = 0.5\n"
                    .lines()
                    .filter(|line| !line.starts_with(key))
                    .map(|line| format!("{line}\n"))
                    .collect();
            let options = replace_section(GENERAL_OPTIONS, "mesh", &mesh_section);

            let message = error_message(parse(&[&options]));

            assert!(
                message.contains(&format!("missing field `{key}`")),
                "{message}"
            );
        }
    }

    #[test]
    fn test_mesh_needs_region_boundaries() {
        let options = replace_section(
            GENERAL_OPTIONS,
            "mesh",
            "[mesh]\nenable = true\nmax_size_factor = 10.0\nslope = 0.5\n",
        );

        let error = input_error(parse(&[&options, ELLIPSE_FAMILY]));

        assert!(
            matches!(error, DfngenError::MeshNeedsRegionBoundaries),
            "{error}"
        );
    }

    #[test]
    fn test_family_without_aperture_is_rejected() {
        let family = ELLIPSE_FAMILY.replace(
            "aperture = { distribution = \"constant\", value = 1e-4 }",
            "",
        );

        let message = error_message(parse(&[GENERAL_OPTIONS, &family]));

        assert!(message.contains("missing field `aperture`"), "{message}");
    }

    #[test]
    fn test_negative_aperture_sd_is_rejected_without_mesh() {
        let family = ELLIPSE_FAMILY.replace(
            "{ distribution = \"constant\", value = 1e-4 }",
            "{ distribution = \"log_normal\", log_mean = -8.0, sd = -1.0 }",
        );

        let error = input_error(parse(&[GENERAL_OPTIONS, &family]));

        assert!(
            matches!(
                error,
                DfngenError::Family {
                    family: 1,
                    source: FamilyError::Aperture(_)
                }
            ),
            "{error}"
        );
    }

    #[test]
    fn test_mesh_slope_of_zero_is_rejected() {
        let options = options_with_mesh_section(
            "[mesh]\nenable = true\nmax_size_factor = 10.0\nslope = 0.0\n",
        );

        let error = input_error(parse(&[&options, ELLIPSE_FAMILY]));

        assert!(matches!(error, DfngenError::MeshSlope), "{error}");
    }

    #[test]
    fn test_mesh_max_size_factor_below_one_is_rejected() {
        let options = options_with_mesh_section(
            "[mesh]\nenable = true\nmax_size_factor = 0.5\nslope = 0.5\n",
        );

        let error = input_error(parse(&[&options, ELLIPSE_FAMILY]));

        assert!(matches!(error, DfngenError::MeshMaxSizeFactor), "{error}");
    }

    #[test]
    fn test_user_fracture_keeps_the_aperture_of_its_polygon() {
        let (input, _) = parse_with_triangles(&[triangle(2e-4, None)]).unwrap();

        assert_eq!(input.user_fractures[0].aperture, 2e-4);
    }

    #[test]
    fn test_zero_and_nan_user_fracture_apertures_are_rejected() {
        for aperture in [0., f64::NAN] {
            let error = input_error(parse_with_triangles(&[triangle(aperture, None)]));

            assert!(
                matches!(
                    error,
                    DfngenError::UserFracture {
                        fracture: 1,
                        source: UserFractureError::Aperture { .. }
                    }
                ),
                "{error}"
            );
        }
    }

    #[test]
    fn test_pieces_with_different_apertures_are_rejected() {
        let error = input_error(parse_with_triangles(&[
            triangle(2e-4, Some(4)),
            triangle(2e-4, Some(5)),
            triangle(3e-4, Some(4)),
        ]));

        assert!(
            matches!(
                error,
                DfngenError::UserFracturePieceApertures {
                    fracture_id: 4,
                    first: 2e-4,
                    other: 3e-4
                }
            ),
            "{error}"
        );
    }

    #[test]
    fn test_zero_min_feature_size_is_rejected() {
        let options = GENERAL_OPTIONS.replace("min_feature_size = 0.1", "min_feature_size = 0.0");

        let error = input_error(parse(&[&options]));

        assert!(matches!(error, DfngenError::MinFeatureSize), "{error}");
    }
}
