use std::path::Path;
use std::rc::Rc;

use parry3d_f64::na::{Point2, Vector3};
use serde::Deserialize;
use tracing::{debug, info, warn};

use crate::cg::{normal_from_dip_strike, normal_from_trend_plunge};
use crate::distribution::{Fisher, LogNormalSampling};
use crate::domain_mesh::{DomainMesh, MeshRegion};
use crate::error::DfngenError;
use crate::fracture::fracture_family::{
    FractureFamily, FractureFamilyBuilder, FractureFamilyCollection,
};
use crate::fracture::poly::DomainBoundaryMode;
use crate::fracture::user_defined_fractures::{FractureDef, FractureNodesDef, UserFracture};
use crate::io::vtu::read_domain_mesh;
use crate::structures::RadiusDistribution;

#[derive(Deserialize, Default, Debug, Clone, Copy, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum StopCondition {
    /// Stop at a fixed number of accepted fractures.
    #[default]
    NPoly,

    /// Stop when each family has its target fracture intensity.
    P32,
}

#[derive(Default, Debug)]
pub struct Input {
    /// DFN generation stop condition.
    pub stop_condition: StopCondition,

    /// Number of polygons to place in the DFN with the `NPoly` stop condition.
    pub n_poly: usize,

    /// Domain size with dimension x*y*z for DFN, centered at the origin.
    /// With a domain mesh, this is the size of its bounding box.
    pub domain_size: Vector3<f64>,

    /// Volume of the domain. With a domain mesh, this is the sum of its cell volumes.
    pub domain_volume: f64,

    /// The cells which make the domain, if the input file gives a mesh.
    pub domain_mesh: Option<Rc<DomainMesh>>,

    /// What to do with the part of a fracture which is outside the domain.
    pub domain_boundary_mode: DomainBoundaryMode,

    /// Minimum feature size, FRAM parameter.
    pub h: f64,

    /// epsilon
    pub eps: f64,

    /// Percent to increase the size of the pre-generated radii lists, per family.
    /// Example: 0.2 will increase the size of the list by %20. See example input files
    /// for more details.
    pub radii_list_increase: f64,

    /// This option disables the FRAM algorithm. There will be no
    /// fracture rejections.
    pub disable_fram: bool,

    /// This option uses a relaxed version of the FRAM algorithm. The mesh may not
    /// be perfectly conforming
    pub r_fram: bool,

    /// Accept or reject triple intersections
    ///     False - Off (Reject)
    ///     True  - On  (Accept)
    pub triple_intersections: bool,

    /// DFN will only keep clusters with connections to
    /// domain boundaries which are set to 1:
    ///
    /// boundaryFaces[0] = +X domain boundary
    /// boundaryFaces[1] = -X domain boundary
    /// boundaryFaces[2] = +Y domain boundary
    /// boundaryFaces[3] = -Y domain boundary
    /// boundaryFaces[4] = +Z domain boundary
    /// boundaryFaces[5] = -Z domain boundary
    pub boundary_faces: [bool; 6],

    /// 0 - Keep any clusters which connects the specified
    ///     boundary faces in boundaryFaces option below
    /// 1 - Keep only the largest cluster which connects
    ///     the specified boundary faces in boundaryFaces option below.
    ///
    /// If ignoreBoundaryFaces is also set to 1, DFNGen will keep the largest
    /// cluster which connects at least any two sides of the domain.
    pub keep_only_largest_cluster: bool,

    /// 0 - remove isolated fractures and clusters
    /// 1 - Keep isolated fractures and clusters
    pub keep_isolated_fractures: bool,

    /// Useful for debugging,
    /// This option will print all fracture rejection reasons as they occur.
    ///     0 - Disable
    ///     1 - Print all rejection reasons to screen
    pub print_reject_reasons: bool,

    /// Inserts the largest possible fracture for each defined fracture family,
    /// defined by the user-defined maxium radius
    ///     0 - Off (Do not force insertion of larest fractures)
    ///     1 - On  (Force insertion of largest fractures)
    pub force_large_fractures: bool,

    /// Seed for random generator.
    pub seed: u64,

    ///    0 - Ignore this option, keep all fractures.
    ///
    /// (>0) - Size of minimum fracture radius. Fractures smaller than
    ///        defined radius will be removed AFTER DFN generation.
    ///
    ///        Minimum and maximum size options under fracture family
    ///        distributions will still be used while generating the DFN.
    pub remove_fractures_less_than: f64,

    /// Caution: Can create a very large file.
    /// True: Write the rejected fractures to fractures_rejected.vtp.
    pub output_rejected_fractures: bool,

    /// If a fracture is rejected, it will be re-translated
    /// to a new position this number of times.
    ///
    /// This helps hit distribution targets for stochastic families
    /// families (Set to 1 to ignore this feature)
    pub rejects_per_fracture: usize,

    /// Array of layers:
    /// e.g. {+z1, -z1, +z2, -z2, ... , +zn, -zn}
    pub layers: Vec<f64>,

    /// Array of volumes for each defined layer, in the same order
    /// which layers were listed.
    pub layer_volumes: Vec<f64>,

    /// Array of regions:
    /// e.g. {+z1, -z1, +z2, -z2, ... , +zn, -zn}
    pub regions: Vec<f64>,

    /// Array of volumes for each defined layer, in the same order
    /// which regions were listed.
    pub region_volumes: Vec<f64>,

    /// Vector of points defining the 2D boundary of the domain polygon.
    /// Empty if the domain has no polygon boundary.
    pub domain_vertices: Vec<Point2<f64>>,

    /// False - Use boundaryFaces option.
    /// True  - Ignore boundaryFaces option, keep all clusters
    ///         and remove fractures with no intersections
    pub ignore_boundary_faces: bool,

    /// User-defined fractures, in the order of the input file.
    pub user_fractures: Vec<UserFracture>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct InputFile {
    seed: u64,
    h: f64,
    domain: DomainSection,
    generation: GenerationSection,
    fram: FramSection,
    clusters: ClustersSection,
    output: OutputSection,
    #[serde(default, rename = "family")]
    families: Vec<FamilySection>,
    #[serde(default, rename = "user_fracture")]
    user_fractures: Vec<UserFractureSection>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct DomainSection {
    /// Necessary if there is no `mesh`
    size: Option<[f64; 3]>,
    /// Path of an ASCII `.vtu` file which gives the domain, in place of `size`.
    /// A relative path starts at the folder of the input file.
    mesh: Option<String>,
    /// Increase of the domain size for the insertion of fracture centers.
    /// Each side moves out by half of the value. Fractures are truncated to `size`.
    size_increase: [f64; 3],
    /// True: cut each fracture at the domain boundary. False: keep it whole.
    truncate_fractures: bool,
    /// {-z, +z} of each layer
    #[serde(default)]
    layers: Vec<[f64; 2]>,
    /// {-x, +x, -y, +y, -z, +z} of each region
    #[serde(default)]
    regions: Vec<[f64; 6]>,
    /// {x, y} of each vertex of the domain polygon
    #[serde(default)]
    polygon_boundary: Vec<[f64; 2]>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct GenerationSection {
    stop_condition: StopCondition,
    n_poly: Option<usize>,
    radii_list_increase: f64,
    rejects_per_fracture: usize,
    force_large_fractures: bool,
    remove_fractures_less_than: f64,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct FramSection {
    disable: bool,
    relaxed: bool,
    triple_intersections: bool,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct ClustersSection {
    boundary_faces: [bool; 6],
    ignore_boundary_faces: bool,
    keep_only_largest: bool,
    keep_isolated_fractures: bool,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct OutputSection {
    print_reject_reasons: bool,
    rejected_fractures: bool,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct FamilySection {
    number_of_points: usize,
    probability: f64,
    aspect_ratio: f64,
    layer: usize,
    region: usize,
    beta: Option<f64>,
    p32_target: Option<f64>,
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

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct UserFractureSection {
    // Definition by parameters
    number_of_points: Option<usize>,
    radius: Option<f64>,
    aspect_ratio: Option<f64>,
    beta: Option<f64>,
    translation: Option<[f64; 3]>,
    normal: Option<[f64; 3]>,
    trend_plunge: Option<[f64; 2]>,
    dip_strike: Option<[f64; 2]>,

    // Definition by vertices
    vertices: Option<Vec<[f64; 3]>>,
}

fn invalid_input(message: String) -> DfngenError {
    DfngenError::InvalidInput(message)
}

/// Reads the TOML input file.
///
/// Returns the input options and the stochastic fracture families.
pub fn read_input(input_file: &str) -> Result<(Input, FractureFamilyCollection), DfngenError> {
    info!("DFN Generator Input File: {}", input_file);

    let text = std::fs::read_to_string(input_file)
        .map_err(|e| invalid_input(format!("cannot read {}: {}", input_file, e)))?;

    let input_folder = Path::new(input_file).parent().unwrap_or(Path::new(""));

    parse_input(&text, input_folder)
}

fn parse_input(
    text: &str,
    input_folder: &Path,
) -> Result<(Input, FractureFamilyCollection), DfngenError> {
    let file: InputFile = toml::from_str(text).map_err(|e| invalid_input(e.to_string()))?;

    let domain_mesh = match &file.domain.mesh {
        Some(mesh_file) => Some(read_domain_mesh(&input_folder.join(mesh_file))?),
        None => None,
    };

    build_input(file, domain_mesh)
}

/// The size, the center and the volume of the domain, from `domain.size` or
/// from the mesh.
fn build_domain_box(
    domain: &DomainSection,
    domain_mesh: Option<&DomainMesh>,
) -> Result<(Vector3<f64>, Vector3<f64>, f64), String> {
    match (domain.size, domain_mesh) {
        (Some(size), None) => Ok((
            Vector3::from(size),
            Vector3::zeros(),
            size[0] * size[1] * size[2],
        )),
        (None, Some(mesh)) => {
            if !domain.layers.is_empty()
                || !domain.regions.is_empty()
                || !domain.polygon_boundary.is_empty()
            {
                return Err(
                    "domain.layers, domain.regions and domain.polygon_boundary cannot be used with domain.mesh. The cell property \"region\" of the mesh gives the regions."
                        .to_string(),
                );
            }
            Ok((mesh.size(), mesh.center(), mesh.volume()))
        }
        (None, None) => Err("domain.size or domain.mesh is necessary".to_string()),
        (Some(_), Some(_)) => {
            Err("domain.size cannot be used with domain.mesh. The mesh gives the size.".to_string())
        }
    }
}

/// Position in the mesh regions (starts at 1) of the region value of a family.
/// 0 is the whole mesh.
fn find_mesh_region_id(region_value: usize, mesh_regions: &[MeshRegion]) -> Result<usize, String> {
    if region_value == 0 {
        return Ok(0);
    }

    mesh_regions
        .iter()
        .position(|region| region.value == region_value as u64)
        .map(|position| position + 1)
        .ok_or(format!(
            "no cell of the domain mesh has the region {}",
            region_value
        ))
}

fn build_input(
    file: InputFile,
    domain_mesh: Option<DomainMesh>,
) -> Result<(Input, FractureFamilyCollection), DfngenError> {
    let eps = file.h * 1e-8;
    let domain_mesh = domain_mesh.map(Rc::new);

    if file.fram.disable {
        info!("FRAM IS DISABLED");
    }

    if file.fram.relaxed {
        info!("Running with relaxed FRAM. Mesh may not be fully conforming");
    }

    let (domain_size, domain_center, domain_volume) =
        build_domain_box(&file.domain, domain_mesh.as_deref()).map_err(invalid_input)?;
    let domain_size_increase = Vector3::from(file.domain.size_increase);

    let layers: Vec<f64> = file.domain.layers.iter().flatten().copied().collect();
    let layer_volumes: Vec<f64> = file
        .domain
        .layers
        .iter()
        .map(|layer| domain_size[0] * domain_size[1] * (layer[1] - layer[0]).abs())
        .collect();

    // {-x, +x, -y, +y, -z, +z} of each region. With a domain mesh, this is the
    // bounding box of the cells of the region.
    let (regions, region_volumes): (Vec<f64>, Vec<f64>) = match &domain_mesh {
        Some(mesh) => (
            mesh.regions()
                .iter()
                .flat_map(|region| {
                    let (mins, maxs) = (region.bounds.mins, region.bounds.maxs);
                    [mins.x, maxs.x, mins.y, maxs.y, mins.z, maxs.z]
                })
                .collect(),
            mesh.regions().iter().map(|region| region.volume).collect(),
        ),
        None => (
            file.domain.regions.iter().flatten().copied().collect(),
            file.domain
                .regions
                .iter()
                .map(|region| {
                    (region[1] - region[0]).abs()
                        * (region[3] - region[2]).abs()
                        * (region[5] - region[4]).abs()
                })
                .collect(),
        ),
    };

    let domain_vertices: Vec<Point2<f64>> = file
        .domain
        .polygon_boundary
        .iter()
        .map(|&[x, y]| Point2::new(x, y))
        .collect();

    debug!(
        "Input: h = {}, layer volumes {:?} m^3, region volumes {:?} m^3, {} polygon boundary vertices",
        file.h,
        layer_volumes,
        region_volumes,
        domain_vertices.len()
    );

    let stop_condition = file.generation.stop_condition;
    let n_poly = match (stop_condition, file.generation.n_poly) {
        (StopCondition::NPoly, None) if !file.families.is_empty() => {
            return Err(invalid_input(
                "generation.n_poly is necessary for stop_condition \"n_poly\"".to_string(),
            ));
        }
        (_, n_poly) => n_poly.unwrap_or(0),
    };

    let user_fractures = file
        .user_fractures
        .iter()
        .enumerate()
        .map(|(i, section)| {
            build_user_fracture(section, eps, &domain_center)
                .map_err(|message| invalid_input(format!("user_fracture {}: {}", i + 1, message)))
        })
        .collect::<Result<Vec<_>, _>>()?;

    let family_context = FamilyContext {
        stop_condition,
        h: file.h,
        eps,
        layers: &layers,
        regions: &regions,
        domain_mesh: domain_mesh.as_ref(),
        domain_size: &domain_size,
        domain_size_increase: &domain_size_increase,
    };

    let mut families = Vec::with_capacity(file.families.len());
    let mut probabilities = Vec::with_capacity(file.families.len());
    for (i, section) in file.families.iter().enumerate() {
        let family = build_family(section, i + 1, &family_context)
            .map_err(|message| invalid_input(format!("family {}: {}", i + 1, message)))?;
        families.push(family);
        probabilities.push(section.probability);
    }

    let input = Input {
        stop_condition,
        n_poly,
        domain_size,
        domain_volume,
        domain_mesh,
        domain_boundary_mode: if file.domain.truncate_fractures {
            DomainBoundaryMode::Truncate
        } else {
            DomainBoundaryMode::KeepWhole
        },
        h: file.h,
        eps,
        radii_list_increase: file.generation.radii_list_increase,
        disable_fram: file.fram.disable,
        r_fram: file.fram.relaxed,
        triple_intersections: file.fram.triple_intersections,
        boundary_faces: file.clusters.boundary_faces,
        keep_only_largest_cluster: file.clusters.keep_only_largest,
        keep_isolated_fractures: file.clusters.keep_isolated_fractures,
        print_reject_reasons: file.output.print_reject_reasons,
        force_large_fractures: file.generation.force_large_fractures,
        seed: file.seed,
        remove_fractures_less_than: file.generation.remove_fractures_less_than,
        output_rejected_fractures: file.output.rejected_fractures,
        rejects_per_fracture: file.generation.rejects_per_fracture,
        layers,
        layer_volumes,
        regions,
        region_volumes,
        domain_vertices,
        ignore_boundary_faces: file.clusters.ignore_boundary_faces,
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

/// Values of the input file which each family needs.
struct FamilyContext<'a> {
    stop_condition: StopCondition,
    h: f64,
    eps: f64,
    layers: &'a [f64],
    regions: &'a [f64],
    domain_mesh: Option<&'a Rc<DomainMesh>>,
    domain_size: &'a Vector3<f64>,
    domain_size_increase: &'a Vector3<f64>,
}

fn build_family(
    section: &FamilySection,
    id: usize,
    context: &FamilyContext,
) -> Result<FractureFamily, String> {
    if section.layer > context.layers.len() / 2 {
        return Err(format!("layer {} is not defined", section.layer));
    }
    let region_id = match context.domain_mesh {
        Some(mesh) => find_mesh_region_id(section.region, mesh.regions())?,
        None => section.region,
    };
    if region_id > context.regions.len() / 6 {
        return Err(format!("region {} is not defined", section.region));
    }

    let mut builder = FractureFamilyBuilder::default();

    builder
        .id(id)
        .number_of_points(valid_number_of_points(section.number_of_points)?)
        .radius(build_radius_distribution(&section.radius, 3.0 * context.h))
        .aspect_ratio(section.aspect_ratio)
        .orientation(build_orientation(&section.orientation, context.eps)?)
        .layer_id(section.layer)
        .region_id(region_id)
        .layers(context.layers)
        .regions(context.regions)
        .domain_size(context.domain_size)
        .domain_size_increase(context.domain_size_increase)
        .domain_mesh(context.domain_mesh.cloned());

    if let Some(beta) = section.beta {
        builder.beta(beta.to_radians());
    }

    if context.stop_condition == StopCondition::P32 {
        let p32_target = section
            .p32_target
            .ok_or("p32_target is necessary for stop_condition \"p32\"")?;
        builder.p32_target(p32_target);
    }

    builder.build()
}

/// Returns `radius`, or `min_radius` with a warning if `radius` is too small for the mesh.
fn limit_to_min_radius(radius: f64, min_radius: f64, description: &str) -> f64 {
    if radius >= min_radius {
        return radius;
    }

    warn!(
        "{} ({}) is less than 3 times h, reset to {}",
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
fn build_orientation(section: &OrientationSection, eps: f64) -> Result<Fisher, String> {
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
        _ => Err(
            "orientation must have one pair: theta and phi, trend and plunge, or dip and strike"
                .to_string(),
        ),
    }
}

fn valid_number_of_points(number_of_points: usize) -> Result<usize, String> {
    if number_of_points < 3 {
        return Err("number_of_points must be 3 or more".to_string());
    }
    Ok(number_of_points)
}

/// The angles of the section are in degrees. The coordinates of the section
/// are moved so that `domain_center` becomes the origin.
fn build_user_fracture(
    section: &UserFractureSection,
    eps: f64,
    domain_center: &Vector3<f64>,
) -> Result<UserFracture, String> {
    let move_to_origin =
        |point: [f64; 3]| -> [f64; 3] { (Vector3::from(point) - domain_center).into() };

    let has_parameters = section.number_of_points.is_some()
        || section.radius.is_some()
        || section.aspect_ratio.is_some()
        || section.beta.is_some()
        || section.translation.is_some()
        || section.normal.is_some()
        || section.trend_plunge.is_some()
        || section.dip_strike.is_some();

    if let Some(vertices) = &section.vertices {
        if has_parameters {
            return Err("a fracture with vertices can have no other parameters".to_string());
        }
        if vertices.len() < 3 {
            return Err("a fracture must have 3 vertices or more".to_string());
        }
        return Ok(UserFracture::ByVertices(FractureNodesDef {
            vertices: vertices.iter().copied().map(move_to_origin).collect(),
        }));
    }

    let number_of_points = valid_number_of_points(
        section
            .number_of_points
            .ok_or("number_of_points or vertices is necessary")?,
    )?;

    let normal = match (section.normal, section.trend_plunge, section.dip_strike) {
        (Some(normal), None, None) => normal,
        (None, Some([trend, plunge]), None) => {
            normal_from_trend_plunge(trend.to_radians(), plunge.to_radians()).into()
        }
        (None, None, Some([dip, strike])) => {
            normal_from_dip_strike(dip.to_radians(), strike.to_radians()).into()
        }
        _ => {
            return Err("one of normal, trend_plunge, or dip_strike is necessary".to_string());
        }
    };

    Ok(UserFracture::ByParameters(FractureDef {
        eps,
        number_of_points,
        radius: section.radius.ok_or("radius is necessary")?,
        beta: section.beta.ok_or("beta is necessary")?.to_radians(),
        aspect_ratio: section.aspect_ratio.ok_or("aspect_ratio is necessary")?,
        translation: move_to_origin(section.translation.ok_or("translation is necessary")?),
        normal,
    }))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain_mesh::test_meshes::l_shaped_mesh;
    use crate::fracture::fracture_family::PositionSampler;
    use crate::structures::RadiusFunction;

    const GENERAL_OPTIONS: &str = r#"
        seed = 12345
        h = 0.1

        [domain]
        size = [10.0, 10.0, 10.0]
        size_increase = [1.0, 1.0, 1.0]
        truncate_fractures = true

        [generation]
        stop_condition = "n_poly"
        n_poly = 60
        radii_list_increase = 0.1
        rejects_per_fracture = 10
        force_large_fractures = false
        remove_fractures_less_than = 0.0

        [fram]
        disable = false
        relaxed = false
        triple_intersections = true

        [clusters]
        boundary_faces = [true, true, false, false, false, false]
        ignore_boundary_faces = true
        keep_only_largest = false
        keep_isolated_fractures = false

        [output]
        print_reject_reasons = false
        rejected_fractures = false
    "#;

    const RECTANGLE_FAMILY: &str = r#"
        [[family]]
        number_of_points = 4
        probability = 0.6
        aspect_ratio = 1.5
        layer = 0
        region = 0
        radius = { distribution = "log_normal", log_mean = 0.5, sd = 0.4, min = 1.0, max = 4.0 }
        orientation = { trend = 90.0, plunge = 30.0, kappa = 10.0 }
    "#;

    const ELLIPSE_FAMILY: &str = r#"
        [[family]]
        number_of_points = 8
        probability = 0.4
        aspect_ratio = 1.0
        layer = 0
        region = 0
        beta = 90.0
        radius = { distribution = "truncated_power_law", alpha = 2.5, min = 1.0, max = 5.0 }
        orientation = { theta = 0.0, phi = 0.0, kappa = 8.0 }
    "#;

    fn parse(sections: &[&str]) -> Result<(Input, FractureFamilyCollection), DfngenError> {
        parse_input(&sections.concat(), Path::new(""))
    }

    /// Builds the input with the L-shaped mesh as the domain. The mesh has
    /// the regions 1 and 2, and its center is (15, 25, 35).
    fn parse_with_mesh(
        sections: &[&str],
    ) -> Result<(Input, FractureFamilyCollection), DfngenError> {
        build_input(
            toml::from_str(&sections.concat()).unwrap(),
            Some(l_shaped_mesh()),
        )
    }

    fn general_options_without_size() -> String {
        GENERAL_OPTIONS.replace("size = [10.0, 10.0, 10.0]", "")
    }

    fn parse_user_fracture(user_fracture: &str) -> Result<UserFracture, DfngenError> {
        let (mut input, _) = parse(&[GENERAL_OPTIONS, "[[user_fracture]]\n", user_fracture])?;
        Ok(input.user_fractures.remove(0))
    }

    fn error_message(result: Result<(Input, FractureFamilyCollection), DfngenError>) -> String {
        result
            .err()
            .expect("the input must be rejected")
            .to_string()
    }

    #[test]
    fn test_general_options_are_read() {
        let (input, _) = parse(&[GENERAL_OPTIONS]).unwrap();

        assert_eq!(
            (input.seed, input.n_poly, input.domain_size),
            (12345, 60, Vector3::new(10., 10., 10.))
        );
    }

    #[test]
    fn test_input_without_size_and_mesh_is_rejected() {
        let message = error_message(parse(&[&general_options_without_size()]));

        assert!(
            message.contains("domain.size or domain.mesh is necessary"),
            "{message}"
        );
    }

    #[test]
    fn test_input_with_size_and_mesh_is_rejected() {
        let message = error_message(parse_with_mesh(&[GENERAL_OPTIONS]));

        assert!(
            message.contains("domain.size cannot be used with domain.mesh"),
            "{message}"
        );
    }

    #[test]
    fn test_input_with_layers_and_mesh_is_rejected() {
        let options = general_options_without_size()
            .replace("[generation]", "layers = [[-5.0, 0.0]]\n[generation]");
        let message = error_message(parse_with_mesh(&[&options]));

        assert!(
            message.contains("cannot be used with domain.mesh"),
            "{message}"
        );
    }

    #[test]
    fn test_mesh_gives_the_domain_size_and_volume() {
        let (input, _) = parse_with_mesh(&[&general_options_without_size()]).unwrap();

        assert_eq!(
            (input.domain_size, input.domain_volume.round()),
            (Vector3::new(10., 10., 10.), 750.)
        );
    }

    #[test]
    fn test_domain_volume_without_mesh_is_the_box_volume() {
        let (input, _) = parse(&[GENERAL_OPTIONS]).unwrap();

        assert_eq!(input.domain_volume, 1000.);
    }

    #[test]
    fn test_family_region_selects_the_mesh_region_with_that_value() {
        let family = RECTANGLE_FAMILY.replace("region = 0", "region = 2");
        let (_, families) = parse_with_mesh(&[&general_options_without_size(), &family]).unwrap();

        assert_eq!(families.families[0].region_id, 2);
    }

    #[test]
    fn test_family_with_mesh_samples_its_positions_in_its_mesh_region() {
        let family = RECTANGLE_FAMILY.replace("region = 0", "region = 2");
        let (_, families) = parse_with_mesh(&[&general_options_without_size(), &family]).unwrap();

        assert!(matches!(
            families.families[0].position_sampler,
            PositionSampler::Mesh { region_id: 2, .. }
        ));
    }

    #[test]
    fn test_family_without_mesh_samples_its_positions_in_a_box() {
        let (_, families) = parse(&[GENERAL_OPTIONS, RECTANGLE_FAMILY]).unwrap();

        assert!(matches!(
            families.families[0].position_sampler,
            PositionSampler::Box(_)
        ));
    }

    #[test]
    fn test_family_with_region_which_the_mesh_does_not_have_is_rejected() {
        let family = RECTANGLE_FAMILY.replace("region = 0", "region = 3");
        let message = error_message(parse_with_mesh(&[&general_options_without_size(), &family]));

        assert!(
            message.contains("no cell of the domain mesh has the region 3"),
            "{message}"
        );
    }

    #[test]
    fn test_user_fracture_coordinates_are_moved_by_the_mesh_center() {
        let (input, _) = parse_with_mesh(&[
            &general_options_without_size(),
            "[[user_fracture]]\nvertices = [[15.0, 25.0, 35.0], [16.0, 25.0, 35.0], [16.0, 26.0, 35.0]]",
        ])
        .unwrap();

        let UserFracture::ByVertices(fracture) = &input.user_fractures[0] else {
            panic!("the fracture must be defined by vertices");
        };
        assert_eq!(fracture.vertices[0], [0., 0., 0.]);
    }

    #[test]
    fn test_eps_is_derived_from_h() {
        let (input, _) = parse(&[GENERAL_OPTIONS]).unwrap();

        assert_eq!(input.eps, 0.1 * 1e-8);
    }

    #[test]
    fn test_layer_volume_uses_domain_cross_section() {
        let options = GENERAL_OPTIONS.replace(
            "[generation]",
            "layers = [[-5.0, 0.0], [0.0, 2.0]]\n[generation]",
        );

        let (input, _) = parse(&[&options]).unwrap();

        assert_eq!(input.layer_volumes, vec![500., 200.]);
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

        assert_eq!(families.families[0].beta, std::f64::consts::FRAC_PI_2);
    }

    #[test]
    fn test_family_without_beta_has_uniform_rotation() {
        let (_, families) = parse(&[GENERAL_OPTIONS, RECTANGLE_FAMILY]).unwrap();

        assert!(!families.families[0].beta_distribution);
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
    fn test_truncate_fractures_true_gives_the_truncate_mode() {
        let (input, _) = parse(&[GENERAL_OPTIONS]).unwrap();

        assert_eq!(input.domain_boundary_mode, DomainBoundaryMode::Truncate);
    }

    #[test]
    fn test_truncate_fractures_false_gives_the_keep_whole_mode() {
        let options =
            GENERAL_OPTIONS.replace("truncate_fractures = true", "truncate_fractures = false");

        let (input, _) = parse(&[&options]).unwrap();

        assert_eq!(input.domain_boundary_mode, DomainBoundaryMode::KeepWhole);
    }

    #[test]
    fn test_input_without_truncate_fractures_is_rejected() {
        let options = GENERAL_OPTIONS.replace("truncate_fractures = true", "");

        let message = error_message(parse(&[&options]));

        assert!(
            message.contains("missing field `truncate_fractures`"),
            "{}",
            message
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

        let message = error_message(parse(&[GENERAL_OPTIONS, &family]));

        assert!(message.contains("family 1: orientation"), "{}", message);
    }

    #[test]
    fn test_p32_stop_condition_needs_family_target() {
        let options = GENERAL_OPTIONS.replace("\"n_poly\"", "\"p32\"");

        let message = error_message(parse(&[&options, ELLIPSE_FAMILY]));

        assert!(message.contains("p32_target is necessary"), "{}", message);
    }

    #[test]
    fn test_family_in_undefined_layer_is_rejected() {
        let family = ELLIPSE_FAMILY.replace("layer = 0", "layer = 1");

        let message = error_message(parse(&[GENERAL_OPTIONS, &family]));

        assert!(message.contains("layer 1 is not defined"), "{}", message);
    }

    #[test]
    fn test_user_fracture_by_normal_is_read() {
        let fracture = parse_user_fracture(
            r#"number_of_points = 4
            radius = 4.0
            aspect_ratio = 1.0
            beta = 0.0
            translation = [0.0, 0.0, 0.0]
            normal = [0.0, 0.0, 1.0]"#,
        )
        .unwrap();

        assert_eq!(fracture.to_poly().normal, Vector3::new(0., 0., 1.));
    }

    #[test]
    fn test_user_ellipse_by_trend_plunge_is_in_degrees() {
        let fracture = parse_user_fracture(
            r#"number_of_points = 8
            radius = 2.0
            aspect_ratio = 1.0
            beta = 0.0
            translation = [0.0, 0.0, 0.0]
            trend_plunge = [90.0, 30.0]"#,
        )
        .unwrap();

        let expected = normal_from_trend_plunge(90f64.to_radians(), 30f64.to_radians());
        assert_eq!(fracture.to_poly().normal, expected);
    }

    #[test]
    fn test_angle_unit_is_rejected() {
        let options = GENERAL_OPTIONS.replace("h = 0.1", "h = 0.1\nangle_unit = \"radian\"");

        let message = error_message(parse(&[&options]));

        assert!(
            message.contains("unknown field `angle_unit`"),
            "{}",
            message
        );
    }

    #[test]
    fn test_user_fracture_by_vertices_is_read() {
        let fracture = parse_user_fracture(
            r#"vertices = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]]"#,
        )
        .unwrap();

        assert_eq!(fracture.to_poly().vertices.len(), 3);
    }

    #[test]
    fn test_user_fracture_with_two_vertices_is_rejected() {
        let result = parse_user_fracture(r#"vertices = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]"#);

        let message = result.expect_err("the input must be rejected").to_string();
        assert!(message.contains("user_fracture 1: a fracture must have 3 vertices or more"));
    }

    #[test]
    fn test_user_fracture_with_vertices_and_parameters_is_rejected() {
        let result = parse_user_fracture(
            r#"number_of_points = 4
            vertices = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]]"#,
        );

        let message = result.expect_err("the input must be rejected").to_string();
        assert!(message.contains("a fracture with vertices can have no other parameters"));
    }

    #[test]
    fn test_family_with_two_points_is_rejected() {
        let family = ELLIPSE_FAMILY.replace("number_of_points = 8", "number_of_points = 2");

        let message = error_message(parse(&[GENERAL_OPTIONS, &family]));

        assert!(
            message.contains("family 1: number_of_points must be 3 or more"),
            "{}",
            message
        );
    }

    #[test]
    fn test_user_fracture_without_orientation_is_rejected() {
        let result = parse_user_fracture(
            r#"number_of_points = 4
            radius = 4.0
            aspect_ratio = 1.0
            beta = 0.0
            translation = [0.0, 0.0, 0.0]"#,
        );

        let message = result.expect_err("the input must be rejected").to_string();
        assert!(message.contains("one of normal, trend_plunge, or dip_strike"));
    }
}
