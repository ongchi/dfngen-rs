use std::path::Path;
use std::rc::Rc;

use parry3d_f64::na::Vector3;
use serde::Deserialize;
use tracing::{debug, info, warn};

use crate::distribution::{Fisher, LogNormalSampling};
use crate::domain_mesh::{DomainMesh, MeshRegion};
use crate::error::DfngenError;
use crate::fracture::fracture_family::{
    FractureFamily, FractureFamilyBuilder, FractureFamilyCollection,
};
use crate::fracture::poly::DomainBoundaryMode;
use crate::fracture::user_defined_fractures::UserFracture;
use crate::io::vtp::{read_user_fractures, FracturePolygon};
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

#[derive(Debug)]
pub struct Input {
    /// DFN generation stop condition.
    pub stop_condition: StopCondition,

    /// Number of polygons to place in the DFN with the `NPoly` stop condition.
    pub n_poly: usize,

    /// Size of the bounding box of the domain mesh. The generation has the
    /// center of this box at the origin.
    pub domain_size: Vector3<f64>,

    /// Volume of the domain: the sum of the cell volumes of the domain mesh.
    pub domain_volume: f64,

    /// The cells which make the domain.
    pub domain_mesh: Rc<DomainMesh>,

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

    /// Volume of each region of the domain mesh, in the order of
    /// `DomainMesh::regions`.
    pub region_volumes: Vec<f64>,

    /// False - Use boundaryFaces option.
    /// True  - Ignore boundaryFaces option, keep all clusters
    ///         and remove fractures with no intersections
    pub ignore_boundary_faces: bool,

    /// User-defined fractures, in the order of their file. The vertices have
    /// the coordinates of the generation.
    pub user_fractures: Vec<UserFracture>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct InputFile {
    seed: u64,
    h: f64,
    /// Path of an ASCII `.vtp` file with the user-defined fractures (optional).
    /// A relative path starts at the folder of the input file.
    user_fractures: Option<String>,
    domain: DomainSection,
    generation: GenerationSection,
    fram: FramSection,
    clusters: ClustersSection,
    output: OutputSection,
    #[serde(default, rename = "family")]
    families: Vec<FamilySection>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct DomainSection {
    /// Path of an ASCII `.vtu` file which gives the domain.
    /// A relative path starts at the folder of the input file.
    mesh: String,
    /// True: cut each fracture at the domain boundary. False: keep it whole.
    truncate_fractures: bool,
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

    let file = parse_input_file(&text)?;
    let input_folder = Path::new(input_file).parent().unwrap_or(Path::new(""));
    let domain_mesh = read_domain_mesh(&input_folder.join(&file.domain.mesh))?;
    let user_fractures = match &file.user_fractures {
        Some(fracture_file) => read_user_fractures(&input_folder.join(fracture_file))?,
        None => Vec::new(),
    };

    build_input(file, domain_mesh, user_fractures)
}

fn parse_input_file(text: &str) -> Result<InputFile, DfngenError> {
    toml::from_str(text).map_err(|e| invalid_input(e.to_string()))
}

/// Position in the mesh regions (starts at 1) of the region value of a family
/// or of a user-defined fracture. 0 is the whole mesh.
fn find_mesh_region_id(region_value: u64, mesh_regions: &[MeshRegion]) -> Result<usize, String> {
    if region_value == 0 {
        return Ok(0);
    }

    mesh_regions
        .iter()
        .position(|region| region.value == region_value)
        .map(|position| position + 1)
        .ok_or(format!(
            "no cell of the domain mesh has the region {}",
            region_value
        ))
}

/// Builds the input options and the stochastic fracture families.
///
/// # Arguments
///
/// * `file` - Content of the TOML input file
/// * `domain_mesh` - Domain mesh of `domain.mesh`
/// * `fracture_polygons` - User-defined fractures, as their file gives them
fn build_input(
    file: InputFile,
    domain_mesh: DomainMesh,
    fracture_polygons: Vec<FracturePolygon>,
) -> Result<(Input, FractureFamilyCollection), DfngenError> {
    let eps = file.h * 1e-8;
    let domain_mesh = Rc::new(domain_mesh);

    if file.fram.disable {
        info!("FRAM IS DISABLED");
    }

    if file.fram.relaxed {
        info!("Running with relaxed FRAM. Mesh may not be fully conforming");
    }

    let region_volumes: Vec<f64> = domain_mesh
        .regions()
        .iter()
        .map(|region| region.volume)
        .collect();

    debug!(
        "Input: h = {}, region volumes {:?} m^3",
        file.h, region_volumes
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

    let user_fractures = fracture_polygons
        .iter()
        .enumerate()
        .map(|(i, polygon)| {
            build_user_fracture(polygon, &domain_mesh).map_err(|message| {
                invalid_input(format!("user-defined fracture {}: {}", i + 1, message))
            })
        })
        .collect::<Result<Vec<_>, _>>()?;

    let family_context = FamilyContext {
        stop_condition,
        h: file.h,
        eps,
        domain_mesh: &domain_mesh,
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
        domain_size: domain_mesh.size(),
        domain_volume: domain_mesh.volume(),
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
        region_volumes,
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

/// Moves a polygon of the user-defined fracture file to the coordinates of
/// the generation (the center of the mesh at the origin), and finds its region
/// in the mesh.
fn build_user_fracture(
    polygon: &FracturePolygon,
    domain_mesh: &DomainMesh,
) -> Result<UserFracture, String> {
    let domain_center = domain_mesh.center();

    Ok(UserFracture {
        vertices: polygon
            .vertices
            .iter()
            .map(|&vertex| (Vector3::from(vertex) - domain_center).into())
            .collect(),
        region_id: find_mesh_region_id(polygon.region, domain_mesh.regions())?,
    })
}

/// Values of the input file which each family needs.
struct FamilyContext<'a> {
    stop_condition: StopCondition,
    h: f64,
    eps: f64,
    domain_mesh: &'a Rc<DomainMesh>,
}

fn build_family(
    section: &FamilySection,
    id: usize,
    context: &FamilyContext,
) -> Result<FractureFamily, String> {
    let region_id = find_mesh_region_id(section.region as u64, context.domain_mesh.regions())?;

    let mut builder = FractureFamilyBuilder::default();

    builder
        .id(id)
        .number_of_points(valid_number_of_points(section.number_of_points)?)
        .radius(build_radius_distribution(&section.radius, 3.0 * context.h))
        .aspect_ratio(section.aspect_ratio)
        .orientation(build_orientation(&section.orientation, context.eps)?)
        .region_id(region_id)
        .domain_mesh(context.domain_mesh.clone());

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

#[cfg(test)]
pub mod test_inputs {
    use super::*;
    use crate::domain_mesh::test_meshes::l_shaped_mesh;

    pub const GENERAL_OPTIONS: &str = r#"
        seed = 12345
        h = 0.1

        [domain]
        mesh = "domain.vtu"
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

    /// Reads the sections with the L-shaped mesh as the domain. The mesh has
    /// the regions 1 and 2, and its center is (15, 25, 35).
    pub fn parse(sections: &[&str]) -> Result<(Input, FractureFamilyCollection), DfngenError> {
        build_input(
            parse_input_file(&sections.concat())?,
            l_shaped_mesh(),
            Vec::new(),
        )
    }

    /// The input of `GENERAL_OPTIONS`
    pub fn general_input() -> Input {
        parse(&[GENERAL_OPTIONS]).unwrap().0
    }
}

#[cfg(test)]
mod tests {
    use super::test_inputs::{parse, GENERAL_OPTIONS};
    use super::*;
    use crate::domain_mesh::test_meshes::l_shaped_mesh;
    use crate::structures::RadiusFunction;

    const RECTANGLE_FAMILY: &str = r#"
        [[family]]
        number_of_points = 4
        probability = 0.6
        aspect_ratio = 1.5
        region = 0
        radius = { distribution = "log_normal", log_mean = 0.5, sd = 0.4, min = 1.0, max = 4.0 }
        orientation = { trend = 90.0, plunge = 30.0, kappa = 10.0 }
    "#;

    const ELLIPSE_FAMILY: &str = r#"
        [[family]]
        number_of_points = 8
        probability = 0.4
        aspect_ratio = 1.0
        region = 0
        beta = 90.0
        radius = { distribution = "truncated_power_law", alpha = 2.5, min = 1.0, max = 5.0 }
        orientation = { theta = 0.0, phi = 0.0, kappa = 8.0 }
    "#;

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
    fn test_input_without_mesh_is_rejected() {
        let options = GENERAL_OPTIONS.replace("mesh = \"domain.vtu\"", "");

        let message = error_message(parse(&[&options]));

        assert!(message.contains("missing field `mesh`"), "{message}");
    }

    #[test]
    fn test_input_with_domain_size_is_rejected() {
        let options = GENERAL_OPTIONS.replace(
            "truncate_fractures = true",
            "truncate_fractures = true\nsize = [10.0, 10.0, 10.0]",
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

    #[test]
    fn test_family_region_selects_the_mesh_region_with_that_value() {
        let family = RECTANGLE_FAMILY.replace("region = 0", "region = 2");
        let (_, families) = parse(&[GENERAL_OPTIONS, &family]).unwrap();

        assert_eq!(families.families[0].region_id, 2);
    }

    #[test]
    fn test_family_with_region_which_the_mesh_does_not_have_is_rejected() {
        let family = RECTANGLE_FAMILY.replace("region = 0", "region = 3");
        let message = error_message(parse(&[GENERAL_OPTIONS, &family]));

        assert!(
            message.contains("no cell of the domain mesh has the region 3"),
            "{message}"
        );
    }

    /// Builds the input with one triangle in the L-shaped mesh as the
    /// user-defined fracture. The first vertex is the center of the mesh.
    fn parse_with_triangle_of_region(
        region: u64,
    ) -> Result<(Input, FractureFamilyCollection), DfngenError> {
        let triangle = FracturePolygon {
            vertices: vec![[15., 25., 35.], [16., 25., 35.], [16., 26., 35.]],
            region,
        };

        build_input(
            parse_input_file(GENERAL_OPTIONS).unwrap(),
            l_shaped_mesh(),
            vec![triangle],
        )
    }

    #[test]
    fn test_user_fracture_coordinates_are_moved_by_the_mesh_center() {
        let (input, _) = parse_with_triangle_of_region(0).unwrap();

        assert_eq!(input.user_fractures[0].vertices[0], [0., 0., 0.]);
    }

    #[test]
    fn test_user_fracture_region_selects_the_mesh_region_with_that_value() {
        let (input, _) = parse_with_triangle_of_region(2).unwrap();

        assert_eq!(input.user_fractures[0].region_id, 2);
    }

    #[test]
    fn test_user_fracture_with_region_which_the_mesh_does_not_have_is_rejected() {
        let message = error_message(parse_with_triangle_of_region(3));

        assert!(
            message
                .contains("user-defined fracture 1: no cell of the domain mesh has the region 3"),
            "{message}"
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
    fn test_family_with_two_points_is_rejected() {
        let family = ELLIPSE_FAMILY.replace("number_of_points = 8", "number_of_points = 2");

        let message = error_message(parse(&[GENERAL_OPTIONS, &family]));

        assert!(
            message.contains("family 1: number_of_points must be 3 or more"),
            "{}",
            message
        );
    }
}
