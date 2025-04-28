use parry3d_f64::na::{Point3, Vector3};
use tracing::{error, info};

use super::read_input_functions::read_domain_vertices;

use crate::fracture::fracture_family::FractureFamilyCollection;
use crate::io::read_input_functions::InputReader;

#[derive(Default, Debug)]
pub struct ExternalFractureFiles {
    /// File name of user polygons defined by coordinates
    pub user_poly_by_coord_file: Option<String>,

    /// File name of user ellipses defined by coordinates
    pub user_ell_by_coord_file: Option<String>,

    /// File name of user retangles defined by coordinates
    pub user_rect_by_coord_file: Option<String>,

    /// File name of user ellipses
    pub user_ell_file: Option<String>,

    /// File name of user rectangles
    pub user_rect_file: Option<String>,
}

#[allow(non_snake_case)]
#[derive(Default, Debug)]
pub struct Input {
    /// DFN generation stop condition. 0 - nPoly option, 1 - P32 option.
    pub stopCondition: u8,

    /// Number of polygons to place in the DFN when uisng nPoly stopCondition option.
    pub nPoly: usize,

    /// Domain size with dimension x*y*z for DFN, centered at the origin.
    pub domainSize: Vector3<f64>,

    /// Minimum feature size, FRAM parameter.
    pub h: f64,

    /// epsilon
    pub eps: f64,

    /// Percent to increase the size of the pre-generated radii lists, per family.
    /// Example: 0.2 will increase the size of the list by %20. See example input files
    /// for more details.
    pub radiiListIncrease: f64,

    /// This option disables the FRAM algorithm. There will be no
    /// fracture rejections or fine mesh. Defaults visualizationMode to 1
    pub disableFram: bool,

    ///  Used during meshing:
    ///      0 - Creates a fine mesh, according to h parameter;
    ///      1 - Produce only first round of triangulations. In this case no
    ///          modeling of flow and transport is possible.
    pub visualizationMode: bool,

    /// This option uses a relaxed version of the FRAM algorithm. The mesh may not
    /// be perfectly conforming
    pub rFram: bool,

    /// Accept or reject triple intersections
    ///     False - Off (Reject)
    ///     True  - On  (Accept)
    pub tripleIntersections: bool,

    /// DFN will only keep clusters with connections to
    /// domain boundaries which are set to 1:
    ///
    /// boundaryFaces[0] = +X domain boundary
    /// boundaryFaces[1] = -X domain boundary
    /// boundaryFaces[2] = +Y domain boundary
    /// boundaryFaces[3] = -Y domain boundary
    /// boundaryFaces[4] = +Z domain boundary
    /// boundaryFaces[5] = -Z domain boundary
    pub boundaryFaces: [bool; 6],

    /// 0 - Keep any clusters which connects the specified
    ///     boundary faces in boundaryFaces option below
    /// 1 - Keep only the largest cluster which connects
    ///     the specified boundary faces in boundaryFaces option below.
    ///
    /// If ignoreBoundaryFaces is also set to 1, DFNGen will keep the largest
    /// cluster which connects at least any two sides of the domain.
    pub keepOnlyLargestCluster: bool,

    /// 0 - remove isolated fractures and clusters
    /// 1 - Keep isolated fractures and clusters
    pub keepIsolatedFractures: bool,

    /// Useful for debugging,
    /// This option will print all fracture rejection reasons as they occur.
    ///     0 - Disable
    ///     1 - Print all rejection reasons to screen
    pub printRejectReasons: bool,

    /// Outputs radii files after isolated fracture removal.
    /// One file per family.
    ///     0: Do not create output files of radii per family
    ///     1: Creates output files per family, containing a list
    ///        of the family's fracture radii that is in the final DFN
    pub outputFinalRadiiPerFamily: bool,

    /// Outputs radii files before isolated fracture removal.
    /// One file per family.
    ///     0: Do not create output files of radii per family
    ///     1: Creates output files per family, containing a list
    ///        of the family's fracture radii in the domain before isolated
    ///        fracture removal.
    pub outputAcceptedRadiiPerFamily: bool,

    /// False - User ellipses will be inserted first
    /// True  - User rectangles will be inserted first
    pub insertUserRectanglesFirst: bool,

    /// Inserts the largest possible fracture for each defined fracture family,
    /// defined by the user-defined maxium radius
    ///     0 - Off (Do not force insertion of larest fractures)
    ///     1 - On  (Force insertion of largest fractures)
    pub forceLargeFractures: bool,

    /// Seed for random generator.
    pub seed: u64,

    /// Size increase for inserting fracture centers outside the domain.
    /// Fracture will be truncated based on domainSize above.
    /// Increases the entire width by this ammount. So, {1,1,1} will increase
    /// the domain by adding .5 to the +x, and subbtracting .5 to the -x, etc
    pub domainSizeIncrease: Vector3<f64>,

    /// Selection of orientation Option
    /// 0 - spherical coordinates
    /// 1 - trend / plunge
    /// 2 - dip / strike
    pub orientationOption: u8,

    /// Number of rectangular families defined below.
    /// Having this option = 0 will ignore all rectangular family variables.
    pub nFamRect: usize,

    /// Number of ellipse families defined below.
    /// Having this option = 0 will ignore all ellipse family variables.
    pub nFamEll: usize,

    ///    0 - Ignore this option, keep all fractures.
    ///
    /// (>0) - Size of minimum fracture radius. Fractures smaller than
    ///        defined radius will be removed AFTER DFN generation.
    ///
    ///        Minimum and maximum size options under fracture family
    ///        distributions will still be used while generating the DFN.
    pub removeFracturesLessThan: f64,

    /// Caution: Can create very large files.
    /// Outputs all fractures which were generated during
    /// DFN generation (Accepted + Rejected).
    ///     False: Do not output all radii file.
    ///     True:  Include file of all raddii, acepted + rejected fractures,
    ///            in output files (radii_All.dat).
    pub outputAllRadii: bool,

    /// If a fracture is rejected, it will be re-translated
    /// to a new position this number of times.
    ///
    /// This helps hit distribution targets for stochastic families
    /// families (Set to 1 to ignore this feature)
    pub rejectsPerFracture: usize,

    /// Array of layers:
    /// e.g. {+z1, -z1, +z2, -z2, ... , +zn, -zn}
    pub layers: Vec<f64>,

    /// Array of volumes for each defined layer, in the same order
    /// which layers were listed.
    pub layerVol: Vec<f64>,

    /// Array of regions:
    /// e.g. {+z1, -z1, +z2, -z2, ... , +zn, -zn}
    pub regions: Vec<f64>,

    /// Array of volumes for each defined layer, in the same order
    /// which regions were listed.
    pub regionVol: Vec<f64>,

    /// flag if the domain is pruned down to a final domain size
    /// bool polygonBoundaryFlag = false;
    pub polygonBoundaryFlag: bool,

    /// Number of points on the 2D boundary of the polygon domain
    pub numOfDomainVertices: usize,

    /// Vector of points defining the 2D boundary of the domain polygon
    pub domainVertices: Vec<Point3<f64>>,

    /// False - Use boundaryFaces option.
    /// True  - Ignore boundaryFaces option, keep all clusters
    ///         and remove fractures with no intersections
    pub ignoreBoundaryFaces: bool,

    pub ext_fracture_files: ExternalFractureFiles,
}

/// Reads in all input variables.
/// Creates Shape structure array from user input if
/// using stochastic fracture families.
///
/// # Arguments
///
/// * `input` - Path to input file
pub fn read_input(input_file: &str) -> (Input, FractureFamilyCollection) {
    let mut input_var = Input::default();
    let mut fracture_families = Vec::new();

    info!("DFN Generator Input File: {}", input_file);
    let mut input_reader = InputReader::new(input_file);

    macro_rules! input_var {
        ($var_name:ident) => {
            input_reader.read_value(
                &format!("{}:", stringify!($var_name)),
                &mut input_var.$var_name,
            );
        };
    }

    input_var!(stopCondition);
    input_var!(printRejectReasons);
    input_var!(domainSize);

    let mut n_layers: usize = 0;
    input_reader.read_value("numOfLayers:", &mut n_layers);
    if n_layers > 0 {
        input_var.layers.reserve(n_layers * 2); // Multiply by 2 for +z and -z for each layer
        input_var!(layers);

        for i in 0..n_layers {
            let idx = i * 2;
            info!(
                "    Layer {}{{-z,+z}}: {:?}, Volume: ",
                i + 1,
                &input_var.layers[idx..idx + 2]
            );
            let vol = input_var.domainSize[0]
                * input_var.domainSize[1]
                * ((input_var.layers[idx + 1] as isize - input_var.layers[idx] as isize).abs()
                    as f64);
            info!("{} m^3", vol);
        }
    }

    let mut n_regions: usize = 0;
    input_reader.read_value("numOfRegions:", &mut n_regions);
    if n_regions > 0 {
        input_var.regions.reserve(n_regions * 6); // Multiply by 6 xmin, xmax, ymin, ymax, zmin, zmax
        input_var!(regions);

        info!("Number of Regions: {}", n_regions);

        for i in 0..n_regions {
            let idx = i * 6;
            info!(
                " Region {}: {{-x,+x,-y,+y,-z,+z}}: {:?}",
                i + 1,
                &input_var.regions[idx..idx + 6]
            );
            let vol = (input_var.regions[idx + 1] - input_var.regions[idx]).abs()
                * (input_var.regions[idx + 3] - input_var.regions[idx + 2]).abs()
                * (input_var.regions[idx + 5] - input_var.regions[idx + 4]).abs();
            info!(" Volume: {} m^3", vol);
        }
    }

    input_var!(h);
    info!("h: {}", input_var.h);

    // Set epsilon
    input_var.eps = input_var.h * 1e-8;

    input_var!(disableFram);

    if input_var.disableFram {
        info!("FRAM IS DISABLED");
    }

    input_var!(rFram);

    if input_var.rFram {
        info!("Running with relaxed FRAM. Mesh may not be fully conforming");
    }

    input_var!(tripleIntersections);
    input_var!(forceLargeFractures);
    input_var!(visualizationMode);

    if input_var.disableFram {
        input_var.visualizationMode = true;
    }

    input_var!(outputAllRadii);
    input_var!(outputFinalRadiiPerFamily);
    input_var!(outputAcceptedRadiiPerFamily);
    input_var!(seed);
    input_var!(domainSizeIncrease);
    input_var!(keepOnlyLargestCluster);
    input_var!(keepIsolatedFractures);
    input_var!(ignoreBoundaryFaces);
    input_var!(boundaryFaces);
    input_var!(rejectsPerFracture);
    input_var!(nFamRect);
    input_var!(nFamEll);
    input_var!(removeFracturesLessThan);
    input_var!(orientationOption);

    if input_var.orientationOption == 0 {
        info!("Expecting Theta and phi for orientations");
    } else if input_var.orientationOption == 1 {
        info!("Expecting Trend and Plunge for orientations");
    } else if input_var.orientationOption == 2 {
        info!("Expecting Dip and Strike (RHR) for orientations");
    }

    input_var!(polygonBoundaryFlag);

    if input_var.polygonBoundaryFlag {
        info!("Expecting Polygon Boundary for domain edges");
        let mut poly_boundary_file: String = String::new();
        input_reader.read_value("polygonBoundaryFile:", &mut poly_boundary_file);
        info!("Polygon Boundary File: {}", &poly_boundary_file);
        read_domain_vertices(&mut input_var, &poly_boundary_file);
        info!(
            "There are {} Vertices on the boundary",
            input_var.numOfDomainVertices
        );

        for i in 0..input_var.numOfDomainVertices {
            info!(
                "Vertex {}: {{{}, {}}}",
                i + 1,
                input_var.domainVertices[i].x,
                input_var.domainVertices[i].y
            );
        }
    }

    fracture_families.extend(input_reader.read_fracture_family(
        "e",
        &input_var.layers,
        &input_var.regions,
        &input_var.domainSize,
        &input_var.domainSizeIncrease,
    ));
    fracture_families.extend(input_reader.read_fracture_family(
        "r",
        &input_var.layers,
        &input_var.regions,
        &input_var.domainSize,
        &input_var.domainSizeIncrease,
    ));

    let fracture_family_probabilities = if !fracture_families.is_empty() {
        input_var!(radiiListIncrease);

        if input_var.stopCondition == 0 {
            // npoly option
            input_var!(nPoly);
        } else if input_var.stopCondition == 1 {
        } else {
            error!("Invalid stopCondition option. Must be 0: nPoly or 1: p32");
            panic!()
        }

        let mut fam_prob = Vec::with_capacity(fracture_families.len());
        input_reader.read_value("famProb:", &mut fam_prob);
        fam_prob
    } else {
        Vec::new()
    };

    //
    // Get external fracture definition files
    //

    let mut user_ell = false;
    input_reader.read_value("userEllipsesOnOff:", &mut user_ell);
    if user_ell {
        let mut path: String = String::new();
        input_reader.read_value("UserEll_Input_File_Path:", &mut path);
        input_var.ext_fracture_files.user_ell_file = Some(path);
    }

    let mut user_rect = false;
    input_reader.read_value("userRectanglesOnOff:", &mut user_rect);
    if user_rect {
        let mut path: String = String::new();
        input_reader.read_value("UserRect_Input_File_Path:", &mut path);
        input_var.ext_fracture_files.user_rect_file = Some(path);
    }

    let mut user_polygon_by_coord = false;
    input_reader.read_value("userPolygonByCoord:", &mut user_polygon_by_coord);
    if user_polygon_by_coord {
        let mut path = String::new();
        input_reader.read_value("PolygonByCoord_Input_File_Path:", &mut path);
        input_var.ext_fracture_files.user_poly_by_coord_file = Some(path);
    }

    let mut user_ell_by_coord = false;
    input_reader.read_value("userEllByCoord:", &mut user_ell_by_coord);
    if user_ell_by_coord {
        let mut path = String::new();
        input_reader.read_value("EllByCoord_Input_File_Path:", &mut path);
        input_var.ext_fracture_files.user_ell_by_coord_file = Some(path);
    }

    let mut user_rect_by_coord = false;
    input_reader.read_value("userRecByCoord:", &mut user_rect_by_coord);
    if user_rect_by_coord {
        let mut path = String::new();
        input_reader.read_value("RectByCoord_Input_File_Path:", &mut path);
        input_var.ext_fracture_files.user_rect_by_coord_file = Some(path);
    }

    if (user_rect || user_rect_by_coord) && (user_ell || user_ell_by_coord) {
        input_var!(insertUserRectanglesFirst);
    } else {
        input_var.insertUserRectanglesFirst = false;
    }

    (
        input_var,
        FractureFamilyCollection {
            families: fracture_families,
            probabilities: fracture_family_probabilities.clone(),
            original_probabilities: fracture_family_probabilities,
        },
    )
}
