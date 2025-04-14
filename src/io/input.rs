use std::fs::File;

use itertools::zip_eq;
use parry3d_f64::na::{Point3, Vector3};

use super::read_input_functions::{
    get_cords, get_rect_coords, read_domain_vertices, search_var, ReadFromTextFile,
};

use crate::{
    distribution::{generating_points::generate_theta, Fisher},
    structures::Shape,
};

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

    // Only output select files for ECPM upscaling
    //     0: Output all files
    //     1: Only output files required for ECPM upscaling.
    //     polygon.dat, radii_final.dat
    // ecpmOutput: bool,
    // ecpmOutput:bool = false,
    /// Beta is the rotation around the polygon's normal vector
    ///     0 - Uniform distribution [0, 2PI)
    ///     1 - Constant angle (specefied below by 'ebeta')
    ebetaDistribution: Vec<bool>,

    /// Beta is the rotation around the polygon's normal vector
    ///     0: Uniform distribution [0, 2PI)
    ///     1: Constant angle (specefied below by 'rbeta')
    rbetaDistribution: Vec<bool>,

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

    /// Each element is the probability of chosing a fracture from
    /// the element's corresponding family to be inserted into the DFN.
    ///
    /// The famProb elements should add up to 1.0 (for %100).
    /// The probabilities are listed in order of families starting with all
    /// stochastic ellipses, and then all stochastic rectangles.
    ///
    /// For example:
    /// If  then there are two ellipse families, each with probabiliy .3,
    /// and two rectangle families, each with probabiliy .2, famProb will be:
    /// famProb: {.3,.3,.2,.2}, famProb elements must add to 1
    pub famProb: Vec<f64>,

    /// Holds a copy of famProb. famProb elements can change as different families
    /// hit their P32 requirement when using the P32 stopCondition option.
    pub famProbOriginal: Vec<f64>,

    /// Mandatory parameter if using statistically generated ellipses.
    /// Statistical distribution options:
    ///
    /// Holds number of elements equal to the number of shape families.
    ///
    ///     1 - Log-normal distribution
    ///     2 - Truncated power law distribution
    ///     3 - Exponential distribution
    ///     4 - Constant
    edistr: Vec<u8>,

    /// Aspect ratio array for stochastic ellipses.
    easpect: Vec<f64>,

    /// Number of vertices used in creating each elliptical
    /// fracture family. Number of elements must match number
    /// of ellipse families
    ///
    /// Holds number of elements equal to the number of ellipse families.
    enumPoints: Vec<usize>,

    /// Rotation around the fractures' normal vector.
    /// Ellipse family parameter.
    ebeta: Vec<f64>,

    /// Parameter for the fisher distribnShaprutions. The
    /// bigger, the more similar (less diverging) are the
    /// elliptical familiy's normal vectors.
    ekappa: Vec<f64>,

    /// Log-normal ellipse parameter. Mean of the underlying normal distribution.
    eLogMean: Vec<f64>,

    /// Log-normal ellipse parameter. Standard deviation of the underlying normal distribution
    esd: Vec<f64>,

    /// Exponential ellipse parameter. Mean values for exponential distributions, defined per family.
    eExpMean: Vec<f64>,

    /// Log-normal rectangle parameter. Minimum radius.
    rLogMin: Vec<f64>,

    /// Log-normal rectangle parameter. Maximum radius.
    rLogMax: Vec<f64>,

    /// Exponential rectangle parameter. Minimum radius.
    rExpMin: Vec<f64>,

    /// Exponential rectangle parameter. Maximum radius.
    rExpMax: Vec<f64>,

    /// Log-normal ellipse parameter. Minimum radius.
    eLogMin: Vec<f64>,

    /// Log-normal ellipse parameter. Maximum radius.
    eLogMax: Vec<f64>,

    /// Exponential ellipse parameter. Minimum radius.
    eExpMin: Vec<f64>,

    /// Exponential ellipse parameter. Maximum radius.
    eExpMax: Vec<f64>,

    /// Contant ellipse parameter. Constant radius.
    econst: Vec<f64>,

    /// Truncated power-law ellipse parameter. Minimum radius.
    emin: Vec<f64>,

    /// Truncated power-law ellipse parameter. Maximum radius.
    emax: Vec<f64>,

    /// Truncated power-law ellipse distribution parameter.
    ealpha: Vec<f64>,

    /// Elliptical families target fracture intensities per family
    /// when using stopCondition = 1, P32 option.
    e_p32Targets: Vec<f64>,

    /// Mandatory parameter if using statistically generated rectangles.
    ///
    /// Holds number of elements equal to the number of shape families.\
    ///
    /// Rectangle statistical distribution options:
    ///     1 - log-normal distribution
    ///     2 - truncated power law distribution
    ///     3 - exponential distribution
    ///     4 - constant
    rdistr: Vec<u8>,

    /// Aspect ratio for stochasic rectangles.
    raspect: Vec<f64>,

    ///    0 - Ignore this option, keep all fractures.
    ///
    /// (>0) - Size of minimum fracture radius. Fractures smaller than
    ///        defined radius will be removed AFTER DFN generation.
    ///
    ///        Minimum and maximum size options under fracture family
    ///        distributions will still be used while generating the DFN.
    pub removeFracturesLessThan: f64,

    /// Rotation around the normal vector.
    rbeta: Vec<f64>,

    /// Parameter for the fisher distribnShaprutions. The
    /// bigger, the more similar (less diverging) are the
    /// rectangle family's normal vectors.
    rkappa: Vec<f64>,

    /// Log-normal rectangle parameter. Standard deviation of the underlying normal distribution
    rLogMean: Vec<f64>,

    /// Log-normal rectangle parameter. Standard deviation of the underlying normal distribution
    rsd: Vec<f64>,

    /// Truncated power-law rectangle parameter. Minimum radius.
    rmin: Vec<f64>,

    /// Truncated power-law rectangle parameter. Maximum radius.
    rmax: Vec<f64>,

    /// Truncated power-law rectangle distribution parameter.
    ralpha: Vec<f64>,

    /// Rectangular families target fracture intensities per family
    /// when using stopCondition = 1, P32 option.
    r_p32Targets: Vec<f64>,

    /// Exponential rectangle parameter. Maximum radius.
    rExpMean: Vec<f64>,

    /// Constant rectangle parameter. Constant radius.
    rconst: Vec<f64>,

    /// True  - The user is using user defined ellipses.
    /// False - No user defined ellipses are being used.
    pub userEllipsesOnOff: bool,

    /// Number of defined, user defined ellipses.
    pub nUserEll: usize,

    /// All angles from input file for stochastic ellipses are in:
    ///     True  - Degrees
    ///     False - Radians
    pub ueAngleOption: bool,

    /// User ellipses radii array.
    pub ueRadii: Vec<f64>,

    /// User ellipses beta array.
    pub ueBeta: Vec<f64>,

    /// User ellipses aspect ratio array.
    pub ueaspect: Vec<f64>,

    /// User ellipses translation array.
    pub uetranslation: Vec<f64>,

    /// User Orientation Option for ellipses
    /// 0 = normal vector
    /// 1 = trend / plunge
    /// 2 = dip / strike
    userEllOrientationOption: u8,

    /// User ellipses normal vector array.
    pub uenormal: Vec<f64>,

    /// User ellipses trend and plunge array.
    ueTrendPlunge: Vec<f64>,

    /// User ellipses dip and strike array.
    ueDipStrike: Vec<f64>,

    /// User ellipses number of points per ellipse array.
    pub uenumPoints: Vec<usize>,

    /// True  - The user is using user defined rectangles.
    /// False - No user defined rectangles are being used.
    pub userRectanglesOnOff: bool,

    /// True  - User rectangles defined by coordinates are being used.
    /// False - No rectangles defined by coordinates are being used.
    pub userRecByCoord: bool,

    /// True  - User ellipses defined by coordinates are being used.
    /// False - No ellpsies defined by coordinates are being used.
    pub userEllByCoord: bool,

    /// True  - User polygons defined by coordinates are being used.
    /// False - No polygons defined by coordinates are being used.
    pub userPolygonByCoord: bool,

    /// Caution: Can create very large files.
    /// Outputs all fractures which were generated during
    /// DFN generation (Accepted + Rejected).
    ///     False: Do not output all radii file.
    ///     True:  Include file of all raddii, acepted + rejected fractures,
    ///            in output files (radii_All.dat).
    pub outputAllRadii: bool,

    /// Number of user defined rectangles.
    pub nUserRect: usize,

    /// User rectangles radii array.
    pub urRadii: Vec<f64>,

    /// All angles from input file for stochastic rectangles are in:
    ///     True  - Degrees
    ///     False - Radians
    pub urAngleOption: bool,

    /// User rectangles beta array.
    pub urBeta: Vec<f64>,

    /// User rectangles aspect ratio array.
    pub uraspect: Vec<f64>,

    /// User rectangles translation array.
    pub urtranslation: Vec<f64>,

    /// User Orientation Option for rectangles
    /// 0 = normal vector
    /// 1 = trend / plunge
    /// 2 = dip / strike
    userRectOrientationOption: u8,

    /// User rectangles normal vector array.
    pub urnormal: Vec<f64>,

    /// User rectangles trend and plunge array.
    urTrendPlunge: Vec<f64>,

    /// User rectangles dip and strike array.
    urDipStrike: Vec<f64>,

    /// Number of user rectangles defined by coordinates.
    pub nRectByCoord: usize,

    /// Number of user ellipses defined by coordinates.
    pub nEllByCoord: usize,

    /// Number of nodes for user defined ellipses by coordinates
    pub nEllNodes: usize,

    /// Array of rectangle coordiates.
    /// Number of elements = 4 * 3 * nRectByCoord
    pub userRectCoordVertices: Vec<f64>,

    /// Array of ellipse coordiates.
    /// Number of elements =  3 * nEllNodes * nEllByCoord
    pub userEllCoordVertices: Vec<f64>,

    /// Name of userPolygon File
    pub polygonFile: String,

    /// If a fracture is rejected, it will be re-translated
    /// to a new position this number of times.
    ///
    /// This helps hit distribution targets for stochastic families
    /// families (Set to 1 to ignore this feature)
    pub rejectsPerFracture: usize,

    /// Z - layers in the DFN
    /// Number of layers defined.
    numOfLayers: usize,

    /// Array of layers:
    /// e.g. {+z1, -z1, +z2, -z2, ... , +zn, -zn}
    pub layers: Vec<f64>,

    /// Array of volumes for each defined layer, in the same order
    /// which layers were listed.
    pub layerVol: Vec<f64>,

    /// Defines which domain, or layer, the family belongs to.
    /// Layer 0 is the entire domain ('domainSize').
    /// Layers numbered > 0 correspond to layers defined above (see 'Layers:').
    /// 1 correspond to the first layer listed, 2 is the next layer listed, etc
    rLayer: Vec<usize>,

    /// Defines which domain, or layer, the family belongs to.
    /// Layer 0 is the entire domain ('domainSize').
    /// Layers numbered > 0 correspond to layers defined above (see 'Layers:').
    /// 1 correspond to the first layer listed, 2 is the next layer listed, etc
    eLayer: Vec<usize>,

    // Regions in the DFN
    /// Number of regions defined.
    numOfRegions: usize,

    /// Array of regions:
    /// e.g. {+z1, -z1, +z2, -z2, ... , +zn, -zn}
    pub regions: Vec<f64>,

    /// Array of volumes for each defined layer, in the same order
    /// which regions were listed.
    pub regionVol: Vec<f64>,

    /// Defines which domain, or regions, the family belongs to.
    /// Regions 0 is the entire domain ('domainSize').
    /// regions numbered > 0 correspond to regions defined above (see 'regions:').
    /// 1 correspond to the first layer listed, 2 is the next layer listed, etc
    rRegion: Vec<usize>,

    /// Defines which domain, or regions, the family belongs to.
    /// Layer 0 is the entire domain ('domainSize').
    /// regions numbered > 0 correspond to regions defined above (see 'regions:').
    /// 1 correspond to the first layer listed, 2 is the next layer listed, etc
    eRegion: Vec<usize>,

    /// flag if the domain is pruned down to a final domain size
    /// bool polygonBoundaryFlag = false;
    pub polygonBoundaryFlag: bool,

    /// Number of points on the 2D boundary of the polygon domain
    pub numOfDomainVertices: usize,

    /// Vector of points defining the 2D boundary of the domain polygon
    pub domainVertices: Vec<Point3<f64>>,

    /// Global boolean array. Used with stopCondition = 1, P32 option.
    /// Number of elements is equal to the number of stochastic shape families.
    /// Elements correspond to families in the same order of the famProb array.
    /// Elements are initialized to false, and are set to true once the families p32
    /// requirement is met.
    /// Once all elements have values all set to true, all families have had their
    /// P32 requirement
    pub p32Status: Vec<bool>,

    /// False - Use boundaryFaces option.
    /// True  - Ignore boundaryFaces option, keep all clusters
    ///         and remove fractures with no intersections
    pub ignoreBoundaryFaces: bool,
}

/// Reads in all input variables.
/// Creates Shape structure array from user input if
/// using stochastic fracture families.
///
/// # Arguments
///
/// * `input` - Path to input file
pub fn read_input(input: &str) -> (Input, Vec<Shape>) {
    let mut input_var = Input::default();
    let mut shape_family = Vec::new();

    let mut e_angle_option = false;
    let mut e_angle1: Vec<f64> = Vec::new();
    let mut e_angle2: Vec<f64> = Vec::new();
    let mut e_orientation = Vec::new();
    let mut r_angle_option = false;
    let mut r_angle1: Vec<f64> = Vec::new();
    let mut r_angle2: Vec<f64> = Vec::new();
    let mut r_orientation = Vec::new();

    println!("DFN Generator Input File: {}\n", input);
    let mut input_file = File::open(input).unwrap();

    macro_rules! input_var {
        ($var_name:ident) => {
            search_var(&mut input_file, &format!("{}:", stringify!($var_name)));
            input_var.$var_name.read_from_text(&mut input_file);
        };

        ($var_name:ident,$n:expr) => {
            search_var(&mut input_file, &format!("{}:", stringify!($var_name)));
            let mut tmp: Vec<_> = Vec::new();
            tmp.read_from_text(&mut input_file);
            for i in 0..input_var.$var_name.len() {
                input_var.$var_name[i] = tmp[i];
            }
        };
    }

    input_var!(stopCondition);
    input_var!(printRejectReasons);
    input_var!(domainSize, 3);
    input_var!(numOfLayers);

    if input_var.numOfLayers > 0 {
        let mut layers: Vec<f64> = Vec::with_capacity(input_var.numOfLayers * 2); // Multiply by 2 for +z and -z for each layer
        let mut layer_vol = vec![0.; input_var.numOfLayers];

        search_var(&mut input_file, "layers:");
        println!("Number of Layers: {}", input_var.numOfLayers);

        layers.read_from_text(&mut input_file);

        for (i, vol) in layer_vol.iter_mut().enumerate().take(input_var.numOfLayers) {
            let idx = i * 2;
            print!(
                "    Layer {}{{-z,+z}}: {:?}, Volume: ",
                i + 1,
                &layers[idx..idx + 2]
            );
            *vol = input_var.domainSize[0]
                * input_var.domainSize[1]
                * ((layers[idx + 1] as isize - layers[idx] as isize).abs() as f64);
            println!("{} m^3", vol);
        }

        println!();
    }

    input_var!(numOfRegions);

    if input_var.numOfRegions > 0 {
        search_var(&mut input_file, "regions:");
        let mut regions: Vec<isize> = Vec::with_capacity(input_var.numOfRegions * 6); // Multiply by 6 xmin, xmax, ymin, ymax, zmin, zmax
        regions.read_from_text(&mut input_file);

        println!("Number of Regions: {}", input_var.numOfRegions);

        let mut region_vol: Vec<isize> = Vec::with_capacity(input_var.numOfRegions);
        for i in 0..input_var.numOfRegions {
            let idx = i * 6;
            println!(
                " Region {}: {{-x,+x,-y,+y,-z,+z}}: {:?}",
                i + 1,
                &regions[idx..idx + 6]
            );
            let vol = (regions[idx + 1] - regions[idx]).abs()
                * (regions[idx + 3] - regions[idx + 2]).abs()
                * (regions[idx + 5] - regions[idx + 4]).abs();
            region_vol.push(vol);
            println!(" Volume: {} m^3", vol);
        }

        println!();
    }

    input_var!(h);

    // Set epsilon
    input_var.eps = input_var.h * 1e-8;

    input_var!(disableFram);

    if input_var.disableFram {
        println!("\nFRAM IS DISABLED");
    }

    input_var!(rFram);

    if input_var.rFram {
        println!("\nRunning with relaxed FRAM. Mesh may not be fully conforming");
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
    input_var!(domainSizeIncrease, 3);
    input_var!(keepOnlyLargestCluster);
    input_var!(keepIsolatedFractures);
    input_var!(ignoreBoundaryFaces);
    input_var!(boundaryFaces, 6);
    input_var!(rejectsPerFracture);
    input_var!(nFamRect);
    input_var!(nFamEll);
    input_var!(removeFracturesLessThan);
    input_var!(orientationOption);

    if input_var.orientationOption == 0 {
        println!("Expecting Theta and phi for orientations");
    } else if input_var.orientationOption == 1 {
        println!("Expecting Trend and Plunge for orientations");
    } else if input_var.orientationOption == 2 {
        println!("Expecting Dip and Strike (RHR) for orientations");
    }

    input_var!(polygonBoundaryFlag);

    if input_var.polygonBoundaryFlag {
        println!("Expecting Polygon Boundary for domain edges");
        search_var(&mut input_file, "polygonBoundaryFile:");
        let mut tempstring: String = String::new();
        tempstring.read_from_text(&mut input_file);
        println!("Polygon Boundary File: {}", &tempstring);
        read_domain_vertices(&mut input_var, &tempstring);
        println!(
            "There are {} Vertices on the boundary",
            input_var.numOfDomainVertices
        );

        for i in 0..input_var.numOfDomainVertices {
            println!(
                "Vertex {}: {{{}, {}}}",
                i + 1,
                input_var.domainVertices[i].x,
                input_var.domainVertices[i].y
            );
        }

        println!();
    }

    if input_var.nFamEll > 0 || input_var.nFamRect > 0 {
        input_var!(famProb);
        input_var.famProbOriginal.extend(input_var.famProb.iter());
        input_var!(radiiListIncrease);
    }

    if input_var.nFamEll > 0 {
        input_var!(ebetaDistribution);
        input_var!(eLayer);
        input_var!(eRegion);
        input_var!(edistr);
        input_var!(easpect);
        input_var!(enumPoints);
        search_var(&mut input_file, "angleOption:"); // eAngleOption
        e_angle_option.read_from_text(&mut input_file);

        if input_var.orientationOption == 0 {
            search_var(&mut input_file, "etheta:");
            e_angle1.read_from_text(&mut input_file);
            search_var(&mut input_file, "ephi:");
            e_angle2.read_from_text(&mut input_file);
        } else if input_var.orientationOption == 1 {
            search_var(&mut input_file, "etrend:");
            e_angle1.read_from_text(&mut input_file);
            search_var(&mut input_file, "eplunge:");
            e_angle2.read_from_text(&mut input_file);
        } else if input_var.orientationOption == 2 {
            search_var(&mut input_file, "edip:");
            e_angle1.read_from_text(&mut input_file);
            search_var(&mut input_file, "estrike:");
            e_angle2.read_from_text(&mut input_file);
        }

        input_var!(ebeta);

        // Convert from degrees to radians
        if e_angle_option {
            let _ = e_angle1
                .iter_mut()
                .map(|v| *v *= std::f64::consts::PI / 180.);
            let _ = e_angle2
                .iter_mut()
                .map(|v| *v *= std::f64::consts::PI / 180.);
            let _ = input_var
                .ebeta
                .iter_mut()
                .map(|v| *v *= std::f64::consts::PI / 180.);
        }

        input_var!(ekappa);

        for ((c1, c2), kappa) in zip_eq(zip_eq(&e_angle1, &e_angle2), &input_var.ekappa) {
            e_orientation.push(match input_var.orientationOption {
                0 => Fisher::new_with_theta_phi(*c1, *c2, *kappa, input_var.eps),
                1 => Fisher::new_with_trend_plunge(*c1, *c2, *kappa, input_var.eps),
                2 => Fisher::new_with_dip_strike(*c1, *c2, *kappa, input_var.eps),
                _ => unreachable!(),
            })
        }

        input_var!(eLogMean);
        input_var!(esd);
        input_var!(eExpMean);
        input_var!(emin);
        input_var!(emax);
        input_var!(ealpha);
        input_var!(econst);
        input_var!(eLogMin);
        input_var!(eLogMax);
        input_var!(eExpMin);
        input_var!(eExpMax);

        if input_var.stopCondition == 1 {
            // Get temp array for ellipse p32 targets
            // Used to simplify initialization of shape family structures below
            input_var!(e_p32Targets);
        }
    }

    // Distribution counters
    let mut shape1 = 0; //longnormal dist counter
    let mut shape2 = 0; //truncated power-law dist counter
    let mut shape3 = 0; //exponential dist counter
    let mut shape4 = 0; //constant dist counter
    let mut beta_count = 0;

    // Create shape structures from data gathered above
    for (i, orientation) in zip_eq(0..input_var.nFamEll, e_orientation) {
        let mut new_shape_fam = Shape {
            shape_family: 0, // shapFam = 0 = ellipse, 1 = rect
            distribution_type: input_var.edistr[i],
            num_points: input_var.enumPoints[i],
            aspect_ratio: input_var.easpect[i],
            orientation: Some(orientation),
            layer: input_var.eLayer[i],
            region: input_var.eRegion[i],
            ..Default::default()
        };

        generate_theta(
            &mut new_shape_fam.theta_list,
            new_shape_fam.aspect_ratio,
            new_shape_fam.num_points,
        );

        if input_var.ebetaDistribution[i] {
            // If constant user defined beta option
            new_shape_fam.beta_distribution = true;
            new_shape_fam.beta = input_var.ebeta[beta_count];
            beta_count += 1;
        } else {
            new_shape_fam.beta_distribution = false;
        }

        // dist options:1 = lognormal, 2= truncated power-law, 3= exponential, 4=constant
        match input_var.edistr[i] {
            // Lognormal
            1 => {
                new_shape_fam.mean = input_var.eLogMean[shape1];
                new_shape_fam.sd = input_var.esd[shape1];
                new_shape_fam.log_min = input_var.eLogMin[shape1];
                new_shape_fam.log_max = input_var.eLogMax[shape1];
                shape1 += 1;
            }

            // Truncated power-law
            2 => {
                new_shape_fam.min = input_var.emin[shape2];
                new_shape_fam.max = input_var.emax[shape2];
                new_shape_fam.alpha = input_var.ealpha[shape2];
                shape2 += 1;
            }

            3 => {
                new_shape_fam.exp_mean = input_var.eExpMean[shape3];
                new_shape_fam.exp_lambda = 1. / input_var.eExpMean[shape3];
                new_shape_fam.exp_min = input_var.eExpMin[shape3];
                new_shape_fam.exp_max = input_var.eExpMax[shape3];
                shape3 += 1;
            }

            4 => {
                new_shape_fam.const_radi = input_var.econst[shape4];
                shape4 += 1;
            }

            _ => {
                unreachable!()
            }
        }

        if input_var.stopCondition == 1 {
            new_shape_fam.p32_target = input_var.e_p32Targets[i];
        }

        // Save shape family to perminant array
        shape_family.push(new_shape_fam);
    }

    if input_var.nFamRect > 0 {
        input_var!(rbetaDistribution, input_var.nFamRect);
        input_var!(rdistr, input_var.nFamRect);
        input_var!(rLayer, input_var.nFamRect);
        input_var!(rRegion, input_var.nFamRect);
        input_var!(raspect, input_var.nFamRect);
        search_var(&mut input_file, "angleOption:"); // rAngleOption
        r_angle_option.read_from_text(&mut input_file);

        if input_var.orientationOption == 0 {
            search_var(&mut input_file, "rtheta:");
            r_angle1.read_from_text(&mut input_file);
            search_var(&mut input_file, "rphi:");
            r_angle2.read_from_text(&mut input_file);
        } else if input_var.orientationOption == 1 {
            search_var(&mut input_file, "rtrend:");
            r_angle1.read_from_text(&mut input_file);
            search_var(&mut input_file, "rplunge:");
            r_angle2.read_from_text(&mut input_file);
        } else if input_var.orientationOption == 2 {
            search_var(&mut input_file, "rdip:");
            r_angle1.read_from_text(&mut input_file);
            search_var(&mut input_file, "rstrike:");
            r_angle2.read_from_text(&mut input_file);
        }

        input_var!(rbeta, input_var.nFamRect);

        if r_angle_option {
            let _ = r_angle1
                .iter_mut()
                .map(|v| *v *= std::f64::consts::PI / 180.);
            let _ = r_angle2
                .iter_mut()
                .map(|v| *v *= std::f64::consts::PI / 180.);
            let _ = input_var
                .rbeta
                .iter_mut()
                .map(|v| *v *= std::f64::consts::PI / 180.);
        }

        input_var!(rkappa, input_var.nFamRect);

        for ((c1, c2), kappa) in zip_eq(zip_eq(&r_angle1, &r_angle2), &input_var.rkappa) {
            r_orientation.push(match input_var.orientationOption {
                0 => Fisher::new_with_theta_phi(*c1, *c2, *kappa, input_var.eps),
                1 => Fisher::new_with_trend_plunge(*c1, *c2, *kappa, input_var.eps),
                2 => Fisher::new_with_dip_strike(*c1, *c2, *kappa, input_var.eps),
                _ => unreachable!(),
            })
        }

        input_var!(rLogMean, input_var.nFamRect);
        input_var!(rsd, input_var.nFamRect);
        input_var!(rmin, input_var.nFamRect);
        input_var!(rmax, input_var.nFamRect);
        input_var!(ralpha, input_var.nFamRect);
        input_var!(rExpMean, input_var.nFamRect);
        input_var!(rconst, input_var.nFamRect);
        input_var!(rLogMin, input_var.nFamRect);
        input_var!(rLogMax, input_var.nFamRect);
        input_var!(rExpMin, input_var.nFamRect);
        input_var!(rExpMax, input_var.nFamRect);

        if input_var.stopCondition == 1 {
            // Get temp array for rectangle p32 targets,
            // Used to simplify initialization of shape family structures below
            input_var!(r_p32Targets, input_var.nFamRect);
        }
    }

    // Set stop condition variables
    if input_var.nFamRect > 0 || input_var.nFamEll > 0 {
        if input_var.stopCondition == 0 {
            // npoly option
            input_var!(nPoly);
        } else if input_var.stopCondition == 1 {
            // Stop program when p32 conditions per fracture family are met
            // Staus array for whether or not the family has reached its p32 requirement
            input_var.p32Status = vec![false; input_var.nFamRect + input_var.nFamEll];
        }
    }

    // Counters, used to place variable into correct array index
    // distribution counters
    shape1 = 0; // Longnormal dist counter
    shape2 = 0; // Truncated power-law dist counter
    shape3 = 0; // Exponential dist counter
    shape4 = 0; // Constant dist counter
    beta_count = 0;

    // Create shape strucutres from data gathered above
    for (i, orientation) in zip_eq(0..input_var.nFamRect, r_orientation) {
        let mut new_shape_fam = Shape {
            shape_family: 1, // shapFam = 0 = ellipse, 1 = rect
            distribution_type: input_var.rdistr[i],
            num_points: 4, // Rectangle
            aspect_ratio: input_var.raspect[i],
            orientation: Some(orientation),
            layer: input_var.rLayer[i],
            region: input_var.rRegion[i],
            ..Default::default()
        };

        if input_var.rbetaDistribution[i] {
            // If constant beta option
            new_shape_fam.beta_distribution = true;
            new_shape_fam.beta = input_var.rbeta[beta_count];
            beta_count += 1;
        } else {
            new_shape_fam.beta_distribution = false;
        }

        // dist options:1 = lognormal, 2= truncated power-law, 3= exponential, 4=constant
        match input_var.rdistr[i] {
            // Lognormal
            1 => {
                new_shape_fam.mean = input_var.rLogMean[shape1];
                new_shape_fam.sd = input_var.rsd[shape1];
                new_shape_fam.log_min = input_var.rLogMin[shape1];
                new_shape_fam.log_max = input_var.rLogMax[shape1];
                shape1 += 1;
            }

            // Truncated power-law
            2 => {
                new_shape_fam.min = input_var.rmin[shape2];
                new_shape_fam.max = input_var.rmax[shape2];
                new_shape_fam.alpha = input_var.ralpha[shape2];
                shape2 += 1;
            }

            3 => {
                new_shape_fam.exp_mean = input_var.rExpMean[shape3];
                new_shape_fam.exp_lambda = 1. / input_var.rExpMean[shape3];
                new_shape_fam.exp_min = input_var.rExpMin[shape3];
                new_shape_fam.exp_max = input_var.rExpMax[shape3];
                shape3 += 1;
            }

            4 => {
                new_shape_fam.const_radi = input_var.rconst[shape4];
                shape4 += 1;
            }

            _ => {
                unreachable!()
            }
        }

        if input_var.stopCondition == 1 {
            new_shape_fam.p32_target = input_var.r_p32Targets[i];
        }

        // Save family to perminant array
        shape_family.push(new_shape_fam);
    }

    input_var!(userEllipsesOnOff);

    if input_var.userEllipsesOnOff {
        search_var(&mut input_file, "UserEll_Input_File_Path:");
        let mut tempstring: String = String::new();
        tempstring.read_from_text(&mut input_file);
        let mut u_ell_file = File::open(&tempstring).unwrap();

        println!("User Defined Ellipses File: {}", &tempstring);
        search_var(&mut u_ell_file, "nUserEll:");
        input_var.nUserEll.read_from_text(&mut u_ell_file);
        search_var(&mut u_ell_file, "Radii:");
        input_var.ueRadii.read_from_text(&mut u_ell_file);
        search_var(&mut u_ell_file, "Aspect_Ratio:");
        input_var.ueaspect.read_from_text(&mut u_ell_file);
        search_var(&mut u_ell_file, "AngleOption:");
        input_var.ueAngleOption.read_from_text(&mut u_ell_file);
        search_var(&mut u_ell_file, "Beta:");
        input_var.ueBeta.read_from_text(&mut u_ell_file);
        search_var(&mut u_ell_file, "Translation:");
        input_var.uetranslation.read_from_text(&mut u_ell_file);
        //
        search_var(&mut u_ell_file, "userOrientationOption:");
        input_var
            .userEllOrientationOption
            .read_from_text(&mut u_ell_file);
        println!(
            "userOrientationOption {}",
            input_var.userEllOrientationOption
        );

        if input_var.userEllOrientationOption == 0 {
            search_var(&mut u_ell_file, "Normal:");
            input_var.uenormal.read_from_text(&mut u_ell_file);
        } else if input_var.userEllOrientationOption == 1 {
            search_var(&mut u_ell_file, "Trend_Plunge:");
            input_var.ueTrendPlunge.read_from_text(&mut u_ell_file);
            // Convert Trend and Plunge into Dip and Strike
            let temp = std::f64::consts::PI / 180.;

            for i in 0..input_var.nUserEll {
                let index1 = i * 2;
                let index2 = i * 3;
                let trend = input_var.ueTrendPlunge[index1] * temp;
                let plunge = input_var.ueTrendPlunge[index1 + 1] * temp;
                input_var.uenormal[index2] = trend.cos() * plunge.cos();
                input_var.uenormal[index2 + 1] = trend.sin() * plunge.cos();
                input_var.uenormal[index2 + 2] = plunge.sin();
            }
        } else if input_var.userEllOrientationOption == 2 {
            search_var(&mut u_ell_file, "Dip_Strike:");
            input_var.ueDipStrike.read_from_text(&mut u_ell_file);
            // Convert Trend and Plunge into Dip and Strike
            let temp = std::f64::consts::PI / 180.;

            for i in 0..input_var.nUserEll {
                let index1 = i * 2;
                let index2 = i * 3;
                // Trend and Plunge
                let dip = input_var.ueDipStrike[index1] * temp;
                let strike = input_var.ueDipStrike[index1 + 1] * temp;
                input_var.uenormal[index2] = dip.sin() * strike.sin();
                input_var.uenormal[index2 + 1] = -dip.sin() * strike.cos();
                input_var.uenormal[index2 + 2] = dip.cos();
            }
        }

        search_var(&mut u_ell_file, "Number_of_Vertices:");
        input_var.uenumPoints.read_from_text(&mut u_ell_file);
    }

    input_var!(userRectanglesOnOff);

    if input_var.userRectanglesOnOff {
        search_var(&mut input_file, "UserRect_Input_File_Path:");
        let mut tempstring: String = String::new();
        tempstring.read_from_text(&mut input_file);
        let mut u_rect_file = File::open(&tempstring).unwrap();
        println!("User Defined Rectangles File: {}", &tempstring);
        search_var(&mut u_rect_file, "nUserRect:");
        input_var.nUserRect.read_from_text(&mut u_rect_file);
        search_var(&mut u_rect_file, "Radii:");
        input_var.urRadii.read_from_text(&mut u_rect_file);
        search_var(&mut u_rect_file, "AngleOption:");
        input_var.urAngleOption.read_from_text(&mut u_rect_file);
        search_var(&mut u_rect_file, "Beta:");
        input_var.urBeta.read_from_text(&mut u_rect_file);
        search_var(&mut u_rect_file, "Aspect_Ratio:");
        input_var.uraspect.read_from_text(&mut u_rect_file);
        search_var(&mut u_rect_file, "Translation:");
        input_var.urtranslation.read_from_text(&mut u_rect_file);
        search_var(&mut u_rect_file, "userOrientationOption:");
        input_var
            .userRectOrientationOption
            .read_from_text(&mut u_rect_file);

        // std::cout << "userRectOrientationOption: " << userRectOrientationOption << " \n";
        if input_var.userRectOrientationOption == 0 {
            search_var(&mut u_rect_file, "Normal:");
            input_var.urnormal.read_from_text(&mut u_rect_file);
        } else if input_var.userRectOrientationOption == 1 {
            search_var(&mut u_rect_file, "Trend_Plunge:");
            input_var.urTrendPlunge.read_from_text(&mut u_rect_file);
            input_var.urnormal.read_from_text(&mut u_rect_file);
            let temp = std::f64::consts::PI / 180.;

            // Convert Trend and Plunge into Dip and Strike
            for i in 0..input_var.nUserRect {
                let index1 = i * 2;
                let index2 = i * 3;
                // Trend and Plunge
                // convert to radians
                let trend = input_var.urTrendPlunge[index1] * temp;
                let plunge = input_var.urTrendPlunge[index1 + 1] * temp;
                input_var.urnormal[index2] = trend.cos() * plunge.cos();
                input_var.urnormal[index2 + 1] = trend.sin() * plunge.cos();
                input_var.urnormal[index2 + 2] = plunge.sin();
            }
        } else if input_var.userRectOrientationOption == 2 {
            search_var(&mut u_rect_file, "Dip_Strike:");
            input_var.urDipStrike.read_from_text(&mut u_rect_file);
            // Convert Dip and Strike into normal vectors
            let temp = std::f64::consts::PI / 180.;

            for i in 0..input_var.nUserRect {
                let index1 = i * 2;
                let index2 = i * 3;
                // dip and strike
                // convert to radians
                let dip = input_var.urDipStrike[index1] * temp;
                let strike = input_var.urDipStrike[index1 + 1] * temp;
                input_var.urnormal[index2] = dip.sin() * strike.sin();
                input_var.urnormal[index2 + 1] = -dip.sin() * strike.cos();
                input_var.urnormal[index2 + 2] = dip.cos();
            }
        }
    }

    input_var!(userEllByCoord);
    input_var!(userRecByCoord);
    input_var!(userPolygonByCoord);

    if (input_var.userRectanglesOnOff || input_var.userRecByCoord)
        && (input_var.userEllipsesOnOff || input_var.userEllByCoord)
    {
        input_var!(insertUserRectanglesFirst);
    } else {
        input_var.insertUserRectanglesFirst = false;
    }

    if !input_var.userPolygonByCoord {
        search_var(&mut input_file, "PolygonByCoord_Input_File_Path:");
        input_var.polygonFile.read_from_text(&mut input_file);
    }

    if input_var.userEllByCoord {
        search_var(&mut input_file, "EllByCoord_Input_File_Path:");
        let mut tempstring: String = String::new();
        tempstring.read_from_text(&mut input_file);
        let mut file = File::open(&tempstring).unwrap();
        println!("User Defined Ellipses by Coordinates File: {}", &tempstring);
        search_var(&mut file, "nEllipses:");
        input_var.nEllByCoord.read_from_text(&mut file);
        search_var(&mut file, "nNodes:");
        input_var.nEllNodes.read_from_text(&mut file);
        search_var(&mut file, "Coordinates:");
        get_cords(
            &mut file,
            &mut input_var.userEllCoordVertices,
            input_var.nEllByCoord,
            input_var.nEllNodes,
        );
    }

    if input_var.userRecByCoord {
        search_var(&mut input_file, "RectByCoord_Input_File_Path:");
        let mut tempstring: String = String::new();
        tempstring.read_from_text(&mut input_file);
        let mut u_coord_file = File::open(&tempstring).unwrap();
        println!(
            "User Defined Rectangles by Coordinates File: {}",
            &tempstring
        );
        search_var(&mut u_coord_file, "nRectangles:");
        input_var.nRectByCoord.read_from_text(&mut u_coord_file);
        search_var(&mut u_coord_file, "Coordinates:");
        get_rect_coords(
            &mut u_coord_file,
            &mut input_var.userRectCoordVertices,
            input_var.nRectByCoord,
        );
    }

    // Error check on stopping parameter
    if input_var.nFamEll + input_var.nFamRect == 0 && input_var.stopCondition != 0 {
        // If no stochastic shapes, use nPoly option with npoly = number of user polygons
        println!("WARNING: You have defined stopCondition = 1 (P32 program stopping condition) but have no stochastic shape families defined. Automatically setting stopCondition to 0 for use with user defined polygons and nPoly.\n");
        input_var.stopCondition = 0;

        if !input_var.userEllipsesOnOff
            && !input_var.userRectanglesOnOff
            && !input_var.userRecByCoord
        {
            panic!("ERROR: All polygon generating options are off or undefined, please check input file for errors.");
        }

        let mut count = 0; // Count of user defined polygons

        if input_var.userEllipsesOnOff {
            count += input_var.nUserEll;
        }

        if input_var.userRectanglesOnOff {
            count += input_var.nUserRect;
        }

        if input_var.userRecByCoord {
            count += input_var.nRectByCoord;
        }

        // Set nPoly to the amount of user defined polygons
        input_var.nPoly = count;
    }

    (input_var, shape_family)
}
