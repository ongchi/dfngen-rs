use std::fs::File;

use parry3d_f64::na::{Point3, Vector3};
use tracing::{error, info, warn};

use super::read_input_functions::{
    get_cords, get_rect_coords, read_domain_vertices, search_var, ReadFromTextFile,
    UserFractureReader,
};

use crate::{io::read_input_functions::InputReader, structures::FractureFamilyOption};

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

    /// True  - The user is using user defined ellipses.
    /// False - No user defined ellipses are being used.
    pub userEllipsesOnOff: bool,

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

    // Regions in the DFN
    /// Number of regions defined.
    numOfRegions: usize,

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

    pub user_defined_ell_fractures: Option<UserDefinedFractures>,
    pub user_defined_rect_fractures: Option<UserDefinedFractures>,
}

/// Reads in all input variables.
/// Creates Shape structure array from user input if
/// using stochastic fracture families.
///
/// # Arguments
///
/// * `input` - Path to input file
pub fn read_input(input_file: &str) -> (Input, FractureFamilyOption) {
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
    input_var!(numOfLayers);

    if input_var.numOfLayers > 0 {
        let mut layers: Vec<f64> = Vec::with_capacity(input_var.numOfLayers * 2); // Multiply by 2 for +z and -z for each layer

        input_reader.read_value("layers:", &mut layers);

        for i in 0..input_var.numOfLayers {
            let idx = i * 2;
            info!(
                "    Layer {}{{-z,+z}}: {:?}, Volume: ",
                i + 1,
                &layers[idx..idx + 2]
            );
            let vol = input_var.domainSize[0]
                * input_var.domainSize[1]
                * ((layers[idx + 1] as isize - layers[idx] as isize).abs() as f64);
            info!("{} m^3", vol);
        }
    }

    input_reader.read_value("numOfRegions:", &mut input_var.numOfRegions);
    if input_var.numOfRegions > 0 {
        let mut regions: Vec<isize> = Vec::with_capacity(input_var.numOfRegions * 6); // Multiply by 6 xmin, xmax, ymin, ymax, zmin, zmax
        input_reader.read_value("regions:", &mut regions);

        info!("Number of Regions: {}", input_var.numOfRegions);

        for i in 0..input_var.numOfRegions {
            let idx = i * 6;
            info!(
                " Region {}: {{-x,+x,-y,+y,-z,+z}}: {:?}",
                i + 1,
                &regions[idx..idx + 6]
            );
            let vol = (regions[idx + 1] - regions[idx]).abs()
                * (regions[idx + 3] - regions[idx + 2]).abs()
                * (regions[idx + 5] - regions[idx + 4]).abs();
            info!(" Volume: {} m^3", vol);
        }
    }

    input_var!(h);

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

    fracture_families.extend(input_reader.read_fracture_family("e"));
    fracture_families.extend(input_reader.read_fracture_family("r"));

    let fracture_family_probabilities = if !fracture_families.is_empty() {
        input_var!(radiiListIncrease);

        if input_var.stopCondition == 0 {
            // npoly option
            input_var!(nPoly);
        } else if input_var.stopCondition == 1 {
            // Stop program when p32 conditions per fracture family are met
            // Staus array for whether or not the family has reached its p32 requirement
            input_var.p32Status = vec![false; fracture_families.len()];
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

    input_var!(userEllipsesOnOff);

    if input_var.userEllipsesOnOff {
        let mut user_ell_file: String = String::new();
        input_reader.read_value("UserEll_Input_File_Path:", &mut user_ell_file);
        info!("User Defined Ellipses File: {}", &user_ell_file);

        input_var.user_defined_ell_fractures = Some(UserDefinedFractures::from_fracture_file(
            &user_ell_file,
            true,
        ));
    }

    input_var!(userRectanglesOnOff);

    if input_var.userRectanglesOnOff {
        let mut user_rect_file: String = String::new();
        input_reader.read_value("UserRect_Input_File_Path:", &mut user_rect_file);
        info!("User Defined Rectangles File: {}", &user_rect_file);

        input_var.user_defined_rect_fractures = Some(UserDefinedFractures::from_fracture_file(
            &user_rect_file,
            false,
        ));
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
        input_reader.read_value(
            "PolygonByCoord_Input_File_Path:",
            &mut input_var.polygonFile,
        );
    }

    if input_var.userEllByCoord {
        let mut ell_file: String = String::new();
        input_reader.read_value("EllByCoord_Input_File_Path:", &mut ell_file);
        let mut file = File::open(&ell_file).unwrap();
        info!("User Defined Ellipses by Coordinates File: {}", &ell_file);
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
        let mut rect_file: String = String::new();
        input_reader.read_value("RectByCoord_Input_File_Path:", &mut rect_file);
        let mut file = File::open(&rect_file).unwrap();
        info!(
            "User Defined Rectangles by Coordinates File: {}",
            &rect_file
        );
        search_var(&mut file, "nRectangles:");
        input_var.nRectByCoord.read_from_text(&mut file);
        search_var(&mut file, "Coordinates:");
        get_rect_coords(
            &mut file,
            &mut input_var.userRectCoordVertices,
            input_var.nRectByCoord,
        );
    }

    // Error check on stopping parameter
    if fracture_families.is_empty() && input_var.stopCondition != 0 {
        // If no stochastic shapes, use nPoly option with npoly = number of user polygons
        warn!("You have defined stopCondition = 1 (P32 program stopping condition) but have no stochastic shape families defined. Automatically setting stopCondition to 0 for use with user defined polygons and nPoly.");
        input_var.stopCondition = 0;

        if !input_var.userEllipsesOnOff
            && !input_var.userRectanglesOnOff
            && !input_var.userRecByCoord
        {
            panic!("ERROR: All polygon generating options are off or undefined, please check input file for errors.");
        }

        let mut count = 0; // Count of user defined polygons

        if let Some(user_ells) = &input_var.user_defined_ell_fractures {
            count += user_ells.n_frac;
        }

        if let Some(user_rects) = &input_var.user_defined_rect_fractures {
            count += user_rects.n_frac;
        }

        if input_var.userRecByCoord {
            count += input_var.nRectByCoord;
        }

        // Set nPoly to the amount of user defined polygons
        input_var.nPoly = count;
    }

    (
        input_var,
        FractureFamilyOption {
            families: fracture_families,
            probabilities: fracture_family_probabilities.clone(),
            original_probabilities: fracture_family_probabilities,
        },
    )
}

#[derive(Debug)]
pub struct UserDefinedFractures {
    pub n_frac: usize,
    pub radii: Vec<f64>,
    pub beta: Vec<f64>,
    pub aspect: Vec<f64>,
    pub translation: Vec<[f64; 3]>,
    pub normal: Vec<Vector3<f64>>,
    pub num_points: Vec<usize>,
}

impl UserDefinedFractures {
    pub fn from_fracture_file(path: &str, is_ell: bool) -> Self {
        let mut frac_reader = if is_ell {
            UserFractureReader::new(path, "nUserEll:")
        } else {
            UserFractureReader::new(path, "nUserRect:")
        };

        let mut radii: Vec<f64> = Vec::new();
        let mut beta: Vec<f64> = Vec::new();
        let mut aspect: Vec<f64> = Vec::new();
        let mut translation: Vec<[f64; 3]> = Vec::new();

        let mut orientation_option: u8 = 0;

        let mut normal: Vec<[f64; 3]> = Vec::new();
        let mut trend_plunge: Vec<[f64; 2]> = Vec::new();
        let mut dip_strike: Vec<[f64; 2]> = Vec::new();

        let mut num_points: Vec<usize> = Vec::new();

        frac_reader.read_vec("Radii:", &mut radii);

        let mut angle_option = String::new();
        frac_reader.read_value("AngleOption:", &mut angle_option);
        let angle_convertion_factor = match angle_option.as_str() {
            "degree" => std::f64::consts::PI / 180.,
            "radian" => 1.,
            _ => panic!("Invalid angle option"),
        };

        frac_reader.read_vec("Beta:", &mut beta);
        beta.iter_mut().for_each(|v| *v *= angle_convertion_factor);

        frac_reader.read_vec("Aspect_Ratio:", &mut aspect);
        frac_reader.read_vec_arr3("Translation:", &mut translation);

        frac_reader.read_value("userOrientationOption:", &mut orientation_option);

        if orientation_option == 0 {
            frac_reader.read_vec_arr3("Normal:", &mut normal);
        } else if orientation_option == 1 {
            frac_reader.read_vec_arr2("Trend_Plunge:", &mut trend_plunge);

            // Convert Trend and Plunge into normal vectors
            for [trend, plunge] in &trend_plunge {
                let trend = trend * angle_convertion_factor;
                let plunge = plunge * angle_convertion_factor;
                normal.push([
                    trend.cos() * plunge.cos(),
                    trend.sin() * plunge.cos(),
                    plunge.sin(),
                ]);
            }
        } else if orientation_option == 2 {
            frac_reader.read_vec_arr2("Dip_Strike:", &mut dip_strike);

            // Convert Dip and Strike into normal vectors
            for [dip, strike] in &dip_strike {
                let dip = dip * angle_convertion_factor;
                let strike = strike * angle_convertion_factor;
                normal.push([
                    dip.sin() * strike.sin(),
                    -dip.sin() * strike.cos(),
                    dip.cos(),
                ]);
            }
        }

        if is_ell {
            frac_reader.read_vec("Number_of_Vertices:", &mut num_points);
        }

        UserDefinedFractures {
            n_frac: frac_reader.n_frac,
            radii,
            beta,
            aspect,
            translation,
            normal: normal
                .into_iter()
                .map(|norm| Vector3::from_row_slice(&norm).normalize())
                .collect(),
            num_points,
        }
    }
}
