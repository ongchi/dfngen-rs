use std::cell::RefCell;
use std::rc::Rc;

use rand::distributions::Uniform;
use rand::Rng;
use rand_distr::LogNormal;
use rand_mt::Mt19937GenRand64;

use crate::distribution::Distribution;
use crate::distribution::generating_points::truncated_power_law;
use crate::io::input::Input;
use crate::{
    computational_geometry::{apply_rotation2_d, apply_rotation3_d, translate},
    distribution::generating_points::{fisher_distribution, random_translation},
    structures::{Poly, Shape},
};

// **************************************************************************
// *****************  Generate Polygon/Fracture  ****************************
// Generates a polygon based on a stochastic fracture shape family
// NOTE: Function does not create bouding box. The bouding box has to be
//       created after fracture truncation
// Arg 1: Shape family to generate fracture from
// Arg 2: Random generator, see std <random> c++ library
// Arg 3: Distributions class, currently used only for exponential dist.
// Arg 4: Index of 'shapeFam' (arg 1) in the shapeFamilies array in main()
// Arg 5: True - Use pre-calculated fracture radii list to pull radii from
//        False - Generate random radii every time (used in dryRun()
//                which estimates number of fractures needed when using
//                p32 option and generates the radii lists)
// Return: Random polygon/fracture based from 'shapeFam'
pub fn generate_poly(
    global: &Input,
    shape_fam: &mut Shape,
    generator: Rc<RefCell<Mt19937GenRand64>>,
    distributions: Rc<RefCell<Distribution>>,
    family_index: isize,
    use_list: bool,
) -> Poly {
    // New polygon to build
    let mut new_poly = Poly::default();
    // Initialize normal to {0,0,1}. ( All polys start on x-y plane )
    new_poly.normal[0] = 0.; // x
    new_poly.normal[1] = 0.; // y
    new_poly.normal[2] = 1.; // z
                             // Assign number of nodes
    new_poly.number_of_nodes = shape_fam.num_points as isize;
    new_poly.vertices = Vec::with_capacity((3 * new_poly.number_of_nodes) as usize); //numPoints*{x,y,z}
                                                                                   // Assign family number (index of array)
    new_poly.family_num = family_index;

    // Switch based on distribution type
    match shape_fam.distribution_type {
        1 => {
            // Lognormal
            let mut radius;
            let mut count = 1;

            if shape_fam.radii_idx >= shape_fam.radii_list.len() || !use_list {
                // If out of radii from list, insert random radius
                let log_distribution = LogNormal::new(shape_fam.mean, shape_fam.sd).unwrap();

                loop {
                    radius = generator.clone().borrow_mut().sample(log_distribution);

                    if count % 1000 == 0 {
                        println!(
                    "\nWarning: Lognormal distribution for {} family {} has been  unable to generate a fracture with radius within set parameters after {} consecutive tries.", shape_type(shape_fam), get_family_number(global, family_index, shape_fam.shape_family), count 
                        );
                        println!("Consider adjusting the lognormal paramerters for this family in the input file.");
                        break;
                    }

                    count += 1;

                    if !(radius < global.h
                        || radius < shape_fam.log_min
                        || radius > shape_fam.log_max)
                    {
                        break;
                    }
                }
            } else {
                // Insert radius from list
                radius = shape_fam.radii_list[shape_fam.radii_idx];
                shape_fam.radii_idx += 1;
            }

            if shape_fam.shape_family == 1 {
                // Rectangle
                // Initialize rectangles vertices using lognormal dist.
                initialize_rect_vertices(&mut new_poly, radius, shape_fam.aspect_ratio);
            } else {
                // Ellipse
                initialize_ell_vertices(
                    &mut new_poly,
                    radius,
                    shape_fam.aspect_ratio,
                    &shape_fam.theta_list,
                    shape_fam.num_points,
                );
            }
        }

        2 => {
            // Truncated power-law
            let radius;

            if shape_fam.radii_idx >= shape_fam.radii_list.len() || !use_list {
                // If out of radii from list, generate random radius
                let uniform_dist = Uniform::new(0., 1.);
                radius = truncated_power_law(
                    generator.clone().borrow_mut().sample(uniform_dist),
                    shape_fam.min,
                    shape_fam.max,
                    shape_fam.alpha,
                );
            } else {
                // Pull radius from list
                radius = shape_fam.radii_list[shape_fam.radii_idx];
                shape_fam.radii_idx += 1;
            }

            if shape_fam.shape_family == 1 {
                initialize_rect_vertices(&mut new_poly, radius, shape_fam.aspect_ratio);
            } else {
                initialize_ell_vertices(
                    &mut new_poly,
                    radius,
                    shape_fam.aspect_ratio,
                    &shape_fam.theta_list,
                    shape_fam.num_points,
                );
            }
        }

        3 => {
            // Exponential
            let mut radius;
            let mut count = 1;

            if shape_fam.radii_idx >= shape_fam.radii_list.len() || !use_list {
                // If out of radii from list, generate random radius
                loop {
                    radius = distributions
                        .clone()
                        .borrow_mut()
                        .exp_dist
                        .get_value_by_min_max_val(
                            shape_fam.exp_lambda,
                            shape_fam.min_dist_input,
                            shape_fam.max_dist_input,
                        );

                    if count % 1000 == 0 {
                        println!(
                            "\nWARNING: Exponential distribution for {} family {} has been unable to generate a fracture with radius within set parameters after {} consecutive tries.",
                            shape_type(shape_fam),
                            get_family_number(global, family_index, shape_fam.shape_family),
                            count 
                        );
                        println!(
                            "Consider adjusting the exponential parameters for this family in the input file."
                        );
                        break;
                    }

                    count += 1;

                    if !(radius < global.h
                        || radius < shape_fam.exp_min
                        || radius > shape_fam.exp_max)
                    {
                        break;
                    }
                }
            } else {
                // Insert radius from list
                radius = shape_fam.radii_list[shape_fam.radii_idx];
                shape_fam.radii_idx += 1;
            }

            if shape_fam.shape_family == 1 {
                // Rectangle
                // Initialize rectangles vertices using exp. dist.
                initialize_rect_vertices(&mut new_poly, radius, shape_fam.aspect_ratio);
            } else {
                // Ellipse
                initialize_ell_vertices(
                    &mut new_poly,
                    radius,
                    shape_fam.aspect_ratio,
                    &shape_fam.theta_list,
                    shape_fam.num_points,
                );
            }
        }

        4 => {
            // Constant
            if shape_fam.shape_family == 1 {
                // Rectangle
                // Initialize rectangles vertices
                initialize_rect_vertices(&mut new_poly, shape_fam.const_radi, shape_fam.aspect_ratio);
            } else {
                // Ellipse
                initialize_ell_vertices(
                    &mut new_poly,
                    shape_fam.const_radi,
                    shape_fam.aspect_ratio,
                    &shape_fam.theta_list,
                    shape_fam.num_points,
                );
            }
        }

        _ => unreachable!(),
    }

    // Initialize beta based on distrubution type: 0 = unifrom on [0,2PI], 1 = constant
    let beta = 
    if !shape_fam.beta_distribution {
        //uniform distribution
        let uniform_dist = Uniform::new(0., 2. * std::f64::consts::PI);
        generator.clone().borrow_mut().sample(uniform_dist)
    } else {
        shape_fam.beta
    };

    // Apply 2d rotation matrix, twist around origin
    // Assumes polygon on x-y plane
    // Angle must be in rad
    apply_rotation2_d(&mut new_poly, beta);
    // Fisher distribution / get normal vector
    let mut norm = fisher_distribution(
        global,
        shape_fam.angle_one,
        shape_fam.angle_two,
        shape_fam.kappa,
        generator.clone(),
    );
    let mag = norm.magnitude();

    if mag < 1. - global.eps || mag > 1. + global.eps {
        norm = norm.normalize(); // Ensure norm is normalized
    }

    apply_rotation3_d(global, &mut new_poly, &norm); // Rotate vertices to norm (new normal)
                                                   // Save newPoly's new normal vector
    new_poly.normal[0] = norm[0];
    new_poly.normal[1] = norm[1];
    new_poly.normal[2] = norm[2];
    let t;

    // HERE
    if shape_fam.layer == 0 && shape_fam.region == 0 {
        // The family layer is the whole domain
        t = random_translation(
            generator.clone(),
            (-global.domainSize[0] - global.domainSizeIncrease[0]) / 2.,
            (global.domainSize[0] + global.domainSizeIncrease[0]) / 2.,
            (-global.domainSize[1] - global.domainSizeIncrease[1]) / 2.,
            (global.domainSize[1] + global.domainSizeIncrease[1]) / 2.,
            (-global.domainSize[2] - global.domainSizeIncrease[2]) / 2.,
            (global.domainSize[2] + global.domainSizeIncrease[2]) / 2.,
        );
    } else if shape_fam.layer > 0 && shape_fam.region == 0 {
        // Family belongs to a certain layer, shapeFam.layer is > zero
        // Layers start at 1, but the array of layers start at 0, hence
        // the subtraction by 1
        // Layer 0 is reservered to be the entire domain
        let layer_idx = (shape_fam.layer - 1) * 2;
        // Layers only apply to z coordinates
        t = random_translation(
            generator.clone(),
            (-global.domainSize[0] - global.domainSizeIncrease[0]) / 2.,
            (global.domainSize[0] + global.domainSizeIncrease[0]) / 2.,
            (-global.domainSize[1] - global.domainSizeIncrease[1]) / 2.,
            (global.domainSize[1] + global.domainSizeIncrease[1]) / 2.,
            global.layers[layer_idx],
            global.layers[layer_idx + 1],
        );
    } else if shape_fam.layer == 0 && shape_fam.region > 0 {
        let region_idx = (shape_fam.region - 1) * 6;
        // Layers only apply to z coordinates
        t = random_translation(
            generator.clone(),
            global.regions[region_idx],
            global.regions[region_idx + 1],
            global.regions[region_idx + 2],
            global.regions[region_idx + 3],
            global.regions[region_idx + 4],
            global.regions[region_idx + 5],
        );
    } else {
        // t = randomTranslation(generator, -1., 1., -1., 1., -1., 1.);
        panic!("ERROR!!!\nLayer and Region both defined for this Family.\nExiting Program\n");
    }

    // Translate - will also set translation vector in poly structure
    translate(&mut new_poly, &t);

    new_poly
}

// **************************************************************************
// *************  Generate Polygon/Fracture With Given Radius  **************
// Similar to generatePoly() except the radius is passed to the function.
// Generates a polygon
// Shape (ell or rect) still comes from the shapes' familiy
// NOTE: Function does not create bouding box. The bouding box has to be
//       created after fracture truncation
// Arg 1: Radius for polygon
// Arg 2: Shape family to generate fracture from
// Arg 3: Random generator, see std <random> c++ library
// Arg 4: Distributions class, currently used only for exponential dist.
// Arg 5: Index of 'shapeFam' (arg 1) in the shapeFamilies array in main()
// Return: Polygond with radius passed in arg 1 and shape based on 'shapeFam'
pub fn generate_poly_with_radius(
    global: &Input,
    radius: f64,
    shape_fam: &Shape,
    generator: Rc<RefCell<Mt19937GenRand64>>,
    family_index: isize,
) -> Poly {
    // New polygon to build
    let mut new_poly = Poly::default();
    // Initialize normal to {0,0,1}. ( All polys start on x-y plane )
    new_poly.normal[0] = 0.; //x
    new_poly.normal[1] = 0.; //y
    new_poly.normal[2] = 1.; //z
                             // Assign number of nodes
    new_poly.number_of_nodes = shape_fam.num_points as isize;
    new_poly.vertices = Vec::with_capacity((3 * new_poly.number_of_nodes) as usize); //numPoints*{x,y,z}
                                                                                   // Assign family number (index of shapeFam array)
    new_poly.family_num = family_index;

    // TODO: Convert any degrees to rad
    // in readInput() to avoid continuous checking

    if shape_fam.shape_family == 1 {
        // If rectangle shape
        // Initialize rectangles vertices
        initialize_rect_vertices(&mut new_poly, radius, shape_fam.aspect_ratio);
    } else {
        // Ellipse
        initialize_ell_vertices(
            &mut new_poly,
            radius,
            shape_fam.aspect_ratio,
            &shape_fam.theta_list,
            shape_fam.num_points,
        );
    }


    // Initialize beta based on distrubution type: 0 = unifrom on [0,2PI], 1 = constant
    let beta = if !shape_fam.beta_distribution {
        // Uniform distribution
        let uniform_dist = Uniform::new(0., 2. * std::f64::consts::PI);
        generator.clone().borrow_mut().sample(uniform_dist)
    } else {
        shape_fam.beta
    };

    // Apply 2d rotation matrix, twist around origin
    // assumes polygon on x-y plane
    // Angle must be in rad
    apply_rotation2_d(&mut new_poly, beta);
    // Fisher distribution / get normal vector
    let mut norm = fisher_distribution(
        global,
        shape_fam.angle_one,
        shape_fam.angle_two,
        shape_fam.kappa,
        generator.clone(),
    );
    let mag = norm.magnitude();

    if mag < 1. - global.eps || mag > 1. + global.eps {
        norm = norm.normalize(); //ensure norm is normalized
    }

    apply_rotation3_d(global, &mut new_poly, &norm); // Rotate vertices to norm (new normal)
                                                   // Save newPoly's new normal vector
    new_poly.normal[0] = norm[0];
    new_poly.normal[1] = norm[1];
    new_poly.normal[2] = norm[2];

    let t;

    if shape_fam.layer == 0 && shape_fam.region == 0 {
        // The family layer is the whole domain
        t = random_translation(
            generator.clone(),
            (-global.domainSize[0] - global.domainSizeIncrease[0]) / 2.,
            (global.domainSize[0] + global.domainSizeIncrease[0]) / 2.,
            (-global.domainSize[1] - global.domainSizeIncrease[1]) / 2.,
            (global.domainSize[1] + global.domainSizeIncrease[1]) / 2.,
            (-global.domainSize[2] - global.domainSizeIncrease[2]) / 2.,
            (global.domainSize[2] + global.domainSizeIncrease[2]) / 2.,
        );
    } else if shape_fam.layer > 0 && shape_fam.region == 0 {
        // Family belongs to a certain layer, shapeFam.layer is > zero
        // Layers start at 1, but the array of layers start at 0, hence
        // the subtraction by 1
        // Layer 0 is reservered to be the entire domain
        let layer_idx = (shape_fam.layer - 1) * 2;
        // Layers only apply to z coordinates
        t = random_translation(
            generator.clone(),
            (-global.domainSize[0] - global.domainSizeIncrease[0]) / 2.,
            (global.domainSize[0] + global.domainSizeIncrease[0]) / 2.,
            (-global.domainSize[1] - global.domainSizeIncrease[1]) / 2.,
            (global.domainSize[1] + global.domainSizeIncrease[1]) / 2.,
            global.layers[layer_idx],
            global.layers[layer_idx + 1],
        );
    } else if shape_fam.layer == 0 && shape_fam.region > 0 {
        let region_idx = (shape_fam.region - 1) * 6;
        // Layers only apply to z coordinates
        t = random_translation(
            generator.clone(),
            global.regions[region_idx],
            global.regions[region_idx + 1],
            global.regions[region_idx + 2],
            global.regions[region_idx + 3],
            global.regions[region_idx + 4],
            global.regions[region_idx + 5],
        );
    } else {
        println!("ERROR!!!");
        println!("Layer and Region both defined for this Family.");
        println!("Exiting Program");
        panic!()
    }

    // Translate - will also set translation vector in poly structure
    translate(&mut new_poly, &t);

    new_poly
}

// **************  Initialize Rectangular Vertices  *****************
// Initializes vertices for rectangular poly using radius (1/2 x length)
// and aspcet ratio. (xradius = radius, yradius = radius * aspectRatio)
// Poly will be on x-y plane
// Arg 1: Polygon to initialize vertices
// Arg 2: Radius (1/2 x dimension length)
// Arg 3: Aspect ratio
pub fn initialize_rect_vertices(new_poly: &mut Poly, radius: f64, aspect_ratio: f64) {
    let x = radius;
    let y = radius * aspect_ratio;
    new_poly.xradius = x;
    new_poly.yradius = y;
    new_poly.aspect_ratio = aspect_ratio;
    // Initialize vertices
    new_poly.vertices[0] = x;
    new_poly.vertices[1] = y;
    new_poly.vertices[2] = 0.;
    new_poly.vertices[3] = -x;
    new_poly.vertices[4] = y;
    new_poly.vertices[5] = 0.;
    new_poly.vertices[6] = -x;
    new_poly.vertices[7] = -y;
    new_poly.vertices[8] = 0.;
    new_poly.vertices[9] = x;
    new_poly.vertices[10] = -y;
    new_poly.vertices[11] = 0.;
}

// ****************************************************************
// ************* Initialize Ellipse Vertices **********************
// Initializes ellipse vertices on x-y plane
// Arg 1: Poly to initialize
// Arg 2: Radius (xradius = radius. yradius = radius * aspectRatio)
// Arg 3: Aspect ratio
pub fn initialize_ell_vertices(
    new_poly: &mut Poly,
    radius: f64,
    aspect_ratio: f64,
    theta_list: &[f64],
    num_points: usize,
) {
    new_poly.xradius = radius;
    new_poly.yradius = radius * aspect_ratio;
    new_poly.aspect_ratio = aspect_ratio;

    for (_, theta) in theta_list.iter().enumerate().take(num_points) {
        new_poly.vertices.push(radius * theta.cos());
        new_poly.vertices.push(radius * aspect_ratio * theta.sin());
        new_poly.vertices.push(0.);
    }
}

// **********************************************************************
// ********************  Retranslate Polygon  ***************************
// Re-translate poly
// Re-initializes/re-builds (if needed) polygon at origin and translates to
// new position preserving its size, shape, and normal vector
// This helps hit target distributions since we reject less
// Arg 1: Polygon
// Arg 2: Shape family structure which Polygon belongs to
// Arg 3: Random Generator
pub fn re_translate_poly(
    global: &Input,
    new_poly: &mut Poly,
    shape_fam: &Shape,
    generator: Rc<RefCell<Mt19937GenRand64>>,
) {
    if !new_poly.truncated {
        // If poly isn't truncated we can skip a lot of steps such
        // as reallocating vertice memory, rotations, etc..
        new_poly.group_num = 0; // Clear cluster group information
        new_poly.intersection_index.clear(); // Clear any saved intersections

        // Move poly back to origin
        for i in 0..new_poly.number_of_nodes {
            let idx = (3 * i) as usize;
            new_poly.vertices[idx] -= new_poly.translation[0]; // x
            new_poly.vertices[idx + 1] -= new_poly.translation[1]; // y
            new_poly.vertices[idx + 2] -= new_poly.translation[2]; // z
        }

        // Translate to new position
        let t;

        if shape_fam.layer == 0 && shape_fam.region == 0 {
            // The family layer is the whole domain
            t = random_translation(
                generator.clone(),
                (-global.domainSize[0] - global.domainSizeIncrease[0]) / 2.,
                (global.domainSize[0] + global.domainSizeIncrease[0]) / 2.,
                (-global.domainSize[1] - global.domainSizeIncrease[1]) / 2.,
                (global.domainSize[1] + global.domainSizeIncrease[1]) / 2.,
                (-global.domainSize[2] - global.domainSizeIncrease[2]) / 2.,
                (global.domainSize[2] + global.domainSizeIncrease[2]) / 2.,
            );
        } else if shape_fam.layer > 0 && shape_fam.region == 0 {
            // Family belongs to a certain layer, shapeFam.layer is > zero
            // Layers start at 1, but the array of layers start at 0, hence
            // the subtraction by 1
            // Layer 0 is reservered to be the entire domain
            let layer_idx = (shape_fam.layer - 1) * 2;
            // Layers only apply to z coordinates
            t = random_translation(
                generator.clone(),
                (-global.domainSize[0] - global.domainSizeIncrease[0]) / 2.,
                (global.domainSize[0] + global.domainSizeIncrease[0]) / 2.,
                (-global.domainSize[1] - global.domainSizeIncrease[1]) / 2.,
                (global.domainSize[1] + global.domainSizeIncrease[1]) / 2.,
                global.layers[layer_idx],
                global.layers[layer_idx + 1],
            );
        } else if shape_fam.layer == 0 && shape_fam.region > 0 {
            let region_idx = (shape_fam.region - 1) * 6;
            // Layers only apply to z coordinates
            t = random_translation(
                generator.clone(),
                global.regions[region_idx],
                global.regions[region_idx + 1],
                global.regions[region_idx + 2],
                global.regions[region_idx + 3],
                global.regions[region_idx + 4],
                global.regions[region_idx + 5],
            );
        } else {
            // you should never get here
            // t = randomTranslation(generator, -1., 1., -1., 1., -1., 1.);
            println!("ERROR!!!");
            println!("Layer and Region both defined for this Family.");
            println!("Exiting Program");
            panic!()
        }

        // Translate - will also set translation vector in poly structure
        translate(new_poly, &t);
    } else {
        // Poly was truncated, need to rebuild the polygon
        new_poly.vertices = Vec::with_capacity(shape_fam.num_points * 3);
        // Reset boundary faces (0 means poly is no longer touching a boundary)
        new_poly.faces[0] = false;
        new_poly.faces[1] = false;
        new_poly.faces[2] = false;
        new_poly.faces[3] = false;
        new_poly.faces[4] = false;
        new_poly.faces[5] = false;
        new_poly.truncated = false; // Set to 0 to mean not truncated
        new_poly.group_num = 0; // Clear cluster group information
        new_poly.intersection_index.clear(); // Clear any saved intersections
        new_poly.number_of_nodes = shape_fam.num_points as isize;

        if shape_fam.shape_family == 1 {
            // 1 is rectanglular families. rebuild rectangle
            // Rebuild poly at origin using previous size
            new_poly.vertices[0] = new_poly.xradius;
            new_poly.vertices[1] = new_poly.yradius;
            new_poly.vertices[2] = 0.;
            new_poly.vertices[3] = -new_poly.xradius;
            new_poly.vertices[4] = new_poly.yradius;
            new_poly.vertices[5] = 0.;
            new_poly.vertices[6] = -new_poly.xradius;
            new_poly.vertices[7] = -new_poly.yradius;
            new_poly.vertices[8] = 0.;
            new_poly.vertices[9] = new_poly.xradius;
            new_poly.vertices[10] = -new_poly.yradius;
            new_poly.vertices[11] = 0.;
        } else {
            // Rebuild ellipse
            initialize_ell_vertices(
                new_poly,
                new_poly.xradius,
                shape_fam.aspect_ratio,
                &shape_fam.theta_list,
                shape_fam.num_points,
            );
        }

        // Save newPoly's previous normal vector and then reset poly normal to {0,0,1} for applyRotation3D function
        let normal_b = new_poly.normal;
        new_poly.normal.x = 0.;
        new_poly.normal.y = 0.;
        new_poly.normal.z = 1.;

        // Initialize beta based on distrubution type: 0 = unifrom on [0,2PI], 1 = constant
        let beta = if !shape_fam.beta_distribution {
            // Uniform distribution
            let uniform_dist = Uniform::new(0., 2. * std::f64::consts::PI);
            generator.clone().borrow_mut().sample(uniform_dist)
        } else {
            // Constant
            shape_fam.beta
        };

        // Apply 2d rotation matrix, twist around origin
        // Assumes polygon on x-y plane
        // Angle must be in rad
        apply_rotation2_d(new_poly, beta);
        // Rotates poly from {0,0,1} to normalB, NEED to save normalB to newPoly.normal afterwards
        apply_rotation3_d(global, new_poly, &normal_b);
        new_poly.normal = normal_b;
        // Translate to new position
        // Translate() will also set translation vector in poly structure
        let t;

        if shape_fam.layer == 0 && shape_fam.region == 0 {
            // The family layer is the whole domain
            t = random_translation(
                generator.clone(),
                (-global.domainSize[0] - global.domainSizeIncrease[0]) / 2.,
                (global.domainSize[0] + global.domainSizeIncrease[0]) / 2.,
                (-global.domainSize[1] - global.domainSizeIncrease[1]) / 2.,
                (global.domainSize[1] + global.domainSizeIncrease[1]) / 2.,
                (-global.domainSize[2] - global.domainSizeIncrease[2]) / 2.,
                (global.domainSize[2] + global.domainSizeIncrease[2]) / 2.,
            );
        } else if shape_fam.layer > 0 && shape_fam.region == 0 {
            // Family belongs to a certain layer, shapeFam.layer is > zero
            // Layers start at 1, but the array of layers start at 0, hence
            // the subtraction by 1
            // Layer 0 is reservered to be the entire domain
            let layer_idx = (shape_fam.layer - 1) * 2;
            // Layers only apply to z coordinates
            t = random_translation(
                generator.clone(),
                (-global.domainSize[0] - global.domainSizeIncrease[0]) / 2.,
                (global.domainSize[0] + global.domainSizeIncrease[0]) / 2.,
                (-global.domainSize[1] - global.domainSizeIncrease[1]) / 2.,
                (global.domainSize[1] + global.domainSizeIncrease[1]) / 2.,
                global.layers[layer_idx],
                global.layers[layer_idx + 1],
            );
        } else if shape_fam.layer == 0 && shape_fam.region > 0 {
            let region_idx = (shape_fam.region - 1) * 6;
            // Layers only apply to z coordinates
            t = random_translation(
                generator.clone(),
                global.regions[region_idx],
                global.regions[region_idx + 1],
                global.regions[region_idx + 2],
                global.regions[region_idx + 3],
                global.regions[region_idx + 4],
                global.regions[region_idx + 5],
            );
        } else {
            // t = randomTranslation(generator, -1., 1., -1., 1., -1., 1.);
            println!("ERROR!!!");
            println!("Layer and Region both defined for this Family.");
            println!("Exiting Program");
            panic!();
        }

        translate(new_poly, &t);
    }
}

// **********************************************************************
// *************  Check if Target P32 Has Been Met  *********************
// For th  P32 stopping option (see input file)
// Function checks the status of the fractures families target P32
// The 'p32Status' array is 1 to 1 with the total number of families
// When a family's P32 requirement has been met, it's corresponding
// element in 'p32Status' is set to '1'
// For example, the status for shapeFamiliy[0] is in p32Status[0].
// If p32Status[0] = 1, shapeFamily[0] has met its p32/fracture intensity target,
// else it has not and will keep inserting more fractures
//
// Arg 1: Size of 'p32Status' array (same size as the shapeFamilies array)
// Return: True once ALL p32 targets have been met for all families
//         False otherwise
pub fn p32_complete(input: &Input, size: usize) -> bool {
    // Check if p32Status array is all 1's, if not return 0
    for i in 0..size {
        if !input.p32Status[i] {
            return false;
        }
    }

    // If function has not returned yet, array is all 1's
    true
}

// ***************************************************************************
// **********************  Print Rejection Reson  ****************************
// Function prints fracture rejection reasons to user based on reject code
// Currtenly only used for user defined fractures.
// Arg 1: Rejection code
// Arg 2: Poly which was rejected
pub fn print_reject_reason(reject_code: i32, new_poly: &Poly) {
    if new_poly.family_num >= 0 {
        println!(
            "Attempted fracture from family {} was rejected:",
            new_poly.family_num
        );
    }

    match reject_code {
        -2 => {
            println!("\trejectCode = -2: Intersection of length < h.");
        }

        -1 => {
            println!("\trejectCode = -1: Fracture too close to a node.\n")
        }

        -6 => {
            println!("\trejectCode = -6: Fracture too close to another fracture's edge.");
        }

        -7 => {
            println!("\trejectCode = -7: Fractures intersecting on same plane");
        }

        -10 => {
            println!("\trejectCode = -10: Rejected triple intersection due to triple intersections being turned off in input file.")
        }

        -11 => {
            println!("\trejectCode = -11: Fracture's intersection landed too close to a previous intersection.");
        }

        -12 => {
            println!("\trejectCode = -12: Fracture created a triple intersection with an angle too small.");
        }

        -13 => {
            println!("\trejectCode = -13: Fracture created a triple intersection with the triple intersection point too close to an intersection's endpoint.");
        }

        -14 => {
            println!("\trejectCode = -14: Fracture created a triple intersection with the triple intersection point too close to another triple intersection point.");
        }

        _ => {
            println!("\trejectCode = {}", reject_code);
        }
    }
}

// ******************************************************************
// *********************  Get Family Number *************************
// Turns the global family number into a number the user can more
// easily understand. The shapeFamily array contains both rectangle
// and ellipse families.
// For example: If shapeFamilies array has 3 families of ellipse
// families and 3 families of rectangle families: {ell, ell, ell, rect, rect, rect}
// If we want the local family number for index 3, it will return family 1, meaning
// the first rectangular family. This function is used in conjuntion with shapeType().
//
// Arg 1: Family index which family belings to
//        in main()'s 'shapeFamilies' array
// Arg 2: Rectangle or ellipse family
//        0 - Ellipse
//        1 - Rectangle
pub fn get_family_number(input: &Input, family_index: isize, family_shape: usize) -> isize {
    if family_shape != 0 {
        // if not ellipse family
        family_index - input.nFamEll as isize + 1
    } else {
        family_index + 1
    }
}

// ******************************************************************
// ********************  Print Shape Type  **************************
// Print type of family (ellipse or rectangle)
// Arg 1: Shape family
pub fn shape_type(shape_fam: &Shape) -> &'static str {
    if shape_fam.shape_family == 0 {
        "Ellipse"
    } else {
        "Rectangular"
    }
}

// ******************************************************************
// ************  Get Max Fracture Radii From Family  ****************
// Returns the largest fracture radii defined by the user for
// a fracture family.
// This function is used for the 'forceLargeFractures' option
// in the input file
// Arg 1: Shape family
// Return: User's maximum radii for shapeFam (arg 1)
pub fn get_largest_fracture_radius(shape_fam: &Shape) -> f64 {
    match shape_fam.distribution_type {
        // Log-normal
        1 => shape_fam.log_max,
        // Power-law
        2 => shape_fam.max,
        // Exponential
        3 => shape_fam.exp_max,
        // Constant
        _ => shape_fam.const_radi,
    }
}
