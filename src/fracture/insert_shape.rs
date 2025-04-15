use std::cell::RefCell;
use std::rc::Rc;

use parry3d_f64::bounding_volume::Aabb;
use parry3d_f64::na::Vector3;
use rand::distr::Uniform;
use rand::Rng;
use rand_mt::Mt64;

use crate::distribution::{TruncExp, TruncLogNormal, TruncPowerLaw};
use crate::structures::{RadiusFunction, ShapeFamily};
use crate::{
    computational_geometry::{apply_rotation2_d, apply_rotation3_d, translate},
    distribution::generating_points::random_translation,
    structures::{Poly, Shape},
};

fn generate_radius(
    h: f64,
    n_fam_ell: usize,
    shape_fam: &mut Shape,
    generator: Rc<RefCell<Mt64>>,
    family_index: isize,
    use_list: bool,
) -> f64 {
    let radius_distribution = shape_fam.radius_distribution.as_ref().unwrap();
    match radius_distribution.function {
        RadiusFunction::Constant(c) => c,
        ref others => {
            if shape_fam.radii_idx >= shape_fam.radii_list.len() || !use_list {
                // If out of radii from list, insert random radius
                let distr = match others {
                    RadiusFunction::LogNormal { mu, sigma } => {
                        let min_val = f64::max(h, radius_distribution.min);
                        let distr =
                            TruncLogNormal::new(min_val, radius_distribution.max, *mu, *sigma)
                                .unwrap();
                        generator.clone().borrow_mut().sample(distr)
                    }
                    RadiusFunction::TruncatedPowerLaw { alpha } => {
                        let distr = TruncPowerLaw::new(
                            radius_distribution.min,
                            radius_distribution.max,
                            *alpha,
                        );
                        Ok(generator.clone().borrow_mut().sample(distr))
                    }
                    RadiusFunction::Exponential { lambda } => {
                        let min_val = f64::max(h, radius_distribution.min);
                        let distr =
                            TruncExp::new(min_val, radius_distribution.max, *lambda).unwrap();

                        generator.clone().borrow_mut().sample(distr)
                    }
                    RadiusFunction::Constant(_) => unreachable!(),
                };

                match distr {
                            Ok(radius) => radius,
                            Err(_) => panic!(
                                    "distribution for {} family {} has been unable to generate a fracture with radius within set parameters after 1000 consecutive tries.",
                                    &shape_fam.shape_family,
                                    get_family_number(n_fam_ell, family_index, shape_fam.shape_family),
                                )
                        }
            } else {
                // Insert radius from list
                let radius = shape_fam.radii_list[shape_fam.radii_idx];
                shape_fam.radii_idx += 1;
                radius
            }
        }
    }
}

/// Generate Polygon/Fracture
///
/// Generates a polygon based on a stochastic fracture shape family
///
/// NOTE: Function does not create bouding box. The bouding box has to be
///       created after fracture truncation
///
/// # Arguments
///
/// * `h` - Minimum feature size
/// * `n_fam_ell` - Number of ellipse families
/// * `domain_size` - Size of the domain
/// * `domain_size_increase` - Increase in domain size
/// * `layers` - Layers in the domain
/// * `regions` - Regions in the domain
/// * `orientation_option` - Orientation option
/// * `eps` - Epsilon value for floating point comparisons
/// * `shape_fam` - Shape family to generate fracture from
/// * `generator` - Random generator
/// * `family_index` - Index of 'shapeFam' in the shapeFamilies array in main()
/// * `use_list` - Use pre-calculated fracture radii list to pull radii from
///              False - Generate random radii every time (used in dryRun()
///              which estimates number of fractures needed when using
///              p32 option and generates the radii lists)
///
/// # Returns
///
/// Random polygon/fracture based from 'shapeFam'
#[allow(clippy::too_many_arguments)]
pub fn generate_poly(
    h: f64,
    n_fam_ell: usize,
    domain_size: &Vector3<f64>,
    domain_size_increase: &Vector3<f64>,
    layers: &[f64],
    regions: &[f64],
    eps: f64,
    shape_fam: &mut Shape,
    generator: Rc<RefCell<Mt64>>,
    family_index: isize,
    use_list: bool,
) -> Poly {
    let radius = generate_radius(
        h,
        n_fam_ell,
        shape_fam,
        generator.clone(),
        family_index,
        use_list,
    );

    let boundary = poly_boundary(
        domain_size,
        domain_size_increase,
        layers,
        regions,
        shape_fam,
    );
    generate_poly_with_radius(eps, radius, shape_fam, boundary, generator, family_index)
}

pub fn poly_boundary(
    domain_size: &Vector3<f64>,
    domain_size_increase: &Vector3<f64>,
    layers: &[f64],
    regions: &[f64],
    shape_fam: &Shape,
) -> Aabb {
    let (mins, maxs) = if shape_fam.layer == 0 && shape_fam.region == 0 {
        // The family layer is the whole domain
        let mins = (-domain_size - domain_size_increase) / 2.;
        let maxs = (domain_size + domain_size_increase) / 2.;
        (mins, maxs)
    } else if shape_fam.layer > 0 && shape_fam.region == 0 {
        // Family belongs to a certain layer, shapeFam.layer is > zero
        // Layers start at 1, but the array of layers start at 0, hence
        // the subtraction by 1
        // Layer 0 is reservered to be the entire domain
        let layer_idx = (shape_fam.layer - 1) * 2;
        let _mins = (-domain_size - domain_size_increase) / 2.;
        let _maxs = (domain_size + domain_size_increase) / 2.;
        // Layers only apply to z coordinates
        let mins = Vector3::new(_mins.x, _mins.y, layers[layer_idx]);
        let maxs = Vector3::new(_maxs.x, _maxs.y, layers[layer_idx + 1]);
        (mins, maxs)
    } else if shape_fam.layer == 0 && shape_fam.region > 0 {
        let region_idx = (shape_fam.region - 1) * 6;
        let mins = Vector3::new(
            regions[region_idx],
            regions[region_idx + 2],
            regions[region_idx + 4],
        );
        let maxs = Vector3::new(
            regions[region_idx + 1],
            regions[region_idx + 3],
            regions[region_idx + 5],
        );
        (mins, maxs)
    } else {
        println!("ERROR!!!");
        println!("Layer and Region both defined for this Family.");
        println!("Exiting Program");
        panic!()
    };

    Aabb::new(mins.into(), maxs.into())
}

/// Generate Polygon/Fracture With Given Radius
///
/// Similar to generatePoly() except the radius is passed to the function.
/// Generates a polygon
/// Shape (ell or rect) still comes from the shapes' familiy
///
/// NOTE: Function does not create bouding box. The bouding box has to be
///       created after fracture truncation
///
/// # Arguments
///
/// * `orientation_option` - Orientation option
/// * `eps` - Epsilon value for floating point comparisons
/// * `radius` - Radius for polygon
/// * `shape_fam` - Shape family to generate fracture from
/// * `bbox` - Bounding box
/// * `generator` - Random generator
/// * `family_index` - Index of 'shapeFam' in the shapeFamilies array in main()
///
/// # Returns
///
/// Polygon with radius passed in arg 1 and shape based on `shape_fam`
pub fn generate_poly_with_radius(
    eps: f64,
    radius: f64,
    shape_fam: &Shape,
    bbox: Aabb,
    generator: Rc<RefCell<Mt64>>,
    family_index: isize,
) -> Poly {
    // New polygon to build
    let mut new_poly = Poly::default();
    // Initialize normal to {0,0,1}. ( All polys start on x-y plane )
    new_poly.normal = Vector3::new(0., 0., 1.);
    new_poly.number_of_nodes = shape_fam.shape_family.number_of_nodes() as isize;
    new_poly.vertices = Vec::with_capacity((3 * new_poly.number_of_nodes) as usize); //numPoints*{x,y,z}
    new_poly.family_num = family_index;

    match shape_fam.shape_family {
        ShapeFamily::Ellipse(n) => {
            initialize_ell_vertices(
                &mut new_poly,
                radius,
                shape_fam.aspect_ratio,
                &shape_fam.theta_list,
                n as usize,
            );
        }
        ShapeFamily::Rectangle => {
            initialize_rect_vertices(&mut new_poly, radius, shape_fam.aspect_ratio);
        }
    }

    // Initialize beta based on distrubution type: 0 = unifrom on [0,2PI], 1 = constant
    let beta = if !shape_fam.beta_distribution {
        let uniform_dist = Uniform::new(0., 2. * std::f64::consts::PI).unwrap();
        generator.borrow_mut().sample(uniform_dist)
    } else {
        shape_fam.beta
    };

    // Apply 2d rotation matrix, twist around origin
    // assumes polygon on x-y plane
    // Angle must be in rad
    apply_rotation2_d(&mut new_poly, beta);
    let mut norm = shape_fam.normal_vector(generator.clone());

    let mag = norm.magnitude();
    if mag < 1. - eps || mag > 1. + eps {
        norm = norm.normalize(); //ensure norm is normalized
    }

    // Rotate vertices to norm (new normal)
    apply_rotation3_d(&mut new_poly, &norm, eps);

    // Save newPoly's new normal vector
    new_poly.normal = norm;

    let t = random_translation(
        generator.clone(),
        bbox.mins.x,
        bbox.maxs.x,
        bbox.mins.y,
        bbox.maxs.y,
        bbox.mins.z,
        bbox.maxs.z,
    );

    // Translate - will also set translation vector in poly structure
    translate(&mut new_poly, t);

    new_poly
}

/// Initialize Rectangular Vertices
///
/// Initializes vertices for rectangular poly using radius (1/2 x length)
/// and aspcet ratio. (xradius = radius, yradius = radius * aspectRatio)
/// Poly will be on x-y plane
///
/// # Arguments
///
/// * `new_poly` - Polygon to initialize vertices
/// * `radius` - Radius (1/2 x dimension length)
/// * `aspect_ratio` - Aspect ratio
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

/// Initialize Ellipse Vertices
///
/// Initializes ellipse vertices on x-y plane
///
/// # Arguments
///
/// * `new_poly` - Polygon to initialize
/// * `radius` - Radius (xradius = radius. yradius = radius * aspectRatio)
/// * `aspect_ratio` - Aspect ratio
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

/// Retranslate Polygon
///
/// Re-translate poly
/// Re-initializes/re-builds (if needed) polygon at origin and translates to
/// new position preserving its size, shape, and normal vector
/// This helps hit target distributions since we reject less
///
/// # Arguments
///
/// * `eps` - Epsilon value for floating point comparisons
/// * `domain_size` - Size of the domain
/// * `domain_size_increase` - Increase in domain size
/// * `layers` - Layers in the domain
/// * `regions` - Regions in the domain
/// * `new_poly` - Polygon
/// * `shape_fam` - Shape family structure which Polygon belongs to
/// * `generator` - Random Generator
#[allow(clippy::too_many_arguments)]
pub fn re_translate_poly(
    eps: f64,
    domain_size: &Vector3<f64>,
    domain_size_increase: &Vector3<f64>,
    layers: &[f64],
    regions: &[f64],
    new_poly: &mut Poly,
    shape_fam: &Shape,
    generator: Rc<RefCell<Mt64>>,
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
        let bbox = poly_boundary(
            domain_size,
            domain_size_increase,
            layers,
            regions,
            shape_fam,
        );
        let t = random_translation(
            generator.clone(),
            bbox.mins.x,
            bbox.maxs.x,
            bbox.mins.y,
            bbox.maxs.y,
            bbox.mins.z,
            bbox.maxs.z,
        );

        // Translate - will also set translation vector in poly structure
        translate(new_poly, t);
    } else {
        // Poly was truncated, need to rebuild the polygon
        new_poly.vertices =
            Vec::with_capacity(shape_fam.shape_family.number_of_nodes() as usize * 3);
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
        new_poly.number_of_nodes = shape_fam.shape_family.number_of_nodes() as isize;

        match shape_fam.shape_family {
            ShapeFamily::Ellipse(n) => {
                initialize_ell_vertices(
                    new_poly,
                    new_poly.xradius,
                    shape_fam.aspect_ratio,
                    &shape_fam.theta_list,
                    n as usize,
                );
            }
            ShapeFamily::Rectangle => {
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
            }
        }

        // Save newPoly's previous normal vector and then reset poly normal to {0,0,1} for applyRotation3D function
        let normal_b = new_poly.normal;
        new_poly.normal.x = 0.;
        new_poly.normal.y = 0.;
        new_poly.normal.z = 1.;

        // Initialize beta based on distrubution type: 0 = unifrom on [0,2PI], 1 = constant
        let beta = if !shape_fam.beta_distribution {
            // Uniform distribution
            let uniform_dist = Uniform::new(0., 2. * std::f64::consts::PI).unwrap();
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
        apply_rotation3_d(new_poly, &normal_b, eps);
        new_poly.normal = normal_b;
        // Translate to new position
        // Translate() will also set translation vector in poly structure
        let bbox = poly_boundary(
            domain_size,
            domain_size_increase,
            layers,
            regions,
            shape_fam,
        );
        let t = random_translation(
            generator.clone(),
            bbox.mins.x,
            bbox.maxs.x,
            bbox.mins.y,
            bbox.maxs.y,
            bbox.mins.z,
            bbox.maxs.z,
        );

        translate(new_poly, t);
    }
}

/// Print Rejection Reson
///
/// Function prints fracture rejection reasons to user based on reject code
/// Currtenly only used for user defined fractures.
///
/// # Arguments
///
/// * `reject_code` - Rejection code
/// * `new_poly` - Poly which was rejected
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

/// Get Family Number
///
/// Turns the global family number into a number the user can more
/// easily understand. The shapeFamily array contains both rectangle
/// and ellipse families.
/// For example: If shapeFamilies array has 3 families of ellipse
/// families and 3 families of rectangle families: {ell, ell, ell, rect, rect, rect}
/// If we want the local family number for index 3, it will return family 1, meaning
/// the first rectangular family. This function is used in conjuntion with shapeType().
///
/// # Arguments
///
/// * `n_fam_ell` - Number of ellipse families
/// * `family_index` - Family index which family belings to in main()'s 'shapeFamilies' array
pub fn get_family_number(
    n_fam_ell: usize,
    family_index: isize,
    family_shape: ShapeFamily,
) -> isize {
    match family_shape {
        ShapeFamily::Ellipse(_) => family_index + 1,
        ShapeFamily::Rectangle => family_index - n_fam_ell as isize + 1,
    }
}
