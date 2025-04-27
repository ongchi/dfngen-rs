use std::cell::RefCell;
use std::rc::Rc;

use rand::distr::Uniform;
use rand::Rng;
use rand_mt::Mt64;
use tracing::warn;

use crate::computational_geometry::{apply_rotation2_d, apply_rotation3_d, translate};
use crate::distribution::generating_points::{generate_theta, random_position};
use crate::fracture::fracture_family::FractureFamily;
use crate::structures::{Poly, Shape};

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
    new_poly.vertices.push(x);
    new_poly.vertices.push(y);
    new_poly.vertices.push(0.);
    new_poly.vertices.push(-x);
    new_poly.vertices.push(y);
    new_poly.vertices.push(0.);
    new_poly.vertices.push(-x);
    new_poly.vertices.push(-y);
    new_poly.vertices.push(0.);
    new_poly.vertices.push(x);
    new_poly.vertices.push(-y);
    new_poly.vertices.push(0.);
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
/// * `frac_fam` - Fracture family structure which Polygon belongs to
/// * `generator` - Random Generator
#[allow(clippy::too_many_arguments)]
pub fn re_translate_poly(
    eps: f64,
    new_poly: &mut Poly,
    frac_fam: &FractureFamily,
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
        let t = random_position(frac_fam.boundary, generator.clone());

        // Translate - will also set translation vector in poly structure
        translate(new_poly, t);
    } else {
        // Poly was truncated, need to rebuild the polygon
        new_poly.vertices = Vec::with_capacity(frac_fam.shape.number_of_nodes() as usize * 3);
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
        new_poly.number_of_nodes = frac_fam.shape.number_of_nodes() as isize;

        match frac_fam.shape {
            Shape::Ellipse(n) => {
                let theta_list = generate_theta(
                    frac_fam.aspect_ratio,
                    match frac_fam.shape {
                        Shape::Ellipse(n) => n as usize,
                        Shape::Rectangle => 4,
                    },
                );
                initialize_ell_vertices(
                    new_poly,
                    new_poly.xradius,
                    frac_fam.aspect_ratio,
                    &theta_list,
                    n as usize,
                );
            }
            Shape::Rectangle => {
                // Rebuild poly at origin using previous size
                new_poly.vertices.push(new_poly.xradius);
                new_poly.vertices.push(new_poly.yradius);
                new_poly.vertices.push(0.);
                new_poly.vertices.push(-new_poly.xradius);
                new_poly.vertices.push(new_poly.yradius);
                new_poly.vertices.push(0.);
                new_poly.vertices.push(-new_poly.xradius);
                new_poly.vertices.push(-new_poly.yradius);
                new_poly.vertices.push(0.);
                new_poly.vertices.push(new_poly.xradius);
                new_poly.vertices.push(-new_poly.yradius);
                new_poly.vertices.push(0.);
            }
        }

        // Save newPoly's previous normal vector and then reset poly normal to {0,0,1} for applyRotation3D function
        let normal_b = new_poly.normal;
        new_poly.normal.x = 0.;
        new_poly.normal.y = 0.;
        new_poly.normal.z = 1.;

        // Initialize beta based on distrubution type: 0 = unifrom on [0,2PI], 1 = constant
        let beta = if !frac_fam.beta_distribution {
            // Uniform distribution
            let uniform_dist = Uniform::new(0., 2. * std::f64::consts::PI).unwrap();
            generator.clone().borrow_mut().sample(uniform_dist)
        } else {
            // Constant
            frac_fam.beta
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
        let t = random_position(frac_fam.boundary, generator.clone());

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
        warn!(
            "Attempted fracture from family {} was rejected:",
            new_poly.family_num
        );
    }

    match reject_code {
        -2 => {
            warn!("\trejectCode = -2: Intersection of length < h.");
        }

        -1 => {
            warn!("\trejectCode = -1: Fracture too close to a node.")
        }

        -6 => {
            warn!("\trejectCode = -6: Fracture too close to another fracture's edge.");
        }

        -7 => {
            warn!("\trejectCode = -7: Fractures intersecting on same plane");
        }

        -10 => {
            warn!("\trejectCode = -10: Rejected triple intersection due to triple intersections being turned off in input file.")
        }

        -11 => {
            warn!("\trejectCode = -11: Fracture's intersection landed too close to a previous intersection.");
        }

        -12 => {
            warn!("\trejectCode = -12: Fracture created a triple intersection with an angle too small.");
        }

        -13 => {
            warn!("\trejectCode = -13: Fracture created a triple intersection with the triple intersection point too close to an intersection's endpoint.");
        }

        -14 => {
            warn!("\trejectCode = -14: Fracture created a triple intersection with the triple intersection point too close to another triple intersection point.");
        }

        _ => {
            warn!("\trejectCode = {}", reject_code);
        }
    }
}

/// Get Family Number
///
/// Turns the global family number into a number the user can more
/// easily understand. The FractureFamily array contains both rectangle
/// and ellipse families.
/// For example: If FractureFamilies array has 3 families of ellipse
/// families and 3 families of rectangle families: {ell, ell, ell, rect, rect, rect}
/// If we want the local family number for index 3, it will return family 1, meaning
/// the first rectangular family. This function is used in conjuntion with shapeType().
///
/// # Arguments
///
/// * `n_fam_ell` - Number of ellipse families
/// * `family_index` - Family index which family belings to in main()'s 'shapeFamilies' array
pub fn get_family_number(n_fam_ell: usize, family_index: isize, family_shape: Shape) -> isize {
    match family_shape {
        Shape::Ellipse(_) => family_index + 1,
        Shape::Rectangle => family_index - n_fam_ell as isize + 1,
    }
}
