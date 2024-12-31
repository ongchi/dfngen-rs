use std::fs::File;

use parry3d_f64::na::{distance, Point3, Vector3};

use super::domain::domain_truncation;
use super::insert_shape::print_reject_reason;

use crate::{
    computational_geometry::{create_bounding_box, intersection_checking},
    io::read_input_functions::{search_var, ReadFromTextFile},
    math_functions::get_area,
    structures::{IntersectionPoints, Poly, RejectedUserFracture, Stats},
};

// **********************************************************************
// **********************************************************************
// Used to read in ellipse coordinates when the user is using
// user ellipses defined by coordinates option.
// Arg 1: ifstream file object
// Arg 2: OUTPUT. Pointer to array to store the coordinates
// Arg 3: Number of ellipses
// Arg 4: Number of points per ellipse
fn get_poly_coords(stream: &mut File, n_vertices: usize) -> Vec<f64> {
    let mut vertices = Vec::with_capacity(n_vertices * 3);

    for _ in 0..n_vertices {
        let mut tmp = [0., 0., 0.];
        tmp.read_from_text(stream);
        vertices.extend(tmp);
    }

    vertices
}

fn create_poly(file: &mut File) -> Poly {
    let mut n_poly_nodes: usize = 0;
    n_poly_nodes.read_from_text(file);

    let mut new_poly = Poly {
        family_num: -3,
        // Set number of nodes. Needed for rotations.
        number_of_nodes: n_poly_nodes as isize,
        vertices: get_poly_coords(file, n_poly_nodes),
        ..Default::default()
    };

    // Get a normal vector
    // Vector from fist node to node across middle of polygon
    let mut pt_idx_12 = 3 * (n_poly_nodes / 2);

    if n_poly_nodes == 3 {
        pt_idx_12 = 8;
    }

    let p1 = Point3::from_slice(&new_poly.vertices[0..3]);
    let p2 = Point3::from_slice(&new_poly.vertices[3..6]);
    let p_12 = Point3::from_slice(&new_poly.vertices[pt_idx_12..pt_idx_12 + 3]);

    let v1 = p_12 - p1;
    let v2 = p2 - p1;

    new_poly.normal = v2.cross(&v1).normalize();
    // Estimate radius
    // across middle if even number of nodes
    // across middle close to perpendicular to xradius magnitude calculation
    new_poly.xradius = 0.5 * v2.magnitude();

    // Get idx for node 1/4 around polygon
    let pt_idx_14 = 3 * (n_poly_nodes / 4);
    // Get idx for node 3/4 around polygon
    let pt_idx_34 = 3 * (3 * n_poly_nodes / 4);

    let p_14 = Point3::from_slice(&new_poly.vertices[pt_idx_14..pt_idx_14 + 3]);
    let p_34 = Point3::from_slice(&new_poly.vertices[pt_idx_34..pt_idx_34 + 3]);

    new_poly.yradius = 0.5 * distance(&p_14, &p_34);
    new_poly.aspect_ratio = new_poly.yradius / new_poly.xradius;

    // Estimate translation (middle of poly)
    // Use midpoint between 1st and and half way around polygon
    // Note: For polygons defined by coordinates, the coordinates
    // themselves provide the translation. We need to estimate the center
    // of the polygon and init. the translation array
    new_poly.translation = 0.5 * Vector3::new(p1.x + p_12.x, p1.y + p_12.y, p1.z + p_12.z);

    new_poly
}

// /****************************************************************/
// /***********  Insert User Polygon By Coord  ********************/
// /*! Inserts user polygon using defined coordinates
//     provided by the user (see input file).
//     Intersection checking, FRAM, and rejection/acceptance are all contained
//     within this function.
//     Arg 1: Array for all accepted polygons
//     Arg 2: Array for all accepted intersections
//     Arg 3: Program statistics structure
//     Arg 4: Array of all triple intersection points */
#[allow(clippy::too_many_arguments)]
pub fn insert_user_polygon_by_coord(
    h: f64,
    eps: f64,
    r_fram: bool,
    disable_fram: bool,
    triple_intersections: bool,
    domain_size: &Vector3<f64>,
    polygon_file: &str,
    accepted_poly: &mut Vec<Poly>,
    intpts: &mut Vec<IntersectionPoints>,
    pstats: &mut Stats,
    triple_points: &mut Vec<Point3<f64>>,
) {
    let family_id = -3;

    println!("Reading User Defined Polygons from {}", &polygon_file);
    let mut file = File::open(polygon_file).unwrap();

    search_var(&mut file, "nPolygons:");
    let mut npoly: usize = 0;
    npoly.read_from_text(&mut file);

    println!("There are {} polygons", npoly);

    for i in 0..npoly {
        let mut new_poly = create_poly(&mut file);

        if domain_truncation(h, eps, &mut new_poly, domain_size) {
            // Poly completely outside domain
            pstats.rejection_reasons.outside += 1;
            pstats.rejected_poly_count += 1;
            println!("User Polygon (defined by coordinates) {} was rejected for being outside the defined domain.", i + 1);
            pstats
                .rejected_user_fracture
                .push(RejectedUserFracture::new(i + 1, family_id));
            continue; // Go to next poly (go to next iteration of for loop)
        }

        create_bounding_box(&mut new_poly);
        // Line of intersection and FRAM
        let reject_code = intersection_checking(
            h,
            eps,
            r_fram,
            disable_fram,
            triple_intersections,
            &mut new_poly,
            accepted_poly,
            intpts,
            pstats,
            triple_points,
        );

        if reject_code == 0 {
            // Incriment counter of accepted polys
            pstats.accepted_poly_count += 1;
            // Calculate poly's area
            new_poly.area = get_area(&new_poly);
            // Add new rejectsPerAttempt counter
            pstats.rejects_per_attempt.push(0);
            println!(
                "User Defined Polygon Fracture (Defined By Coordinates) {} Accepted",
                i + 1
            );
            accepted_poly.push(new_poly); // Save newPoly to accepted polys list
        } else {
            pstats.rejects_per_attempt[pstats.accepted_poly_count] += 1;
            pstats.rejected_poly_count += 1;
            println!(
                "Rejected User Defined Polygon Fracture (Defined By Coordinates) {}",
                i + 1
            );
            print_reject_reason(reject_code, &new_poly);
            pstats
                .rejected_user_fracture
                .push(RejectedUserFracture::new(i + 1, family_id));
        }
    }
}
