use parry3d_f64::na::{distance, Point3, Vector3};
use tracing::info;

use super::domain::domain_truncation;
use super::insert_shape::print_reject_reason;

use crate::{
    computational_geometry::{create_bounding_box, intersection_checking},
    math_functions::get_area,
    structures::{IntersectionPoints, Poly, RejectedUserFracture, Stats},
};

fn create_poly(n_ell_nodes: usize, user_ell_coord_vertices: &[f64], idx: usize) -> Poly {
    let mut new_poly = Poly {
        family_num: -1,
        // Set number of nodes. Needed for rotations.
        number_of_nodes: n_ell_nodes as isize,
        // Initialize normal to {0,0,1}. need initialized for 3D rotation
        normal: Vector3::new(0., 0., 1.),
        ..Default::default()
    };

    new_poly.vertices.reserve(n_ell_nodes * 3);

    let poly_vert_idx = idx * 3 * n_ell_nodes; // Each polygon has nEllNodes * 3 vertices

    // Initialize vertices
    for j in 0..n_ell_nodes {
        let v_idx = j * 3;
        new_poly.vertices[v_idx] = user_ell_coord_vertices[poly_vert_idx + v_idx];
        new_poly.vertices[v_idx + 1] = user_ell_coord_vertices[poly_vert_idx + 1 + v_idx];
        new_poly.vertices[v_idx + 2] = user_ell_coord_vertices[poly_vert_idx + 2 + v_idx];
    }

    // Get a normal vector
    // Vector from fist node to node accross middle of polygon
    let pt_idx_12 = 3 * (n_ell_nodes / 2);
    let p1 = Point3::from_slice(&new_poly.vertices[0..3]);
    let p2 = Point3::from_slice(&new_poly.vertices[3..6]);
    let p_12 = Point3::from_slice(&new_poly.vertices[pt_idx_12..pt_idx_12 + 3]);

    let v1 = p_12 - p1;
    let v2 = p2 - p1;

    new_poly.normal = v2.cross(&v1).normalize();
    // Estimate radius
    // across middle if even number of nodes
    new_poly.xradius = 0.5 * v2.magnitude();

    // Get idx for node 1/4 around polygon
    let pt_idx_14 = 3 * (pt_idx_12 / 2);
    // Get idx for node 3/4 around polygon
    // across middle close to perpendicular to xradius magnitude calculation
    let pt_idx_34 = 3 * (pt_idx_14 + pt_idx_12);

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

/// Insert User Ellipses By Coord
///
/// Inserts user ellipses using defined coordinates
/// provided by the user (see input file).
/// Intersection checking, FRAM, and rejection/accptance are all contained
/// within this function.
///
/// # Arguments
///
/// * `h` - Minimum feature size
/// * `eps` - Epsilon value for floating point comparisons
/// * `n_ell_nodes` - Number of nodes for user defined ellipses by coordinates
/// * `n_ell_by_coord` - Number of user ellipses defined by coordinates
/// * `r_fram` - Uses a relaxed version of the FRAM algorithm. The mesh may not be perfectly conforming
/// * `disable_fram` - If true, FRAM is disabled
/// * `triple_intersections` - If true, triple intersections are accepted
/// * `domain_size` - Size of the domain
/// * `user_ell_coord_vertices` - Array of ellipse coordiates.
/// * `accepted_poly` - Array of all accepted polygons
/// * `intpts` - Array of all accepted intersections
/// * `pstats` - Program statistics structure
/// * `triple_points` - Array of all triple intersection points
#[allow(clippy::too_many_arguments)]
pub fn insert_user_ell_by_coord(
    h: f64,
    eps: f64,
    n_ell_nodes: usize,
    n_ell_by_coord: usize,
    r_fram: bool,
    disable_fram: bool,
    triple_intersections: bool,
    domain_size: &Vector3<f64>,
    user_ell_coord_vertices: &mut Vec<f64>,
    accepted_poly: &mut Vec<Poly>,
    intpts: &mut Vec<IntersectionPoints>,
    pstats: &mut Stats,
    triple_points: &mut Vec<Point3<f64>>,
) {
    let family_id = -1;
    let npoly = n_ell_by_coord;
    info!("{} User Ellipses By Coordinates Defined", npoly);

    for i in 0..npoly {
        let mut new_poly = create_poly(n_ell_nodes, user_ell_coord_vertices, i);

        if domain_truncation(h, eps, &mut new_poly, domain_size) {
            // Poly completely outside domain
            pstats.rejection_reasons.outside += 1;
            pstats.rejected_poly_count += 1;
            info!("User Ellipse (defined by coordinates) {} was rejected for being outside the defined domain.", i + 1);
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
            // If intersection is ok (FRAM passed all tests)
            if new_poly.truncated {
                pstats.truncated += 1;
            }

            // Incriment counter of accepted polys
            pstats.accepted_poly_count += 1;
            // Calculate poly's area
            new_poly.area = get_area(&new_poly);
            // Add new rejectsPerAttempt counter
            pstats.rejects_per_attempt.push(0);
            info!(
                "User Defined Elliptical Fracture (Defined By Coordinates) {} Acceptedn",
                i + 1
            );
            accepted_poly.push(new_poly); // Save newPoly to accepted polys list
        } else {
            pstats.rejects_per_attempt[pstats.accepted_poly_count] += 1;
            pstats.rejected_poly_count += 1;
            info!(
                "Rejected Eser Defined Elliptical Fracture (Defined By Coordinates) {}",
                i + 1
            );
            print_reject_reason(reject_code, &new_poly);
            pstats
                .rejected_user_fracture
                .push(RejectedUserFracture::new(i + 1, family_id));
        }
    } // End loop

    user_ell_coord_vertices.clear();
}
