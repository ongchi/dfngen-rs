use parry3d_f64::na::{Point3, Vector3};

use super::domain::domain_truncation;
use super::insert_shape::print_reject_reason;

use crate::{
    computational_geometry::{create_bounding_box, intersection_checking},
    math_functions::get_area,
    structures::{IntersectionPoints, Poly, RejectedUserFracture, Stats},
};

fn create_poly(user_rect_coord_vertices: &[f64], idx: usize) -> Poly {
    let mut new_poly = Poly {
        family_num: -2,
        // Set number of nodes. Needed for rotations.
        number_of_nodes: 4,
        ..Default::default()
    };

    new_poly.vertices.reserve(12); // 4 * {x,y,z}
    let poly_vert_idx = idx * 12; // Each polygon has 4 vertices (12 elements, 4*{x,y,z}))

    // Initialize vertices
    for j in 0..4 {
        let v_idx = j * 3;
        new_poly.vertices[v_idx] = user_rect_coord_vertices[poly_vert_idx + v_idx];
        new_poly.vertices[v_idx + 1] = user_rect_coord_vertices[poly_vert_idx + 1 + v_idx];
        new_poly.vertices[v_idx + 2] = user_rect_coord_vertices[poly_vert_idx + 2 + v_idx];
    }

    // Check that rectangle lays one a single plane:
    // let xProd1 = cross Product vector (1st node to 2nd node) with vector(1st node to 3rd node)
    // and xProd2 = cross product vector (1st node to 3th node) with vector (1st node to 4th node)
    // Then, cross product xProd1 and xProd2, if this produces zero vector, all coords are on the same plane
    // v1 is vector from first vertice to third vertice
    // Vector from fist node to 3rd node (vector through middle of sqare)
    let p1 = Point3::from_slice(&new_poly.vertices[0..3]);
    let p2 = Point3::from_slice(&new_poly.vertices[3..6]);
    let p3 = Point3::from_slice(&new_poly.vertices[6..9]);
    let p4 = Point3::from_slice(&new_poly.vertices[9..12]);

    let v1 = p3 - p1;
    let v2 = p2 - p1;
    let v3 = p4 - p1;
    let x_prod1 = v2.cross(&v1).normalize();
    // let x_prod2 = crossProduct(&v3, &v1);
    // let x_prod3 = crossProduct(&x_prod1, &x_prod2);
    //will be zero vector if all vertices are on the same plane
    //TODO: Error check below is too sensitive. Adjust it.
    // Error check for points not on the same plane
    //        if (std::abs(magnitude(xProd3[0],xProd3[1],xProd3[2])) > eps) { //points do not lay on the same plane. reject poly else meshing will fail
    // if (!(std::abs(xProd3[0]) < eps && std::abs(xProd3[1]) < eps && std::abs(xProd3[2]) < eps)) {
    //     delete[] newPoly.vertices;
    //     pstats.rejectedPolyCount++;
    //     std::cout << "\nUser Rectangle (defined by coordinates) " << i+1 << " was rejected. The defined vertices are not co-planar.\n";
    //     std::cout << "Please check user defined coordinates for rectanle " << i+1 << " in input file\n";
    //     delete[] xProd1;
    //     delete[] xProd2;
    //     delete[] xProd3;
    //     continue; //go to next poly
    // }

    // Set normal vector
    new_poly.normal = x_prod1;

    // Set radius (x and y radii might be switched based on order of users coordinates)
    new_poly.xradius = 0.5 * v2.magnitude();
    new_poly.yradius = 0.5 * v3.magnitude();
    new_poly.aspect_ratio = new_poly.yradius / new_poly.xradius;

    // Estimate translation
    // Use midpoint between 1st and 3rd vertices
    // Note: For polygons defined by coordinates, the coordinates
    // themselves provide the translation. We are just filling the
    // translation array for completeness even though the translation
    // array might not be used
    new_poly.translation = 0.5 * Vector3::new(p1.x + p3.x, p1.y + p3.y, p1.z + p3.z);

    new_poly
}

/// Insert User Rects By Coord
///
/// Inserts user rectangles using defined coordinates
/// provided by the user (see input file).
/// Intersection checking, FRAM, and rejection/accptance are all contained
/// within this function.
///
/// # Arguments
///
/// * `h` - Minimum feature size
/// * `eps` - Epsilon value for floating point comparisons
/// * `r_fram` - Uses a relaxed version of the FRAM algorithm. The mesh may not be perfectly conforming
/// * `disable_fram` - If true, FRAM is disabled
/// * `triple_intersections` - If true, triple intersections are accepted
/// * `n_rect_by_coord` - Number of user rectangles defined by coordinates
/// * `domain_size` - Size of the domain
/// * `user_rect_coord_vertices` - User defined rectangle coordinates
/// * `accepted_poly` - Array for all accepted polygons
/// * `intpts` - Array for all accepted intersections
/// * `pstats` - Program statistics structure
/// * `triple_points` - Array of all triple intersection points
#[allow(clippy::too_many_arguments)]
pub fn insert_user_rects_by_coord(
    h: f64,
    eps: f64,
    r_fram: bool,
    disable_fram: bool,
    triple_intersections: bool,
    n_rect_by_coord: usize,
    domain_size: &Vector3<f64>,
    user_rect_coord_vertices: &mut Vec<f64>,
    accepted_poly: &mut Vec<Poly>,
    intpts: &mut Vec<IntersectionPoints>,
    pstats: &mut Stats,
    triple_points: &mut Vec<Point3<f64>>,
) {
    let family_id = -2;
    let npoly = n_rect_by_coord;
    println!("{} User Rectangles By Coordinates Defined", npoly);

    for i in 0..npoly {
        let mut new_poly = create_poly(user_rect_coord_vertices, i);

        if domain_truncation(h, eps, &mut new_poly, domain_size) {
            // Poly completely outside domain
            pstats.rejection_reasons.outside += 1;
            pstats.rejected_poly_count += 1;
            println!("User Rectangle (defined by coordinates) {} was rejected for being outside the defined domain.", i + 1);
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
            println!(
                "User Defined Rectangular Fracture (Defined By Coordinates) {} Accepted",
                i + 1
            );
            accepted_poly.push(new_poly); // Save newPoly to accepted polys list
        } else {
            pstats.rejects_per_attempt[pstats.accepted_poly_count] += 1;
            pstats.rejected_poly_count += 1;
            println!(
                "Rejected User Defined Rectangular Fracture (Defined By Coordinates) {}",
                i + 1
            );
            print_reject_reason(reject_code, &new_poly);
            pstats
                .rejected_user_fracture
                .push(RejectedUserFracture::new(i + 1, family_id));
        }
    }

    user_rect_coord_vertices.clear();
}
