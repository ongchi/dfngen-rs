use std::fs::File;

use crate::{
    computational_geometry::{create_bounding_box, intersection_checking},
    domain::domain_truncation,
    insert_shape::print_reject_reason,
    math_functions::get_area,
    read_input::Input,
    read_input_functions::{read_var, search_var},
    structures::{IntPoints, Point, Poly, RejectedUserFracture, Stats},
    vector_functions::{cross_product, euclidean_distance, magnitude, normalize},
};

// **********************************************************************
// **********************************************************************
// Used to read in ellipse coordinates when the user is using
// user ellipses defined by coordinates option.
// Arg 1: ifstream file object
// Arg 2: OUTPUT. Pointer to array to store the coordinates
// Arg 3: Number of ellipses
// Arg 4: Number of points per ellipse
fn get_poly_coords(stream: &mut File, out_ary: &mut [f64], n_vertices: usize) {
    for i in 0..n_vertices {
        let x = i * 3;
        out_ary[x] = read_var(stream);
        out_ary[x + 1] = read_var(stream);
        out_ary[x + 2] = read_var(stream);
    }
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
pub fn insert_user_polygon_by_coord(
    input: &Input,
    accepted_poly: &mut Vec<Poly>,
    intpts: &mut Vec<IntPoints>,
    pstats: &mut Stats,
    triple_points: &mut Vec<Point>,
) {
    let mut n_poly_nodes: usize;
    println!(
        "Domain Size {} {} {}",
        input.domainSize[0], input.domainSize[1], input.domainSize[2]
    );
    println!("Reading User Defined Polygons from {}", input.polygonFile);
    let mut file = File::open(&input.polygonFile).unwrap();
    search_var(&mut file, "nPolygons:");
    let n_polygon_by_coord: usize = read_var(&mut file);
    println!("There are {} polygons", n_polygon_by_coord);
    accepted_poly.reserve(n_polygon_by_coord);

    for i in 0..n_polygon_by_coord {
        let mut new_poly = Poly::default();
        let mut rejected_user_fracture = RejectedUserFracture::default();
        new_poly.family_num = -3;
        n_poly_nodes = read_var(&mut file);
        new_poly.number_of_nodes = n_poly_nodes as isize;
        new_poly.vertices.reserve(3 * n_poly_nodes); // 3 * number of nodes
        get_poly_coords(&mut file, &mut new_poly.vertices, n_poly_nodes);

        // Get a normal vector
        // Vector from fist node to node across middle of polygon
        let mut mid_pt_idx = 3 * (n_poly_nodes / 2);

        if n_poly_nodes == 3 {
            mid_pt_idx = 8;
        }

        let v1 = [
            new_poly.vertices[mid_pt_idx] - new_poly.vertices[0],
            new_poly.vertices[mid_pt_idx + 1] - new_poly.vertices[1],
            new_poly.vertices[mid_pt_idx + 2] - new_poly.vertices[2],
        ];
        // Vector from first node to 2nd node
        let v2 = [
            new_poly.vertices[3] - new_poly.vertices[0],
            new_poly.vertices[4] - new_poly.vertices[1],
            new_poly.vertices[5] - new_poly.vertices[2],
        ];
        let x_prod1 = cross_product(&v2, &v1);
        // Set normal vector
        new_poly.normal[0] = x_prod1[0]; //x
        new_poly.normal[1] = x_prod1[1]; //y
        new_poly.normal[2] = x_prod1[2]; //z
        normalize(&mut new_poly.normal);
        // Estimate radius
        new_poly.xradius = 0.5 * magnitude(v2[0], v2[1], v2[2]); // across middle if even number of nodes
                                                                 // across middle close to perpendicular to xradius magnitude calculation
        let temp_idx1 = 3 * (n_poly_nodes / 4); // Get idx for node 1/4 around polygon
        let temp_idx2 = 3 * (3 * n_poly_nodes / 4); // Get idx for node 1/4 around polygon
        new_poly.yradius = 0.5
            * euclidean_distance(
                &new_poly.vertices[temp_idx1..temp_idx1 + 3]
                    .try_into()
                    .unwrap(),
                &new_poly.vertices[temp_idx2..temp_idx2 + 3]
                    .try_into()
                    .unwrap(),
            );
        new_poly.aspect_ratio = new_poly.yradius / new_poly.xradius;
        // Estimate translation (middle of poly)
        // Use midpoint between 1st and and half way around polygon
        // Note: For polygons defined by coordinates, the coordinates
        // themselves provide the translation. We need to estimate the center
        // of the polygon and init. the translation array
        new_poly.translation[0] = 0.5 * (new_poly.vertices[0] + new_poly.vertices[mid_pt_idx]);
        new_poly.translation[1] = 0.5 * (new_poly.vertices[1] + new_poly.vertices[mid_pt_idx + 1]);
        new_poly.translation[2] = 0.5 * (new_poly.vertices[2] + new_poly.vertices[mid_pt_idx + 2]);

        if domain_truncation(input, &mut new_poly, &input.domainSize) {
            // Poly completely outside domain
            pstats.rejection_reasons.outside += 1;
            pstats.rejected_poly_count += 1;
            println!("User Polygon (defined by coordinates) {} was rejected for being outside the defined domain.", i + 1);
            let mut idx;

            for j in 0..n_poly_nodes {
                idx = j * 3;
                println!("{:?}", &new_poly.vertices[idx..idx + 3]);
            }

            rejected_user_fracture.id = (i + 1) as isize;
            rejected_user_fracture.user_fracture_type = -3;
            pstats.rejected_user_fracture.push(rejected_user_fracture);
            new_poly.vertices.clear();
            continue; // Go to next poly (go to next iteration of for loop)
        }

        create_bounding_box(&mut new_poly);
        // Line of intersection and FRAM
        let reject_code = intersection_checking(
            input,
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
                                          //std::cout << "size of accepted Poly " << acceptedPoly.size() << std::endl;
        } else {
            new_poly.vertices.clear(); // Need to delete manually, created with new[]
            pstats.rejects_per_attempt[pstats.accepted_poly_count] += 1;
            pstats.rejected_poly_count += 1;
            println!(
                "Rejected User Defined Polygon Fracture (Defined By Coordinates) {}",
                i + 1
            );
            print_reject_reason(reject_code, &new_poly);
            rejected_user_fracture.id = (i + 1) as isize;
            rejected_user_fracture.user_fracture_type = -3;
            pstats.rejected_user_fracture.push(rejected_user_fracture);
        }
    }
}
