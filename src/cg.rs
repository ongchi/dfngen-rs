use parry3d_f64::bounding_volume::BoundingVolume;
use parry3d_f64::na::{distance, Point3, Rotation3, Unit, Vector3};
use parry3d_f64::query::details::closest_points_segment_segment_with_locations_nD;
use parry3d_f64::query::PointQuery;
use parry3d_f64::shape::Segment;

use crate::{
    fracture::{
        cluster_groups::{assign_group, update_groups},
        poly::Poly,
        user_defined_fractures::USER_FRACTURE_FAMILY,
    },
    math_functions::{max_elmt_idx, sorted_index, sum_deviation},
    structures::{IntersectionPoints, PolyOptions, Stats, TriplePtTempData},
};

// Check if two vectors are parallel
pub fn is_parallel(v1: &Vector3<f64>, v2: &Vector3<f64>, eps: f64) -> bool {
    let dot_prod = v1.dot(v2);
    1. - eps < dot_prod && dot_prod < 1. + eps
}

/// Upward-pointing normal vector of a plane given by the trend and plunge of its pole.
///
/// Angles are in radians. Trend is clockwise from North (+y), plunge is positive downward.
/// x is East, y is North, z is up.
pub fn normal_from_trend_plunge(trend: f64, plunge: f64) -> Vector3<f64> {
    Vector3::new(
        -trend.sin() * plunge.cos(),
        -trend.cos() * plunge.cos(),
        plunge.sin(),
    )
}

/// Normal vector of a plane given by its dip and strike (right-hand rule).
///
/// Angles are in radians. Strike is clockwise from North (+y).
/// x is East, y is North, z is up.
pub fn normal_from_dip_strike(dip: f64, strike: f64) -> Vector3<f64> {
    Vector3::new(
        dip.sin() * strike.cos(),
        -dip.sin() * strike.sin(),
        dip.cos(),
    )
}

/// Rotation which turns the unit vector `from` to the unit vector `to`.
pub fn rotation_between(from: &Vector3<f64>, to: &Vector3<f64>) -> Rotation3<f64> {
    Rotation3::rotation_between(from, to).unwrap_or_else(|| {
        // Opposite vectors: half turn about an axis perpendicular to them
        let not_parallel = if from.x.abs() < 0.9 {
            Vector3::x()
        } else {
            Vector3::y()
        };
        let axis = Unit::new_normalize(from.cross(&not_parallel));
        Rotation3::from_axis_angle(&axis, std::f64::consts::PI)
    })
}

/// Check Bounding Box***************************
///
/// Compares two polygons' bounding boxes, returns 1 if bounding boxes intersect
///
/// # Arguments
///
/// * `poly1` - Poly 1
/// * `poly2` - Poly 2
///
/// # Returns
///
/// * `bool` - True if bounding boxes intersect, false otherwise
fn check_bounding_box(poly1: &Poly, poly2: &Poly) -> bool {
    poly1.bounding_box.intersects(&poly2.bounding_box)
}

/// Find Intersections
/// Finds intersection end points of two intersecting polygons (Poly 1 and Poly 2)
/// Or, finds that polygons do not intersect (flag will = 0 )
///
/// # Arguments
///
/// * `flag` - The only flag which is currently used is '0'
///   0 - no intersection,
///   1 - intersection is completely inside poly 1,
///   2 - intersection is completely inside poly 2,
///   3 - intersection on both polys edges
/// * `poly1` - Poly 1
/// * `poly2` - Poly 2
///
/// Return: Intersection end points, Valid only if flag != 0 */
fn find_intersections(flag: &mut i32, poly1: &Poly, poly2: &Poly, eps: f64) -> IntersectionPoints {
    // This code is mostly converted directly from the mathematica version.
    // Re-write may be worth doing for increased performance and code clarity
    *flag = 0;

    let mut int_pts = IntersectionPoints::new();
    // The 2 points where the boundary of poly2 crosses the plane of poly1,
    // then the 2 points where the boundary of poly1 crosses the plane of poly2
    let mut inters = [Point3::origin(); 4];

    for (jj, (f1, f2)) in [(poly1, poly2), (poly2, poly1)].into_iter().enumerate() {
        let mut crossings = [Point3::origin(); 2];
        let mut count = 0;
        let plane_point = f1.vertices[0];
        let mut previous = &f2.vertices[f2.vertices.len() - 1];
        // Signed distance from the plane of F1 to the previous vertex of F2
        let mut prevdist = (previous - plane_point).dot(&f1.normal);

        for vertex in &f2.vertices {
            let currdist = (vertex - plane_point).dot(&f1.normal);

            if prevdist.abs() < eps {
                // Previous vertex is on the plane
                crossings[count] = *previous;
                count += 1;
            } else {
                let mut curr_times_prev = currdist * prevdist;

                if curr_times_prev.abs() < eps {
                    curr_times_prev = 0.;
                }

                if curr_times_prev < 0. {
                    // Consecutive vertices of F2 are at opposide sides of the plane,
                    // compute the point where the edge crosses the plane
                    let c = prevdist.abs() / (currdist.abs() + prevdist.abs());
                    crossings[count] = previous + (vertex - previous) * c;
                    count += 1;
                }
            }

            prevdist = currdist;
            previous = vertex;

            if count == 2 {
                break;
            }
        }

        match count {
            // No intersection
            0 => return int_pts,
            // Happens only when a single vertex of F2 is on the plane
            1 => crossings[1] = crossings[0],
            _ => {}
        }

        for crossing in &mut crossings {
            for coordinate in crossing.iter_mut() {
                if coordinate.abs() < eps {
                    *coordinate = 0.;
                }
            }
        }

        inters[2 * jj] = crossings[0];
        inters[2 * jj + 1] = crossings[1];
    }

    // Sort the 4 points along the axis on which they are spread the most
    let coordinates = |axis: usize| inters.map(|point| point[axis]);
    let stdev = [0, 1, 2].map(|axis| sum_deviation(&coordinates(axis)));
    let s = sorted_index(&coordinates(max_elmt_idx(&stdev)));

    // Points 0 and 1 are on the boundary of poly2, points 2 and 3 on the boundary of poly1.
    // If the smallest two points are on the boundary of the same poly, there is no intersection.
    if s[0] + s[1] == 1 || s[0] + s[1] == 5 {
        return int_pts;
    }

    // The polygons intersect. The middle two points form the intersection.
    int_pts.p1 = inters[s[1]];
    int_pts.p2 = inters[s[2]];

    *flag = match s[1] + s[2] {
        1 => 1, // Intersection inside poly1
        5 => 2, // Intersection inside poly2
        _ => 3, // Intersection on edges
    };

    int_pts
}

/// FRAM - Feature Rejection Algorithm for Meshing
///
/// Checks new poly and new intersection against other intersecting polygons
/// for violation of the minimum feature size
/// In some cases, the intersection may be shortened in order to accepted the fracture.
///
/// # Arguments
///
/// * `min_feature_size` - Minimum feature size
/// * `eps` - Epsilon value for floating point comparisons
/// * `r_fram` - Uses a relaxed version of the FRAM algorithm. The mesh may not be perfectly conforming
/// * `is_fram_enabled` - If false, FRAM is disabled
/// * `triple_intersections` - If true, triple intersections are accepted
/// * `int_pts` - Newest intersection found on new poly
/// * `count` - counter of number of intersections on new poly
/// * `int_pts_list` - Intersection endpoints list for entire DFN
/// * `new_poly` - New poly, poly being checked with FRAM
/// * `poly2` - Poly which new poly intersects with
/// * `pstats` - Stats structure (Should only be one, singleton)
/// * `temp_data` - Temp triple point data. Must keep intersections and triple points as temp data untill newPoly has been accepted
/// * `triple_points` - Triple points for entire DFN
/// * `temp_int_pts` - Temp intersection points. Must keep intersections and triple points as temp
///   data untill newPoly has been accepted
/// * `accepted_poly` - Array of all accepted polygons
///
/// # Returns
///
/// * `i32` - 0 if accepted, 1 if rejected
#[allow(clippy::too_many_arguments)]
fn fram(
    min_feature_size: f64,
    eps: f64,
    r_fram: bool,
    triple_intersections: bool,
    int_pts: &mut IntersectionPoints,
    count: usize,
    int_pts_list: &[IntersectionPoints],
    new_poly: &Poly,
    poly2: &Poly,
    pstats: &mut Stats,
    temp_data: &mut Vec<TriplePtTempData>,
    triple_points: &[Point3<f64>],
    temp_int_pts: &[IntersectionPoints],
    accepted_poly: &[Poly],
) -> i32 {
    /******* Check for intersection of length less than min_feature_size *******/
    if (int_pts.p1 - int_pts.p2).magnitude() < min_feature_size {
        //std::cout<<"\nrejectCode = -2: Intersection of length <= min_feature_size.\n";
        pstats.rejection_reasons.short_intersection += 1;
        return -2;
    }

    if !r_fram {
        /******************* distance to edges *****************/
        // Reject if intersection shirnks < 'shrinkLimit'
        let shrink_limit = 0.9 * (int_pts.p1 - int_pts.p2).magnitude();

        if check_close_edge(
            new_poly,
            int_pts,
            shrink_limit,
            pstats,
            min_feature_size,
            eps,
        ) {
            // std::cout<<"\nrejectCode = -6: Fracture too close to another fracture's edge.\n";
            pstats.rejection_reasons.close_to_edge += 1;
            return -6;
        }

        if check_close_edge(poly2, int_pts, shrink_limit, pstats, min_feature_size, eps) {
            // std::cout<<"\nrejectCode = -6: Fracture too close to another fracture's edge.\n";
            pstats.rejection_reasons.close_to_edge += 1;
            return -6;
        }

        /******* Intersection to Intersection Distance Checks ********/
        // Check distance from new intersection to other intersections on
        // poly2 (fracture newPoly is intersecting with)

        if check_dist_to_old_intersections(int_pts_list, int_pts, poly2, min_feature_size, eps) {
            pstats.rejection_reasons.inter_close_to_inter += 1;
            return -5;
        }

        // Check distance from new intersection to intersections already
        // existing on newPoly
        // Also checks for undetected triple points
        if check_dist_to_new_intersections(
            temp_int_pts,
            int_pts,
            temp_data,
            poly2,
            accepted_poly,
            min_feature_size,
            eps,
        ) {
            pstats.rejection_reasons.inter_close_to_inter += 1;
            return -5;
        }
    }

    /*************** Triple Intersection Checks *************/
    // NOTE: for debugging, there are several rejection codes for triple intersections
    // -14 <= rejCode <= -10 are for triple intersection rejections
    let rej_code = check_for_triple_intersections(
        min_feature_size,
        eps,
        triple_intersections,
        int_pts,
        count,
        int_pts_list,
        poly2,
        temp_data,
        triple_points,
    );

    if rej_code != 0 && !r_fram {
        pstats.rejection_reasons.triple += 1;
        return rej_code;
    }

    // Check if polys intersect on same plane
    if ((new_poly.normal[0] - poly2.normal[0]).abs()) < eps // If the normals are the same
                && (new_poly.normal[1] - poly2.normal[1]).abs() < eps
                && (new_poly.normal[2] - poly2.normal[2]).abs() < eps
    {
        // The intersection has already been found so we know that if the
        // normals are the same they must be on the same plane
        return -7;
    }

    0
}

/// FRAM CHECK
///
/// New intersction to Old Intersections Check
///
/// Checks distance of new intersection to intersections on poly2
///
/// # Arguments
///
/// * `int_pts_list` - Intersections arry for entire DFN
/// * `int_pts` - Current intersection being checked, intersection between newPoly and poly2
/// * `poly2` - Poly2
/// * `min_distance` - Minimum distance allowed if not a triple intersection
/// * `eps` - Epsilon value for floating point comparisons
///
/// # Returns
///
/// * `bool`
///   False - if not all distances are larger than minDistance or minDistance = 0 with triple intersection point
///   True - otherwise
fn check_dist_to_old_intersections(
    int_pts_list: &[IntersectionPoints],
    int_pts: &IntersectionPoints,
    poly2: &Poly,
    min_distance: f64,
    eps: f64,
) -> bool {
    let intersection = [int_pts.p1, int_pts.p2];

    for i in 0..poly2.intersection_index.len() {
        let int2 = [
            int_pts_list[poly2.intersection_index[i]].p1,
            int_pts_list[poly2.intersection_index[i]].p2,
        ];

        if are_segments_farther_than(&intersection, &int2, min_distance) {
            continue;
        }

        let (dist, _) = line_seg_to_line_seg(&intersection, &int2);

        if dist < (min_distance - eps) && dist > eps {
            return true;
        }
    }

    false
}

/// FRAM Check
///
/// New intersction to New Intersections Check
///
/// Checks distance of new intersection to other intersections on newPoly
///
/// # Arguments
///
/// * `temp_int_pts` - Array of intersections previously found on newPoly
/// * `int_pts` - Current intersection being checked, intersection between newPoly and poly2
/// * `temp_tri_pts` - Temp triple points, triple points found on newPoly
/// * `poly2` - Poly2
/// * `accepted_poly` - Array of all accepted polygons
/// * `min_distance` - Minimum distance allowed if not a triple intersection
/// * `eps` - Epsilon value for floating point comparisons
///
/// # Returns
///
/// * `bool`
///   False - if not all distances are larger than minDistance or minDistance = 0 with triple intersection point
///   True - otherwise
///
/// NOTE: If the distance between two intersections is 0, this function will verify the that
///       the triple interersection exists in 'tempTriPts'. If not found, fracture will be rejected
///       Due to the shrinkIntersection algorithm, it may be possible for a triple intersection point
///       to exist on only one fracture. This check resolves this issue.
fn check_dist_to_new_intersections(
    temp_int_pts: &[IntersectionPoints],
    int_pts: &IntersectionPoints,
    temp_tri_pts: &[TriplePtTempData],
    poly2: &Poly,
    accepted_poly: &[Poly],
    min_distance: f64,
    eps: f64,
) -> bool {
    let intersection = [int_pts.p1, int_pts.p2];

    for tmp_int in temp_int_pts {
        // `fract1` of an intersection on newPoly is the index of the other polygon.
        // The intersections with two pieces of one user-defined fracture can
        // touch, at the shared edge of the pieces.
        if accepted_poly[tmp_int.fract1 as usize].is_piece_of_same_fracture(poly2) {
            continue;
        }

        let int2 = [tmp_int.p1, tmp_int.p2];

        if are_segments_farther_than(&intersection, &int2, min_distance) {
            continue;
        }

        let (dist, pt) = line_seg_to_line_seg(&intersection, &int2);

        if dist < min_distance && dist > eps {
            return true;
        } else if dist < eps {
            // Make sure there is a triple intersection point
            // before accepting
            let mut reject = true;

            for tmp_tri in temp_tri_pts {
                // If the triple point is found, continue with checks, else reject
                if (pt.x - tmp_tri.triple_point.x).abs() < eps
                    && (pt.y - tmp_tri.triple_point.y).abs() < eps
                    && (pt.z - tmp_tri.triple_point.z).abs() < eps
                {
                    reject = false;
                    break;
                }
            }

            if reject {
                return true;
            }
        }
    }

    false
}

/// Shrink Intersection
/// Shrinks intersection untill the intersection is greater than 'minDist
/// to 'edge', or intersection shrinks to length < 'shrinkLimit'
///
/// 'firstNodeMinDist' can be used to allow a shoter first discretized node
/// distance. This allows for slight angles for intersections starting on the
/// edges of polygons without there the intersection being shortened.
/// If the first node is of distnace smaller than 'firstNodeMinDist', the
/// 'minDist' will be used from this point on to shorten the intersection.
///
/// # Arguments
///
/// * `int_pts` - Intersection being shrunk
/// * `edge` - Double array[6] of two end points which intersection is being tested against
/// * `shrink_limit` - Minimum length intersection is allowed to shrink
/// * `first_node_min_dist` - Fist node minimum distance
/// * `global_min_dist` - Minimum allowed distance between intersection and edge
///
/// # Returns
///
/// * `bool`
///   False - If intersection successfully shortened and minDist <= dist to edge && shrinkLimit <= intersection length
///   True - If intersection length shrinks to less than shrinkLimit
fn shrink_intersection(
    int_pts: &mut IntersectionPoints,
    edge: &[Point3<f64>; 2],
    shrink_limit: f64,
    first_node_min_dist: f64,
    global_min_dist: f64,
    eps: f64,
) -> bool {
    let vect = int_pts.p2 - int_pts.p1;
    let dist = vect.magnitude();
    // n is number of discrete points on intersection
    let n = (2. * dist / global_min_dist).ceil();
    let step_size = 1. / n;
    let pt = int_pts.p1;
    // Start step at first descrete point
    let mut step;

    // Check both sides of intersection against edge
    for i in 0..2 {
        let mut node_count = 0;

        if i == 0 {
            step = 0.;
        } else {
            step = 1.;
        }

        let point = pt + vect * step;
        let first_pt_dist_to_edge = point_to_line_seg(&point, edge, eps);
        let mut first_pt = true;

        while node_count <= n as i32 {
            node_count += 1;

            if i == 0 {
                step += step_size;
            } else {
                step -= step_size;
            }

            let pt_on_intersection = pt + vect * step;
            let dist = point_to_line_seg(&pt_on_intersection, edge, eps);

            if first_pt && (dist > first_node_min_dist) {
                // || (stepSize == 1 && dist < eps)))  {
                if first_pt_dist_to_edge < eps || first_pt_dist_to_edge >= global_min_dist {
                    // Leave intersection end point un-modified
                    break;
                }
            }

            first_pt = false;

            if dist > global_min_dist {
                if i == 0 {
                    int_pts.p1 = pt_on_intersection;
                } else {
                    int_pts.p2 = pt_on_intersection;
                }

                int_pts.intersection_shortened = true;
                break;
            }

            if node_count >= n as i32 {
                // All nodes are bad, reject
                return true;
            }
        }
    }

    // If intersection length shrank to less than shrinkLimit, reject
    (int_pts.p2 - int_pts.p1).magnitude() < shrink_limit
}

/// INTERSECTION CHECKING
///
/// This function will check for intersections with all polys whos bounding boxes intersect. It also will
/// run FRAM on the intersections.
/// 'newPoly' will be checked with FRAM one intersection at a time. At the first FRAM rejection, further
/// intersection checking will be aborted and the poly will be retranslated or thrown away.
/// This function saves intersections, and updates cluster groups when a poly is accepted.
/// This functino returns 0 if the poly was accepted and 1 if rejected. User needs to push the newPoly
/// into the accepted poly array if this function returns 0
///
/// # Arguments
///
/// * `new_poly` - Polygon being tested (newest poly to come into the DFN)
/// * `accepted_poly` - Array of all accepted polygons
/// * `int_pts_list` - Array of all accepted intersections
/// * `pstats` - Program statistics structure
/// * `triple_points` - Array of all accepted triple intersection points
///
/// # Returns
///
/// * `i32`
///   0 - Fracture had no intersections or features violating the minimum feature size (Passed all FRAM tests)
///   1 - Otherwise
pub fn intersection_checking(
    opts: &PolyOptions,
    new_poly: &mut Poly,
    accepted_poly: &mut [Poly],
    int_pts_list: &mut Vec<IntersectionPoints>,
    pstats: &mut Stats,
    triple_points: &mut Vec<Point3<f64>>,
) -> i32 {
    // List of fractures which new fracture intersected.
    // Used to update fractures intersections and
    // intersection count if newPoly is accepted
    let mut temp_intersect_list = Vec::new();
    let mut temp_int_pts = Vec::new();
    let mut temp_original_intersection = Vec::new();
    // Index to newPoly's position if accepted
    let new_poly_index = accepted_poly.len();
    // Index to intpts if newPoly intersections
    let int_pts_index = int_pts_list.len();
    let mut encountered_groups = Vec::new();
    // Counts number of accepted intersections on newPoly.
    let mut count = 0;
    let mut temp_data = Vec::new();

    for (ii, poly) in accepted_poly.iter().enumerate() {
        let mut flag = 0;
        let mut intersection;

        if new_poly.is_piece_of_same_fracture(poly) {
            // Pieces of one user-defined fracture are always in the same group.
            // Their contact is not an intersection, and FRAM does not check it.
            join_group(new_poly, poly, &mut encountered_groups);
            continue;
        }

        // NOTE: findIntersections() searches bounding boxes
        // Bounding box search
        if check_bounding_box(new_poly, poly) {
            intersection = find_intersections(&mut flag, new_poly, poly, opts.eps);

            if flag != 0 {
                // If flag != 0, intersection exists
                // Holds origintal intersection, used to update
                // stats on how much intersections were shortened
                temp_original_intersection.push(intersection.clone());
                // FRAM returns 0 if no intersection problems.
                // 'count' is number of already accepted intersections on new poly
                let is_user_fracture_pair = new_poly.family_num == USER_FRACTURE_FAMILY
                    && poly.family_num == USER_FRACTURE_FAMILY;

                // The option can turn off FRAM for two user-defined fractures
                let reject_code = if !opts.is_fram_enabled
                    || (is_user_fracture_pair && !opts.check_user_fractures)
                {
                    0
                } else {
                    fram(
                        opts.min_feature_size,
                        opts.eps,
                        opts.r_fram,
                        opts.triple_intersections,
                        &mut intersection,
                        count,
                        int_pts_list,
                        new_poly,
                        poly,
                        pstats,
                        &mut temp_data,
                        triple_points,
                        &temp_int_pts,
                        accepted_poly,
                    )
                };

                // If intersection is NOT rejected
                if reject_code == 0 {
                    // If FRAM returned 0, everything is OK
                    // Update group numbers, intersection indexes
                    intersection.fract1 = ii as isize; // ii is the intersecting fracture
                    intersection.fract2 = new_poly_index as isize; // newPolys ID/index in array of accpted polys, if accepted
                                                                   // 'tempIntersectList' keeps all fracture #'s intersecting with newPoly
                    temp_intersect_list.push(ii); // Save fracture index to update if newPoly accepted
                                                  // Save index to polys intersections (intPts) array
                    new_poly.intersection_index.push(int_pts_index + count);
                    count += 1; // Increment counter of number of fractures intersecting with newPoly
                    temp_int_pts.push(intersection); // Save intersection

                    join_group(new_poly, poly, &mut encountered_groups);
                } else {
                    // 'newPoly' rejected
                    // SAVE REJECTED POLYS HERE IF THIS FUNCTINALITY IS NEEDED
                    return reject_code; // Break loop/function and return 1, poly is rejected
                }
            }
        }
    } // If it makes it here, no problematic intersections with new polygon. SAVE POLY AND UPDATE GROUP/CLUSTER INFO

    // Done searching intersections. All FRAM tests have passed. Polygon is accetepd.

    // After searching for intersections, if newPoly still has no intersection or group:
    if new_poly.group_num == 0 {
        // 'newPoly' had no intersections. Assign it to its own/new group number
        assign_group(new_poly, pstats, new_poly_index);
    } else {
        // Intersections exist and were accepted, newPoly already has group number
        // Save temp. intersections to intPts (permanent array),
        // Update all intersected polygons intersection-points index lists and intersection count
        // Update polygon group lists
        // Append temp intersection points array to permanent one
        int_pts_list.extend(temp_int_pts.clone());

        // Update poly's indexs to intersections list (intpts)
        for i in 0..temp_intersect_list.len() {
            // Update each intersected poly's intersection index
            accepted_poly[temp_intersect_list[i]]
                .intersection_index
                .push(int_pts_index + i);
            // intPtsIndex+i will be the index position of the intersection once it is saved to the intersections array
        }

        // Fracture is now accepted.
        // Update intersection structures with triple intersection points.
        // triple intersection points will be found 3 times (3 fractures make up one triple int point)
        // We need only to save the point to the permanent triplePoints array once, and then give each intersection a
        // reference to it.
        if opts.triple_intersections {
            let trip_index = triple_points.len();

            for (j, tri_pt_tmp) in temp_data.iter().enumerate() {
                // Loop through newly found triple points
                triple_points.push(tri_pt_tmp.triple_point);

                // Update index pointers to the triple intersection points
                for ii in 0..tri_pt_tmp.int_index.len() {
                    let idx = tri_pt_tmp.int_index[ii];
                    int_pts_list[idx].triple_points_idx.push(trip_index + j);
                }
            }
        }

        // ***********************
        // Update group numbers **
        // ***********************
        update_groups(
            new_poly,
            accepted_poly,
            &encountered_groups,
            pstats,
            new_poly_index,
        );
    }

    // Keep track of how much intersection length we are looising from 'shrinkIntersection()'
    // Calculate and store total original intersection length (all intersections)
    // and actual intersection length, after intersection has been shortened.
    for i in 0..temp_int_pts.len() {
        let length =
            (temp_original_intersection[i].p1 - temp_original_intersection[i].p2).magnitude();
        pstats.original_length += length;

        if temp_int_pts[i].intersection_shortened {
            pstats.intersections_shortened += 1;
            let new_length = (temp_int_pts[i].p1 - temp_int_pts[i].p2).magnitude();
            pstats.discarded_length += length - new_length;
        }
    }

    0
}

/// Puts `new_poly` into the group of `poly`, a fracture which it is connected to.
///
/// If `new_poly` has a group already and the group of `poly` is another one,
/// the group of `poly` goes into `encountered_groups`, for the merge of the
/// groups when `new_poly` is accepted.
fn join_group(new_poly: &mut Poly, poly: &Poly, encountered_groups: &mut Vec<usize>) {
    if new_poly.group_num == 0 {
        new_poly.group_num = poly.group_num;
    } else if new_poly.group_num != poly.group_num {
        encountered_groups.push(poly.group_num);
    }
}

/// Disatance from intersection line to Nodes/Vertices
///
/// Checks distance from line of intersection to poly vertices
///
/// # Arguments
///
/// * `poly` - Poly to check dist to intersection
/// * `int_pts` - Intersection structure (intersection end points)
/// * `min_dist` - Minimum distance allowed
/// * `pstats` - Program statistics structure
/// * `eps` - Epsilon value for floating point comparisons
///
/// # Returns
///
/// * `bool`
///   False - Distance from intersection to vertices are all > min_dist
///   True - No distance less than minDist was found
fn check_distance_from_nodes(
    poly: &Poly,
    int_pts: &IntersectionPoints,
    min_dist: f64,
    pstats: &mut Stats,
    eps: f64,
) -> bool {
    // Intersection dist to poly vertices
    let line = [int_pts.p1, int_pts.p2];
    let mut dist = point_to_line_seg(&poly.vertices[0], &line, eps);

    for vertex in &poly.vertices[1..] {
        let temp = point_to_line_seg(vertex, &line, eps);

        // If new distance is less than the last calculated distance...
        if temp < dist {
            dist = temp;
        }

        if dist < min_dist && dist > eps {
            pstats.rejection_reasons.close_to_node += 1;
            return true;
        }
    }

    // If it makes it here, no distance less than 'minDist' found
    false
}

/// Shortest Disatance, point to line seg
///
/// # Arguments
///
/// * `point` - Point in 3d space
/// * `line` - Line defined by two end points
///
/// # Returns
///
/// * `f64` - Shortest distance between point and line segment
fn point_to_line_seg(point: &Point3<f64>, line: &[Point3<f64>; 2], eps: f64) -> f64 {
    if (line[0] - line[1]).magnitude_squared() < eps {
        // Line endpoints are equal to each other
        return distance(point, &line[0]);
    }

    Segment::new(line[0], line[1]).distance_to_local_point(point, true)
}

/// Check if nodes are too close to edge
///
/// Checks distances from intersection to poly edges. If the distance is less than
/// min_feature_size, the intersection is allowed to shrink by %10 of its original length. If
/// the intersection is still closer than min_feature_size to a poly edge, the polygon is rejected.
///
/// # Arguments
///
/// * `poly1` - Poly to be tested
/// * `int_pts` - Intersection points to be tested
/// * `shrink_limit` - Minimum length the intersection is allowed to shrink to
/// * `pstats` - Program statistics structure, used to report stats on how much
///   intersection length is being reduced by from shrinkIntersection()
/// * `min_dist` - Minimum distance allowed from an end point to the edge of a polygon
///   if the intersection does not land accross a poly's edge
/// * `eps` - Epsilon value for floating point comparisons
fn check_close_edge(
    poly1: &Poly,
    int_pts: &mut IntersectionPoints,
    shrink_limit: f64,
    pstats: &mut Stats,
    min_dist: f64,
    eps: f64,
) -> bool {
    // 'line' is newest intersection end points
    let line = [int_pts.p1, int_pts.p2];
    // Counts how many endPoints are on the polys edge.
    // If both end points of intersection are on polys edge,
    // we must check the distance from end points to
    // vertices
    let mut on_edge_count = 0;

    let number_of_nodes = poly1.vertices.len();

    for i in 0..number_of_nodes {
        // Edge from the vertex to the next one. The last edge closes the polygon.
        let edge = [poly1.vertices[i], poly1.vertices[(i + 1) % number_of_nodes]];

        if are_segments_farther_than(&edge, &line, min_dist) {
            // All distances between this edge and the intersection are larger than min_dist
            continue;
        }

        let mut end_pts_to_edge = [
            point_to_line_seg(&line[0], &edge, eps),
            point_to_line_seg(&line[1], &edge, eps),
            point_to_line_seg(&edge[0], &line, eps),
            point_to_line_seg(&edge[1], &line, eps),
        ];
        // Sort smallest to largest
        end_pts_to_edge.sort_by(|a, b| a.partial_cmp(b).unwrap());

        // If two smallest distances are < min_feature_size,
        // the line is almost parallel and closer to edge than min_feature_size, reject it
        if (end_pts_to_edge[0] < min_dist && end_pts_to_edge[1] < min_dist)
            && end_pts_to_edge[0] > eps
        {
            return true;
        }

        // Minimum dist from poly edge segment to intersection segment
        let (dist, _) = line_seg_to_line_seg(&edge, &line);

        if dist < min_dist && dist > eps {
            // Try to shrink the intersection slightly in order to
            // not reject the polygon
            if shrink_intersection(int_pts, &edge, shrink_limit, min_dist, min_dist, eps) {
                // Returns one if insterection shrinks to less than .9*min_feature_size
                return true;
            }
        } else if dist < eps {
            // Endpoint is almost exactly on poly's edge, must check
            // whether the discretized nodes will ne closer
            // than the minimum allowed distance
            // shrinkIntersection() will check if the descretized intersection points violate any
            // distance to edge rules
            // NOTE: Intersections discretize with set size = .5*min_feature_size
            // Minimum distance to edge must be less than .5*min_feature_size to allow for angles
            let min_dist2 = 0.4 * min_dist;

            if shrink_intersection(int_pts, &edge, shrink_limit, min_dist2, min_dist, eps) {
                //returns one if insterection shrinks to less than .9*min_feature_size
                return true;
            }

            on_edge_count += 1;

            // If both end points are on the edge we must also check the distnace
            // of the intersection to the poly vertices.
            // IF the intersecion is within the polygon, distances less than min_feature_size to
            // edges and vertices will be caught by lineSegToLineSeg() in this function.
            if (on_edge_count >= 2)
                & check_distance_from_nodes(poly1, int_pts, min_dist, pstats, eps)
            {
                return true;
            }
        }
    }

    #[cfg(feature = "disable_shortening_int")]
    if int_pts.intersection_shortened {
        return true;
    }

    false
}

/// Fast check for two line segments which are far from each other
///
/// The distance between the bounding boxes of two segments is not larger than
/// the distance between the segments. Thus this check finds most of the far
/// pairs without the calculation of the exact distance.
///
/// # Arguments
///
/// * `line1` - Line 1 defined by two end points
/// * `line2` - Line 2 defined by two end points
/// * `limit` - Distance limit
///
/// # Returns
///
/// true - The distance between the segments is larger than `limit`
/// false - The distance can be larger or smaller than `limit`
fn are_segments_farther_than(
    line1: &[Point3<f64>; 2],
    line2: &[Point3<f64>; 2],
    limit: f64,
) -> bool {
    // The margin makes the check safe against rounding of the exact distance
    const RELATIVE_MARGIN: f64 = 1e-6;

    let mut squared_box_distance = 0.;

    for axis in 0..3 {
        let (min1, max1) = (
            line1[0][axis].min(line1[1][axis]),
            line1[0][axis].max(line1[1][axis]),
        );
        let (min2, max2) = (
            line2[0][axis].min(line2[1][axis]),
            line2[0][axis].max(line2[1][axis]),
        );
        let gap = (min1 - max2).max(min2 - max1).max(0.);
        squared_box_distance += gap * gap;
    }

    let limit_with_margin = limit * (1. + RELATIVE_MARGIN);
    squared_box_distance > limit_with_margin * limit_with_margin
}

/// Closest Distance from Line Seg to Line Seg
///
/// # Arguments
///
/// * `line1` - Line 1 defined by two end points
/// * `line2` - Line 2 defined by two end points
///
/// # Returns
///
/// The minimum distance between line 1 and line 2, and the point on line 1
/// which is closest to line 2. If the lines intersect, this is the intersection point.
fn line_seg_to_line_seg(line1: &[Point3<f64>; 2], line2: &[Point3<f64>; 2]) -> (f64, Point3<f64>) {
    let (location1, location2) = closest_points_segment_segment_with_locations_nD(
        (&line1[0], &line1[1]),
        (&line2[0], &line2[1]),
    );
    let point1 = Segment::new(line1[0], line1[1]).point_at(&location1);
    let point2 = Segment::new(line2[0], line2[1]).point_at(&location2);

    (distance(&point1, &point2), point1)
}

/// Check for Triple Intersection , get int. point
///
/// Check for triple intersection features of less than min_feature_size.
/// Returns rejection code if fracture is rejected,  zero if accepted
///
/// Rejection codes:
///     0 = poly accepted
///     -10 = triple_intersectionsNotAllowed (rejected triple intersections
///           due to triple intersections not allowed in input file)
///     -11 = triple_closeToIntersection (newPoly's intersection landed too close to a previous intersection)
///     -12 = triple_smallIntersectionAngle
///     -13 = triple_closeEndPoint   (triple intersection point too close to an endpoint)
///     -14 = triple_closeToTriplePt  (new triple point too close to previous triple point) */
#[allow(clippy::too_many_arguments)]
fn check_for_triple_intersections(
    min_feature_size: f64,
    eps: f64,
    triple_intersections: bool,
    int_pts: &IntersectionPoints,
    count: usize,
    int_pts_list: &[IntersectionPoints],
    poly2: &Poly,
    temp_data: &mut Vec<TriplePtTempData>,
    triple_points: &[Point3<f64>],
) -> i32 {
    let min_dist = 1.5 * min_feature_size;
    let int_end_pts = [int_pts.p1, int_pts.p2]; //newest intersection
                                                // Number of intersections already on poly2
    let n = poly2.intersection_index.len();

    // Check new fracure's new intesrsection against previous intersections on poly2
    for i in 0..n {
        let intersection_index = poly2.intersection_index[i]; // Index to previous intersecction (old intersection)
        let intersection_index2 = int_pts_list.len() + count; // Index to current intersection, if accepted (new intersection)
        let line = [
            int_pts_list[intersection_index].p1,
            int_pts_list[intersection_index].p2,
        ];

        if are_segments_farther_than(&int_end_pts, &line, min_feature_size) {
            continue;
        }

        let (mut dist1, pt) = line_seg_to_line_seg(&int_end_pts, &line); // distance and pt of intersection

        if dist1 >= min_feature_size {
            continue;
        }

        if !triple_intersections && dist1 < min_feature_size {
            return -10; // Triple intersections not allowed (user input option)
        }

        if dist1 > eps && dist1 < min_feature_size {
            return -11; // Point too close to other intersection
        }

        // Overlaping intersections
        if dist1 <= eps {
            // ANGLE CHECK, using definition of dot product: A dot B = Mag(A)*Mag(B) * Cos(angle)
            // Normalize  A and B first. Compare to precalculated values of cos(47deg)= .681998 and cos(133deg)= -.681998
            let u = (line[1] - line[0]).normalize();
            let v = (int_pts.p2 - int_pts.p1).normalize();
            let dot_prod = u.dot(&v);

            if !(-0.68199836..=0.68199836).contains(&dot_prod) {
                //if angle is less than 47 deg or greater than 133, reject
                return -12;
            }

            // Check that the triple intersection point isn't too close to an endpoint
            dist1 = distance(&pt, &int_end_pts[0]);
            let mut dist2 = distance(&pt, &int_end_pts[1]);

            if dist1 < min_feature_size || dist2 < min_feature_size {
                return -13;
            }

            dist1 = distance(&pt, &line[0]);
            dist2 = distance(&pt, &line[1]);

            if dist1 < min_feature_size || dist2 < min_feature_size {
                return -13;
            }

            // Check that the triple intersection point isn't too close to previous triple intersection points
            for k in 0..int_pts_list[intersection_index].triple_points_idx.len() {
                let idx = int_pts_list[intersection_index].triple_points_idx[k];
                dist1 = distance(&pt, &triple_points[idx]);

                if dist1 < min_dist {
                    return -14;
                }
            }

            // Look at previously found triple points on current line of intersection,
            // do not save duplicate triple ponits
            // If duplicate point is found, a different intersection may still be needing an update
            // (three intersections per single triple intersection point )
            let mut duplicate = false;
            let mut k = 0;

            // See if the found triple pt is already saved
            while k < temp_data.len() {
                if (pt.x - temp_data[k].triple_point.x).abs() < eps
                    && (pt.y - temp_data[k].triple_point.y).abs() < eps
                    && (pt.z - temp_data[k].triple_point.z).abs() < eps
                {
                    duplicate = true;
                    break;
                }
                k += 1;
            }

            // If it is a duplicate triple pt
            if duplicate {
                // The old fracture's intersection needs
                // a reference to the new triple int pt
                // Save the index to the old fractures intersection
                // which needs updating. If newPoly is accpted
                // the old fracture will be updated
                temp_data[k].int_index.push(intersection_index2);
            } else {
                // Not a duplicate point
                // Save temp triple point data, store to permanent at end if fracture is accepted
                temp_data.push(TriplePtTempData {
                    triple_point: pt,
                    int_index: vec![intersection_index, intersection_index2],
                });
            }
        }
    }

    // Test distances to triple points on own line of intersection
    // If any triple points are found to be less than 'minDist' to each other, reject.
    let trip_size = temp_data.len();
    // Each tempData structure contains 1 triple pt and the intersection index of intersection it is on

    if trip_size != 0 {
        let y = trip_size - 1;

        for k in 0..y {
            let point1 = &temp_data[k].triple_point;

            for tri_pt_tmp in temp_data.iter().take(trip_size).skip(k + 1) {
                let point2 = &tri_pt_tmp.triple_point;
                let dist = distance(point1, point2);

                if dist < min_dist && dist > eps {
                    return -14;
                }
            }
        }
    }

    0
}

#[cfg(test)]
mod tests {
    use std::f64::consts::FRAC_PI_2;

    use super::*;

    fn assert_vector_eq(actual: Vector3<f64>, expected: [f64; 3]) {
        assert!(
            (actual - Vector3::from(expected)).magnitude() < 1e-12,
            "expected {expected:?}, found {:?}",
            actual.as_slice()
        );
    }

    /// A segment of length 2 along x, and a segment along y at the given
    /// distance from it in z
    fn crossed_segments(z_distance: f64) -> ([Point3<f64>; 2], [Point3<f64>; 2]) {
        (
            [Point3::new(-1., 0., 0.), Point3::new(1., 0., 0.)],
            [
                Point3::new(0., -1., z_distance),
                Point3::new(0., 1., z_distance),
            ],
        )
    }

    #[test]
    fn test_segments_with_a_box_distance_above_the_limit_are_far() {
        let (line1, line2) = crossed_segments(2.);

        assert!(are_segments_farther_than(&line1, &line2, 1.));
    }

    #[test]
    fn test_segments_with_a_distance_below_the_limit_are_not_far() {
        let (line1, line2) = crossed_segments(0.5);

        assert!(!are_segments_farther_than(&line1, &line2, 1.));
    }

    #[test]
    fn test_segments_with_a_distance_equal_to_the_limit_are_not_far() {
        let (line1, line2) = crossed_segments(1.);

        assert!(!are_segments_farther_than(&line1, &line2, 1.));
    }

    #[test]
    fn test_far_segments_with_touching_boxes_are_not_found_as_far() {
        // The distance is sqrt(2), but the bounding boxes have the corner (1, 1, 0) in common
        let line1 = [Point3::new(0., 0., 0.), Point3::new(1., 1., 0.)];
        let line2 = [Point3::new(1., 3., 0.), Point3::new(3., 1., 0.)];

        assert!(!are_segments_farther_than(&line1, &line2, 1.));
    }

    #[test]
    fn test_far_check_agrees_with_the_exact_distance() {
        let (line1, line2) = crossed_segments(3.);
        let (exact_distance, _) = line_seg_to_line_seg(&line1, &line2);

        assert!(are_segments_farther_than(&line1, &line2, 2.) && exact_distance > 2.);
    }

    #[test]
    fn test_vertical_plunge_gives_upward_normal() {
        assert_vector_eq(normal_from_trend_plunge(0.3, FRAC_PI_2), [0., 0., 1.]);
    }

    #[test]
    fn test_horizontal_pole_with_trend_north_gives_normal_to_south() {
        assert_vector_eq(normal_from_trend_plunge(0., 0.), [0., -1., 0.]);
    }

    #[test]
    fn test_horizontal_pole_with_trend_east_gives_normal_to_west() {
        assert_vector_eq(normal_from_trend_plunge(FRAC_PI_2, 0.), [-1., 0., 0.]);
    }

    #[test]
    fn test_zero_dip_gives_upward_normal() {
        assert_vector_eq(normal_from_dip_strike(0., 1.1), [0., 0., 1.]);
    }

    #[test]
    fn test_vertical_plane_with_strike_north_gives_normal_to_east() {
        assert_vector_eq(normal_from_dip_strike(FRAC_PI_2, 0.), [1., 0., 0.]);
    }

    #[test]
    fn test_vertical_plane_with_strike_east_gives_normal_to_south() {
        assert_vector_eq(normal_from_dip_strike(FRAC_PI_2, FRAC_PI_2), [0., -1., 0.]);
    }

    #[test]
    fn test_normals_from_angles_have_unit_length() {
        assert!((normal_from_trend_plunge(0.7, 0.4).magnitude() - 1.).abs() < 1e-12);
        assert!((normal_from_dip_strike(0.7, 0.4).magnitude() - 1.).abs() < 1e-12);
    }
}
