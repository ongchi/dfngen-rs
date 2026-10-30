use parry3d_f64::bounding_volume::BoundingVolume;
use parry3d_f64::na::{distance, Point3};
use serde::{Deserialize, Serialize};
use tracing::debug;

use crate::{
    dfn::{Dfn, Intersection, PolyOptions},
    fracture::poly::{FractureOrigin, Poly},
    geometry::{are_segments_farther_than, is_same_point, line_seg_to_line_seg, point_to_line_seg},
    region_boundaries::{check_boundary_lines, is_intersection_near_boundary},
    stats::Stats,
};

/// Options of FRAM. The input has them only when FRAM is on.
#[derive(Debug, Clone, PartialEq)]
pub struct FramOptions {
    /// Relaxed FRAM: only short intersections and fractures in the same plane
    /// are rejected. The mesh may not be perfectly conforming.
    pub relaxed: bool,

    /// Accept (true) or reject (false) triple intersections
    pub triple_intersections: bool,

    /// False: FRAM does not check an intersection of two user-defined fractures
    pub check_user_fractures: bool,

    /// Keep the rejected fractures and write them to `fractures_rejected.vtp`.
    /// Caution: the file can be very large.
    pub write_rejected_fractures: bool,

    /// Number of FRAM rejects of a stochastic fracture before it is rejected.
    /// After each reject, the fracture moves to a new position.
    pub rejects_per_fracture: usize,
}

/// Reason why FRAM rejected a fracture
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum RejectReason {
    /// Intersection shorter than `min_feature_size`
    ShortIntersection,
    /// Intersection too near to another intersection
    CloseToIntersection,
    /// Fracture too near to an edge of another fracture
    CloseToEdge,
    /// Intersecting fractures are in the same plane
    SamePlane,
    /// Triple intersection, and the input does not permit triple intersections
    TripleIntersectionNotAllowed,
    /// Intersection too near to a previous intersection
    CloseToPreviousIntersection,
    /// Triple intersection with an angle that is too small
    TripleAngleTooSmall,
    /// Triple intersection point too near to an end point of an intersection
    TriplePointCloseToEndpoint,
    /// Triple intersection point too near to another triple intersection point
    TriplePointCloseToTriplePoint,
    /// Boundary line of the region shorter than `min_feature_size`
    ShortBoundaryLine,
    /// Boundary line of the region too near to the edge of its fracture
    BoundaryCloseToEdge,
    /// Intersection too near to a boundary line of the region, or crossing it badly
    IntersectionCloseToBoundary,
}

impl RejectReason {
    /// Name of the reason in `fractures_rejected.vtp`: the name of the item
    pub fn name(self) -> &'static str {
        match self {
            Self::ShortIntersection => "ShortIntersection",
            Self::CloseToIntersection => "CloseToIntersection",
            Self::CloseToEdge => "CloseToEdge",
            Self::SamePlane => "SamePlane",
            Self::TripleIntersectionNotAllowed => "TripleIntersectionNotAllowed",
            Self::CloseToPreviousIntersection => "CloseToPreviousIntersection",
            Self::TripleAngleTooSmall => "TripleAngleTooSmall",
            Self::TriplePointCloseToEndpoint => "TriplePointCloseToEndpoint",
            Self::TriplePointCloseToTriplePoint => "TriplePointCloseToTriplePoint",
            Self::ShortBoundaryLine => "ShortBoundaryLine",
            Self::BoundaryCloseToEdge => "BoundaryCloseToEdge",
            Self::IntersectionCloseToBoundary => "IntersectionCloseToBoundary",
        }
    }

    /// Text of the reason for the log
    pub fn description(self) -> &'static str {
        match self {
            Self::ShortIntersection => "Intersection of length < min_feature_size.",
            Self::CloseToIntersection => "Intersection too close to another intersection.",
            Self::CloseToEdge => "Fracture too close to another fracture's edge.",
            Self::SamePlane => "Fractures intersecting on same plane",
            Self::TripleIntersectionNotAllowed => "Rejected triple intersection due to triple intersections being turned off in input file.",
            Self::CloseToPreviousIntersection => "Fracture's intersection landed too close to a previous intersection.",
            Self::TripleAngleTooSmall => "Fracture created a triple intersection with an angle too small.",
            Self::TriplePointCloseToEndpoint => "Fracture created a triple intersection with the triple intersection point too close to an intersection's endpoint.",
            Self::TriplePointCloseToTriplePoint => "Fracture created a triple intersection with the triple intersection point too close to another triple intersection point.",
            Self::ShortBoundaryLine => "Boundary line of the region on the fracture is shorter than min_feature_size.",
            Self::BoundaryCloseToEdge => "Boundary line of the region is too close to the edge of the fracture.",
            Self::IntersectionCloseToBoundary => "Intersection is too close to a boundary line of the region, or crosses it at a small angle or near its end.",
        }
    }
}

/// Triple intersection point of a new fracture which FRAM has not accepted yet
///
/// It goes into the triple points of the DFN when the fracture is accepted, and
/// it is discarded when the fracture is rejected.
pub struct CandidateTriplePoint {
    /// Position of the triple intersection point
    pub triple_point: Point3<f64>,
    /// Indices of the intersections which the point is on. They get a reference to
    /// the point when the fracture is accepted.
    pub intersection_indices: Vec<usize>,
}

/// Writes the reject reason of an attempt of a stochastic fracture to the log,
/// at the level `debug`.
///
/// # Arguments
///
/// * `reject_reason` - Reason why FRAM rejected the fracture
/// * `family_index` - Index of the family of the fracture
pub fn log_reject_reason(reject_reason: RejectReason, family_index: usize) {
    debug!(
        "Fracture of family {} rejected: {} ({})",
        family_index + 1,
        reject_reason.name(),
        reject_reason.description()
    );
}

/// Sum of the squared deviations from the mean
///
/// It gives the same order of spread as the standard deviation, without the
/// division and the square root.
fn sum_of_squared_deviations(vec: &[f64]) -> f64 {
    let mean = vec.iter().sum::<f64>() / (vec.len() as f64);
    vec.iter().map(|&v| (v - mean).powi(2)).sum()
}

/// Index of the largest element. If elements are equal, the last one wins.
fn max_element_index(vec: &[f64]) -> usize {
    vec.iter()
        .enumerate()
        .reduce(|max, el| if el.1 < max.1 { max } else { el })
        .map(|(i, _)| i)
        .unwrap_or(0)
}

/// Indices of the elements in the order of their values, smallest first
fn sorted_index(vec: &[f64]) -> Vec<usize> {
    let mut v = vec.iter().enumerate().collect::<Vec<(usize, &f64)>>();
    // Not total_cmp: it puts -0.0 before 0.0, and that changes the order of equal values
    v.sort_by(|a, b| (a.1).partial_cmp(b.1).expect("the values are not NaN"));
    v.into_iter().map(|(i, _)| i).collect()
}

/// Finds the line of intersection of two polygons.
///
/// # Returns
///
/// The intersection with its end points, or `None` if the polygons do not intersect.
/// The fractures of the intersection are not set.
fn find_intersections(poly: &Poly, other_poly: &Poly, eps: f64) -> Option<Intersection> {
    // This code is mostly converted directly from the mathematica version.
    // Re-write may be worth doing for increased performance and code clarity
    // The 2 points where the boundary of other_poly crosses the plane of poly,
    // then the 2 points where the boundary of poly crosses the plane of other_poly
    let mut inters = [Point3::origin(); 4];

    for (jj, (f1, f2)) in [(poly, other_poly), (other_poly, poly)]
        .into_iter()
        .enumerate()
    {
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
            0 => return None,
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
    let stdev = [0, 1, 2].map(|axis| sum_of_squared_deviations(&coordinates(axis)));
    let s = sorted_index(&coordinates(max_element_index(&stdev)));

    // Points 0 and 1 are on the boundary of other_poly, points 2 and 3 on the boundary of poly.
    // If the smallest two points are on the boundary of the same poly, there is no intersection.
    if s[0] + s[1] == 1 || s[0] + s[1] == 5 {
        return None;
    }

    // The polygons intersect. The middle two points form the intersection.
    Some(Intersection {
        p1: inters[s[1]],
        p2: inters[s[2]],
        ..Default::default()
    })
}

/// The accepted fractures of the DFN, with their intersections and triple points
struct AcceptedDfn<'a> {
    polys: &'a [Poly],
    intersections: &'a [Intersection],
    triple_points: &'a [Point3<f64>],
}

/// Intersections and triple points of a new fracture which FRAM accepted
///
/// They are kept apart from the DFN until all intersections of the new
/// fracture are accepted. If one is rejected, they are discarded.
#[derive(Default)]
struct CandidateIntersections {
    intersections: Vec<Intersection>,
    triple_points: Vec<CandidateTriplePoint>,
}

/// FRAM - Feature Rejection Algorithm for Meshing
///
/// Checks a new intersection of a new fracture for features smaller than the
/// minimum feature size. The intersection can be shortened to accept the fracture.
///
/// # Arguments
///
/// * `options` - Limits and options of FRAM
/// * `accepted` - Accepted fractures of the DFN, with their intersections and triple points
/// * `candidate` - Intersections and triple points which FRAM accepted before for the new fracture
/// * `intersection` - New intersection between `new_poly` and `other_poly`
/// * `new_poly` - New fracture
/// * `other_poly` - Accepted fracture which the new fracture intersects
/// * `stats` - Statistics of the generation
///
/// # Returns
///
/// The reject reason, if the intersection is rejected
fn fram(
    options: &PolyOptions,
    accepted: &AcceptedDfn,
    candidate: &mut CandidateIntersections,
    intersection: &mut Intersection,
    new_poly: &Poly,
    other_poly: &Poly,
    stats: &mut Stats,
) -> Result<(), RejectReason> {
    let Some(fram_options) = &options.fram else {
        // FRAM is off
        return Ok(());
    };
    let min_feature_size = options.min_feature_size;
    let eps = options.eps;

    // Check for intersection of length less than min_feature_size
    if (intersection.p1 - intersection.p2).magnitude() < min_feature_size {
        stats.rejection_reasons.short_intersection += 1;
        return Err(RejectReason::ShortIntersection);
    }

    if !fram_options.relaxed {
        // distance to edges
        // Reject if intersection shirnks < `shrink_limit`
        let shrink_limit = 0.9 * (intersection.p1 - intersection.p2).magnitude();

        if check_close_edge(
            new_poly,
            intersection,
            shrink_limit,
            stats,
            min_feature_size,
            eps,
        ) {
            stats.rejection_reasons.close_to_edge += 1;
            return Err(RejectReason::CloseToEdge);
        }

        if check_close_edge(
            other_poly,
            intersection,
            shrink_limit,
            stats,
            min_feature_size,
            eps,
        ) {
            stats.rejection_reasons.close_to_edge += 1;
            return Err(RejectReason::CloseToEdge);
        }

        // Intersection to boundary line checks, after the shortening at the edges
        let segment = [intersection.p1, intersection.p2];
        if [new_poly, other_poly].iter().any(|poly| {
            is_intersection_near_boundary(&poly.boundary_lines, &segment, min_feature_size, eps)
        }) {
            stats.rejection_reasons.boundary_line += 1;
            return Err(RejectReason::IntersectionCloseToBoundary);
        }

        // Intersection to Intersection Distance Checks
        // Check distance from new intersection to other intersections on
        // other_poly (fracture the new fracture is intersecting with)

        if is_too_close_to_old_intersections(
            accepted.intersections,
            intersection,
            other_poly,
            min_feature_size,
            eps,
        ) {
            stats.rejection_reasons.inter_close_to_inter += 1;
            return Err(RejectReason::CloseToIntersection);
        }

        // Check distance from new intersection to intersections already
        // existing on the new fracture
        // Also checks for undetected triple points
        if is_too_close_to_new_intersections(
            &candidate.intersections,
            intersection,
            &candidate.triple_points,
            other_poly,
            accepted.polys,
            min_feature_size,
            eps,
        ) {
            stats.rejection_reasons.inter_close_to_inter += 1;
            return Err(RejectReason::CloseToIntersection);
        }
    }

    // Triple Intersection Checks
    // NOTE: for debugging, there are several reject reasons for triple intersections
    let triple_result = check_for_triple_intersections(
        options,
        fram_options,
        accepted,
        candidate,
        intersection,
        other_poly,
    );

    if triple_result.is_err() && !fram_options.relaxed {
        stats.rejection_reasons.triple += 1;
        return triple_result;
    }

    // Check if polys intersect on same plane
    if ((new_poly.normal[0] - other_poly.normal[0]).abs()) < eps // If the normals are the same
                && (new_poly.normal[1] - other_poly.normal[1]).abs() < eps
                && (new_poly.normal[2] - other_poly.normal[2]).abs() < eps
    {
        // The intersection has already been found so we know that if the
        // normals are the same they must be on the same plane
        return Err(RejectReason::SamePlane);
    }

    Ok(())
}

/// True if the new intersection is too near to an intersection which the other
/// fracture has already.
///
/// # Arguments
///
/// * `intersections` - Intersections of the DFN
/// * `intersection` - New intersection between the new fracture and `other_poly`
/// * `other_poly` - Accepted fracture which the new fracture intersects
/// * `min_distance` - Minimum distance between two intersections which do not touch
/// * `eps` - Epsilon value for floating point comparisons
fn is_too_close_to_old_intersections(
    intersections: &[Intersection],
    intersection: &Intersection,
    other_poly: &Poly,
    min_distance: f64,
    eps: f64,
) -> bool {
    let new_segment = [intersection.p1, intersection.p2];

    for &index in &other_poly.intersection_indices {
        let other_segment = [intersections[index].p1, intersections[index].p2];

        if are_segments_farther_than(&new_segment, &other_segment, min_distance) {
            continue;
        }

        let (dist, _) = line_seg_to_line_seg(&new_segment, &other_segment);

        if dist < (min_distance - eps) && dist > eps {
            return true;
        }
    }

    false
}

/// True if the new intersection is too near to another new intersection of the
/// new fracture.
///
/// Two intersections which touch must have a triple intersection point in
/// `candidate_triple_points`. If they have none, the result is true: a shortened
/// intersection can give a triple intersection point which only one fracture has.
///
/// # Arguments
///
/// * `candidate_intersections` - Intersections which FRAM accepted before for the new fracture
/// * `intersection` - New intersection between the new fracture and `other_poly`
/// * `candidate_triple_points` - Triple points which FRAM found before for the new fracture
/// * `other_poly` - Accepted fracture which the new fracture intersects
/// * `accepted_polys` - Accepted fractures of the DFN
/// * `min_distance` - Minimum distance between two intersections which do not touch
/// * `eps` - Epsilon value for floating point comparisons
fn is_too_close_to_new_intersections(
    temp_int_pts: &[Intersection],
    intersection: &Intersection,
    candidate_triple_points: &[CandidateTriplePoint],
    other_poly: &Poly,
    accepted_polys: &[Poly],
    min_distance: f64,
    eps: f64,
) -> bool {
    let new_segment = [intersection.p1, intersection.p2];

    for candidate_intersection in temp_int_pts {
        // `fracture1` of an intersection on the new fracture is the index of the other polygon.
        // The intersections with two pieces of one user-defined fracture can
        // touch, at the shared edge of the pieces.
        if accepted_polys[candidate_intersection.fracture1].is_piece_of_same_fracture(other_poly) {
            continue;
        }

        let other_segment = [candidate_intersection.p1, candidate_intersection.p2];

        if are_segments_farther_than(&new_segment, &other_segment, min_distance) {
            continue;
        }

        let (dist, pt) = line_seg_to_line_seg(&new_segment, &other_segment);

        if dist < min_distance && dist > eps {
            return true;
        } else if dist < eps {
            // Make sure there is a triple intersection point
            // before accepting
            let has_triple_point = candidate_triple_points
                .iter()
                .any(|candidate_triple_point| {
                    is_same_point(&pt, &candidate_triple_point.triple_point, eps)
                });

            if !has_triple_point {
                return true;
            }
        }
    }

    false
}

/// Shortens an intersection at its two ends until each end is farther than
/// `global_min_dist` from an edge.
///
/// The intersection is examined at discrete points. For the first point from an
/// end, `first_node_min_dist` is the limit: a smaller limit there permits an
/// intersection which starts on the edge at a small angle. After the first point,
/// `global_min_dist` is the limit.
///
/// # Arguments
///
/// * `intersection` - Intersection to shorten
/// * `edge` - The two end points of the edge
/// * `shrink_limit` - Minimum length of the shortened intersection
/// * `first_node_min_dist` - Minimum distance to the edge for the first point
/// * `global_min_dist` - Minimum distance between the intersection and the edge
/// * `eps` - Epsilon value for floating point comparisons
///
/// # Returns
///
/// True if the intersection must be rejected: all its points are too near to the
/// edge, or it is shorter than `shrink_limit`.
fn shrink_intersection(
    intersection: &mut Intersection,
    edge: &[Point3<f64>; 2],
    shrink_limit: f64,
    first_node_min_dist: f64,
    global_min_dist: f64,
    eps: f64,
) -> bool {
    let vect = intersection.p2 - intersection.p1;
    let dist = vect.magnitude();
    // n is number of discrete points on intersection
    let n = (2. * dist / global_min_dist).ceil();
    let step_size = 1. / n;
    let pt = intersection.p1;
    let node_limit = n as i32;
    let Intersection {
        p1,
        p2,
        is_shortened,
        ..
    } = intersection;

    // Check both sides of intersection against edge: from p1 in the direction
    // of p2, then from p2 in the direction of p1
    for (start, step_to_next, end_point) in [(0., step_size, p1), (1., -step_size, p2)] {
        let mut node_count = 0;
        // Start step at first descrete point
        let mut step = start;

        let point = pt + vect * step;
        let first_pt_dist_to_edge = point_to_line_seg(&point, edge, eps);
        let mut first_pt = true;

        while node_count <= node_limit {
            node_count += 1;
            step += step_to_next;

            let pt_on_intersection = pt + vect * step;
            let dist = point_to_line_seg(&pt_on_intersection, edge, eps);

            if first_pt
                && dist > first_node_min_dist
                && (first_pt_dist_to_edge < eps || first_pt_dist_to_edge >= global_min_dist)
            {
                // Leave intersection end point un-modified
                break;
            }

            first_pt = false;

            if dist > global_min_dist {
                *end_point = pt_on_intersection;
                *is_shortened = true;
                break;
            }

            if node_count >= node_limit {
                // All nodes are bad, reject
                return true;
            }
        }
    }

    // If intersection length shrank to less than `shrink_limit`, reject
    (intersection.p2 - intersection.p1).magnitude() < shrink_limit
}

/// Checks a new fracture against the accepted fractures.
///
/// The function finds the intersections with each fracture whose bounding box
/// intersects that of the new fracture, and runs FRAM on each intersection. It
/// stops at the first intersection which FRAM rejects. If all intersections are
/// accepted, it saves them in the DFN and updates the cluster groups and the
/// statistics. The caller then adds the new fracture to the accepted fractures.
///
/// # Arguments
///
/// * `options` - Limits and options of FRAM
/// * `new_poly` - New fracture
/// * `dfn` - The DFN: accepted fractures, intersections, triple points, cluster groups and statistics
///
/// # Returns
///
/// `Ok` if the fracture has no feature smaller than the minimum feature size, or
/// the reject reason of FRAM
pub fn check_intersections(
    options: &PolyOptions,
    new_poly: &mut Poly,
    dfn: &mut Dfn,
) -> Result<(), RejectReason> {
    let Dfn {
        stats,
        clusters,
        accepted_polys,
        intersections,
        triple_points,
        ..
    } = dfn;
    // List of fractures which new fracture intersected.
    // Used to update fractures intersections and
    // intersection count if the new fracture is accepted
    let mut intersected_poly_indices = Vec::new();
    let mut candidate = CandidateIntersections::default();
    // Length of each intersection before FRAM shortens it
    let mut original_lengths = Vec::new();
    // Index to the new fracture's position if accepted
    let new_poly_index = accepted_polys.len();
    // Index to intersections if the new fracture intersections
    let first_intersection_index = intersections.len();
    let mut encountered_groups = Vec::new();

    // Boundary lines of the region on the new fracture. A user-defined
    // fracture is checked only if FRAM checks the user-defined fractures.
    if let Some(fram_options) = &options.fram {
        let is_checked =
            new_poly.origin != FractureOrigin::UserDefined || fram_options.check_user_fractures;

        if is_checked {
            check_boundary_lines(
                &new_poly.boundary_lines,
                &new_poly.vertices,
                options.min_feature_size,
                fram_options.relaxed,
            )
            .inspect_err(|_| stats.rejection_reasons.boundary_line += 1)?;
        }
    }

    for (ii, poly) in accepted_polys.iter().enumerate() {
        if new_poly.is_piece_of_same_fracture(poly) {
            // Pieces of one user-defined fracture are always in the same group.
            // Their contact is not an intersection, and FRAM does not check it.
            join_group(new_poly, poly, &mut encountered_groups);
            continue;
        }

        // NOTE: `find_intersections` searches bounding boxes
        // Bounding box search
        if new_poly.bounding_box.intersects(&poly.bounding_box) {
            if let Some(mut intersection) = find_intersections(new_poly, poly, options.eps) {
                // Used to update the statistics on how much the intersections were shortened
                original_lengths.push((intersection.p1 - intersection.p2).magnitude());
                // FRAM returns Ok if no intersection problems.
                let is_user_fracture_pair = new_poly.origin == FractureOrigin::UserDefined
                    && poly.origin == FractureOrigin::UserDefined;

                // The option can turn off FRAM for two user-defined fractures
                let checks_user_fractures = options
                    .fram
                    .as_ref()
                    .is_some_and(|fram_options| fram_options.check_user_fractures);
                let fram_result = if is_user_fracture_pair && !checks_user_fractures {
                    Ok(())
                } else {
                    let accepted = AcceptedDfn {
                        polys: accepted_polys,
                        intersections,
                        triple_points,
                    };

                    fram(
                        options,
                        &accepted,
                        &mut candidate,
                        &mut intersection,
                        new_poly,
                        poly,
                        stats,
                    )
                };

                // If intersection is NOT rejected
                if fram_result.is_ok() {
                    // If FRAM returned Ok, everything is OK
                    // Update group numbers, intersection indexes
                    intersection.fracture1 = ii; // ii is the intersecting fracture
                    intersection.fracture2 = new_poly_index; // the new fracture's ID/index in array of accpted polys, if accepted
                                                             // `intersected_poly_indices` keeps all fracture #'s intersecting with the new fracture
                    intersected_poly_indices.push(ii); // Save fracture index to update if the new fracture accepted
                                                       // Save index to polys intersections intersections of the DFN
                    new_poly
                        .intersection_indices
                        .push(first_intersection_index + candidate.intersections.len());
                    candidate.intersections.push(intersection); // Save intersection

                    join_group(new_poly, poly, &mut encountered_groups);
                } else {
                    // the new fracture rejected
                    // SAVE REJECTED POLYS HERE IF THIS FUNCTINALITY IS NEEDED
                    return fram_result; // Break loop/function, poly is rejected
                }
            }
        }
    } // If it makes it here, no problematic intersections with new polygon. SAVE POLY AND UPDATE GROUP/CLUSTER INFO

    // Done searching intersections. All FRAM tests have passed. Polygon is accetepd.

    // After searching for intersections, if the new fracture still has no intersection or group:
    let Some(group) = new_poly.group else {
        // the new fracture had no intersections. Assign it to its own/new group number.
        // It has no intersection lengths for the statistics.
        new_poly.group = Some(clusters.add_isolated(new_poly_index));
        return Ok(());
    };

    // Intersections exist and were accepted, the new fracture already has group number
    // Save temp. intersections to intersections of the DFN (permanent array),
    // Update all intersected polygons intersection-points index lists and intersection count
    // Update polygon group lists

    // Keep track of how much intersection length we are looising from `shrink_intersection`'
    // Calculate and store total original intersection length (all intersections)
    // and actual intersection length, after intersection has been shortened.
    for (&length, shortened) in original_lengths.iter().zip(&candidate.intersections) {
        stats.original_length += length;

        if shortened.is_shortened {
            stats.intersections_shortened += 1;
            let new_length = (shortened.p1 - shortened.p2).magnitude();
            stats.discarded_length += length - new_length;
        }
    }

    // Append temp intersection points array to permanent one
    intersections.extend(candidate.intersections);

    // Update poly's indexs to intersections list (intersections)
    for (i, &poly_index) in intersected_poly_indices.iter().enumerate() {
        // first_intersection_index + i is the position of the intersection in the permanent array
        accepted_polys[poly_index]
            .intersection_indices
            .push(first_intersection_index + i);
    }

    // Fracture is now accepted.
    // Give the intersections their triple intersection points.
    // triple intersection points will be found 3 times (3 fractures make up one triple int point)
    // We need only to save the point to the permanent triple points of the DFN once, and then give each intersection a
    // reference to it.
    if options
        .fram
        .as_ref()
        .is_some_and(|fram_options| fram_options.triple_intersections)
    {
        let first_triple_point_index = triple_points.len();

        for (j, candidate_triple_point) in candidate.triple_points.iter().enumerate() {
            // Loop through newly found triple points
            triple_points.push(candidate_triple_point.triple_point);

            // Each intersection of the point gets the index of the point
            for &idx in &candidate_triple_point.intersection_indices {
                intersections[idx]
                    .triple_point_indices
                    .push(first_triple_point_index + j);
            }
        }
    }

    // Update group numbers **
    clusters.add_connected(group, new_poly_index, &encountered_groups, accepted_polys);

    Ok(())
}

/// Puts `new_poly` into the group of `poly`, a fracture which it is connected to.
///
/// If `new_poly` has a group already and the group of `poly` is another one,
/// the group of `poly` goes into `encountered_groups`, for the merge of the
/// groups when `new_poly` is accepted.
fn join_group(new_poly: &mut Poly, poly: &Poly, encountered_groups: &mut Vec<usize>) {
    let group = poly
        .group
        .expect("an accepted fracture has a cluster group");

    match new_poly.group {
        None => new_poly.group = Some(group),
        Some(new_poly_group) if new_poly_group != group => encountered_groups.push(group),
        Some(_) => {}
    }
}

/// True if the intersection is too near to a vertex of the polygon. The function
/// then also counts a reject in `stats.rejection_reasons.close_to_node`.
///
/// # Arguments
///
/// * `poly` - Polygon
/// * `intersection` - Intersection
/// * `min_dist` - Minimum distance between the intersection and a vertex
/// * `stats` - Statistics of the generation
/// * `eps` - Epsilon value for floating point comparisons
fn check_distance_from_nodes(
    poly: &Poly,
    intersection: &Intersection,
    min_dist: f64,
    stats: &mut Stats,
    eps: f64,
) -> bool {
    // Intersection dist to poly vertices
    let line = [intersection.p1, intersection.p2];
    let mut dist = point_to_line_seg(&poly.vertices[0], &line, eps);

    for vertex in &poly.vertices[1..] {
        let temp = point_to_line_seg(vertex, &line, eps);

        // If new distance is less than the last calculated distance...
        if temp < dist {
            dist = temp;
        }

        if dist < min_dist && dist > eps {
            stats.rejection_reasons.close_to_node += 1;
            return true;
        }
    }

    // If it makes it here, no distance less than `min_dist` found
    false
}

/// True if the intersection is too near to an edge of the polygon.
///
/// An intersection which is nearer than `min_dist` to an edge is first shortened,
/// down to `shrink_limit`. If it is then still too near to the edge, the result
/// is true.
///
/// # Arguments
///
/// * `poly` - Polygon
/// * `intersection` - Intersection. The function can shorten it.
/// * `shrink_limit` - Minimum length of the shortened intersection
/// * `stats` - Statistics of the generation
/// * `min_dist` - Minimum distance between an end point of the intersection and
///   an edge which the intersection does not cross
/// * `eps` - Epsilon value for floating point comparisons
fn check_close_edge(
    poly: &Poly,
    intersection: &mut Intersection,
    shrink_limit: f64,
    stats: &mut Stats,
    min_dist: f64,
    eps: f64,
) -> bool {
    // End points of the newest intersection
    let line = [intersection.p1, intersection.p2];
    // Counts how many end points are on the polys edge.
    // If both end points of intersection are on polys edge,
    // we must check the distance from end points to
    // vertices
    let mut on_edge_count = 0;

    // Edge from each vertex to the next one. The last edge closes the polygon.
    let next_vertices = poly.vertices.iter().cycle().skip(1);

    for (&vertex, &next_vertex) in poly.vertices.iter().zip(next_vertices) {
        let edge = [vertex, next_vertex];

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
        end_pts_to_edge.sort_by(f64::total_cmp);

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
            if shrink_intersection(intersection, &edge, shrink_limit, min_dist, min_dist, eps) {
                // Returns one if insterection shrinks to less than .9*min_feature_size
                return true;
            }
        } else if dist < eps {
            // Endpoint is almost exactly on poly's edge, must check
            // whether the discretized nodes will ne closer
            // than the minimum allowed distance
            // `shrink_intersection` will check if the descretized intersection points violate any
            // distance to edge rules
            // NOTE: Intersections discretize with set size = .5*min_feature_size
            // Minimum distance to edge must be less than .5*min_feature_size to allow for angles
            let min_dist2 = 0.4 * min_dist;

            if shrink_intersection(intersection, &edge, shrink_limit, min_dist2, min_dist, eps) {
                //returns one if insterection shrinks to less than .9*min_feature_size
                return true;
            }

            on_edge_count += 1;

            // If both end points are on the edge we must also check the distnace
            // of the intersection to the poly vertices.
            // IF the intersecion is within the polygon, distances less than min_feature_size to
            // edges and vertices will be caught by `line_seg_to_line_seg` in this function.
            // No short circuit: the node check also counts its rejects in the statistics
            if (on_edge_count >= 2)
                & check_distance_from_nodes(poly, intersection, min_dist, stats, eps)
            {
                return true;
            }
        }
    }

    #[cfg(feature = "disable_shortening_int")]
    if intersection.is_shortened {
        return true;
    }

    false
}

/// Checks the triple intersection points which a new intersection makes with the
/// intersections of the other fracture. The points go into `candidate`.
///
/// # Returns
///
/// The reject reason, if the fracture is rejected:
///
/// * `TripleIntersectionNotAllowed` - The input does not permit triple intersections
/// * `CloseToPreviousIntersection` - The intersection is too near to a previous intersection
/// * `TripleAngleTooSmall` - The angle between the two intersections is too small
/// * `TriplePointCloseToEndpoint` - The triple point is too near to an end point
/// * `TriplePointCloseToTriplePoint` - The triple point is too near to another triple point
fn check_for_triple_intersections(
    options: &PolyOptions,
    fram_options: &FramOptions,
    accepted: &AcceptedDfn,
    candidate: &mut CandidateIntersections,
    intersection: &Intersection,
    other_poly: &Poly,
) -> Result<(), RejectReason> {
    let min_feature_size = options.min_feature_size;
    let eps = options.eps;
    let intersections = accepted.intersections;
    let triple_points = accepted.triple_points;
    // Index which the new intersection gets, if it is accepted
    let new_intersection_index = intersections.len() + candidate.intersections.len();
    let candidate_triple_points = &mut candidate.triple_points;
    let min_dist = 1.5 * min_feature_size;
    let new_segment = [intersection.p1, intersection.p2]; //newest intersection

    // Check the new intersection of the new fracture against the previous intersections on other_poly
    for &old_intersection_index in &other_poly.intersection_indices {
        let line = [
            intersections[old_intersection_index].p1,
            intersections[old_intersection_index].p2,
        ];

        if are_segments_farther_than(&new_segment, &line, min_feature_size) {
            continue;
        }

        let (mut dist1, pt) = line_seg_to_line_seg(&new_segment, &line); // distance and pt of intersection

        if dist1 >= min_feature_size {
            continue;
        }

        if !fram_options.triple_intersections && dist1 < min_feature_size {
            return Err(RejectReason::TripleIntersectionNotAllowed); // Triple intersections not allowed (user input option)
        }

        if dist1 > eps && dist1 < min_feature_size {
            return Err(RejectReason::CloseToPreviousIntersection); // Point too close to other intersection
        }

        // Overlaping intersections
        if dist1 <= eps {
            // ANGLE CHECK, using definition of dot product: A dot B = Mag(A)*Mag(B) * Cos(angle)
            // Normalize  A and B first. Compare to precalculated values of cos(47deg)= .681998 and cos(133deg)= -.681998
            let u = (line[1] - line[0]).normalize();
            let v = (intersection.p2 - intersection.p1).normalize();
            let dot_prod = u.dot(&v);

            if !(-0.68199836..=0.68199836).contains(&dot_prod) {
                //if angle is less than 47 deg or greater than 133, reject
                return Err(RejectReason::TripleAngleTooSmall);
            }

            // Check that the triple intersection point isn't too close to an endpoint
            dist1 = distance(&pt, &new_segment[0]);
            let mut dist2 = distance(&pt, &new_segment[1]);

            if dist1 < min_feature_size || dist2 < min_feature_size {
                return Err(RejectReason::TriplePointCloseToEndpoint);
            }

            dist1 = distance(&pt, &line[0]);
            dist2 = distance(&pt, &line[1]);

            if dist1 < min_feature_size || dist2 < min_feature_size {
                return Err(RejectReason::TriplePointCloseToEndpoint);
            }

            // Check that the triple intersection point isn't too close to previous triple intersection points
            for &idx in &intersections[old_intersection_index].triple_point_indices {
                dist1 = distance(&pt, &triple_points[idx]);

                if dist1 < min_dist {
                    return Err(RejectReason::TriplePointCloseToTriplePoint);
                }
            }

            // Look at previously found triple points on current line of intersection,
            // do not save duplicate triple ponits
            // If duplicate point is found, a different intersection may still be needing an update
            // (three intersections per single triple intersection point )
            let duplicate = candidate_triple_points
                .iter_mut()
                .find(|saved| is_same_point(&pt, &saved.triple_point, eps));

            match duplicate {
                // The intersection of the old fracture needs a reference to the
                // triple point. It gets the reference if the new fracture is accepted.
                Some(saved) => saved.intersection_indices.push(new_intersection_index),
                // Save the triple point. It goes into the DFN if the new fracture is accepted.
                None => candidate_triple_points.push(CandidateTriplePoint {
                    triple_point: pt,
                    intersection_indices: vec![old_intersection_index, new_intersection_index],
                }),
            }
        }
    }

    // Test distances to triple points on own line of intersection
    // If any triple points are found to be less than `min_dist` to each other, reject.
    for (k, first) in candidate_triple_points.iter().enumerate() {
        for second in &candidate_triple_points[k + 1..] {
            let dist = distance(&first.triple_point, &second.triple_point);

            if dist < min_dist && dist > eps {
                return Err(RejectReason::TriplePointCloseToTriplePoint);
            }
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    const ALL_REJECT_REASONS: [RejectReason; 12] = [
        RejectReason::ShortIntersection,
        RejectReason::CloseToIntersection,
        RejectReason::CloseToEdge,
        RejectReason::SamePlane,
        RejectReason::TripleIntersectionNotAllowed,
        RejectReason::CloseToPreviousIntersection,
        RejectReason::TripleAngleTooSmall,
        RejectReason::TriplePointCloseToEndpoint,
        RejectReason::TriplePointCloseToTriplePoint,
        RejectReason::ShortBoundaryLine,
        RejectReason::BoundaryCloseToEdge,
        RejectReason::IntersectionCloseToBoundary,
    ];

    #[test]
    fn test_reject_reason_name_is_the_name_of_the_item() {
        for reason in ALL_REJECT_REASONS {
            assert_eq!(reason.name(), format!("{reason:?}"));
        }
    }
}
