use std::rc::Rc;

use parry3d_f64::na::{Point3, Vector3};
use tracing::{error, info, warn};

use crate::domain_mesh::DomainMesh;
use crate::fracture::cluster_groups::ClusterGroups;
use crate::fracture::poly::Poly;
use crate::fracture::user_defined_fractures::DomainCheck;
use crate::fram::{check_intersections, FramOptions, RejectReason};
use crate::region_boundaries::RegionBoundaryOptions;
use crate::stats::Stats;

/// A fracture which was not accepted, with the FRAM reject reason of its last attempt.
pub struct RejectedFracture {
    pub poly: Poly,
    pub reason: RejectReason,
}

/// Line of intersection between two accepted fractures
#[derive(Clone, Default)]
pub struct Intersection {
    /// Index of one of the two fractures in the accepted fractures
    pub fracture1: usize,

    /// Index of the other fracture in the accepted fractures
    pub fracture2: usize,

    /// Intersection endpoint 1
    pub p1: Point3<f64>,
    /// Intersection endpoint 2
    pub p2: Point3<f64>,

    /// Indices of the triple intersection points on the intersection
    pub triple_point_indices: Vec<usize>,

    /// True if FRAM shortened the intersection (see `Stats::discarded_length`)
    pub is_shortened: bool,
}

/// The discrete fracture network which is generated
#[derive(Default)]
pub struct Dfn {
    /// Statistics of the generation
    pub stats: Stats,

    /// Clusters of connected fractures
    pub clusters: ClusterGroups,

    /// Accepted fractures, in the order in which they were accepted
    pub accepted_polys: Vec<Poly>,

    /// Rejected fractures. Stays empty unless FRAM is on with `write_rejected_fractures`.
    pub rejected_polys: Vec<RejectedFracture>,

    /// Intersections between the accepted fractures
    pub intersections: Vec<Intersection>,

    /// Triple intersection points
    pub triple_points: Vec<Point3<f64>>,
}

impl Dfn {
    pub fn new() -> Self {
        Self::default()
    }

    /// Inserts a user-defined fracture, or rejects it.
    ///
    /// # Arguments
    ///
    /// * `new_poly` - The fracture
    /// * `poly_id` - Number of the fracture in its file, from 1
    /// * `domain_check` - Domain check of the fracture
    /// * `options` - Options of the generation
    pub fn insert_poly(
        &mut self,
        mut new_poly: Poly,
        poly_id: usize,
        domain_check: DomainCheck,
        options: &PolyOptions,
    ) {
        let is_outside_box_domain = new_poly.is_outside_box_domain(
            options.min_feature_size,
            options.eps,
            &options.domain_size,
        );
        let is_rejected = match domain_check {
            // Poly completely outside domain, or it intersects no cell of its region of the domain mesh
            DomainCheck::Regions(region_ids) => {
                is_outside_box_domain
                    || !options
                        .domain_mesh
                        .intersects_polygon_in_any_region(&new_poly.vertices, region_ids)
            }
            DomainCheck::FractureIsInDomain => false,
        };

        if is_rejected {
            self.stats.rejection_reasons.outside += 1;
            self.stats.rejected_count += 1;
            self.stats.rejected_user_fractures.push(poly_id);
            match domain_check {
                DomainCheck::Regions(region_ids) if !region_ids.is_empty() => warn!(
                    "User-defined fracture {}: rejected for being outside its regions of the domain.",
                    poly_id
                ),
                _ => warn!(
                    "User-defined fracture {}: rejected for being outside the defined domain.",
                    poly_id
                ),
            }
            return;
        }

        new_poly.assign_bounding_box();
        // Line of intersection and FRAM
        let fram_result = check_intersections(options, &mut new_poly, self);

        if let Err(reject_reason) = fram_result {
            self.stats.rejects_per_attempt[self.stats.accepted_count] += 1;
            self.stats.rejected_count += 1;
            self.stats.rejected_user_fractures.push(poly_id);

            warn!(
                "User-defined fracture {}: rejected by FRAM: {} ({})",
                poly_id,
                reject_reason.name(),
                reject_reason.description()
            );

            if options.keeps_rejected_fractures() {
                self.rejected_polys.push(RejectedFracture {
                    poly: new_poly,
                    reason: reject_reason,
                });
            }
        } else {
            // If intersection is ok (FRAM passed all tests)
            // Incriment counter of accepted polys
            self.stats.accepted_count += 1;
            // Calculate poly's area
            new_poly.assign_area();
            // Add new counter of `rejects_per_attempt`
            self.stats.rejects_per_attempt.push(0);
            self.accepted_polys.push(new_poly); // Save the new fracture to accepted polys list
        }
    }

    /// Resets intersection, triple point, and cluster group state.
    /// Used before rebuilding the DFN (e.g. after removing fractures).
    fn reset_state(&mut self) {
        self.clusters.clear();
        self.triple_points.clear();
        self.intersections.clear();
    }

    /// Removes the accepted fractures with an x radius smaller than `min_size`,
    /// and rebuilds the DFN from the fractures which stay.
    ///
    /// Use it after the generation and before the cluster analysis. The minimum
    /// radius of each family still applies during the generation, thus `min_size`
    /// has an effect only if it is larger than that minimum.
    pub fn remove_small_fractures(&mut self, options: &PolyOptions, min_size: f64) {
        self.rebuild_with(options, |poly| poly.x_radius >= min_size);
    }

    /// Rebuilds the DFN (intersections, triple points, cluster groups) from
    /// the accepted fractures for which `is_kept` is true.
    fn rebuild_with(&mut self, options: &PolyOptions, is_kept: impl Fn(&Poly) -> bool) {
        let previous_polys = std::mem::take(&mut self.accepted_polys);
        self.reset_state();

        for mut new_poly in previous_polys.into_iter().filter(|poly| is_kept(poly)) {
            new_poly.group = None; // Reset cluster group number
            new_poly.intersection_indices.clear(); // Remove ref to old intersections

            // Find line of intersection and FRAM check
            let fram_result = check_intersections(options, &mut new_poly, self);

            if fram_result.is_ok() {
                // Intersections and triple points are saved within check_intersections()
                self.accepted_polys.push(new_poly);
            } else {
                error!("Error rebuilding dfn, previously accepted fracture was rejected during DFN rebuild.");
            }
        }

        info!("Rebuilding DFN complete.");
    }
}

pub struct PolyOptions {
    pub min_feature_size: f64,
    pub eps: f64,
    pub domain_size: Vector3<f64>,
    /// The cells which make the domain.
    pub domain_mesh: Rc<DomainMesh>,
    /// Options of FRAM. `None`: FRAM is off.
    pub fram: Option<FramOptions>,
    /// Options of the boundary lines of the regions. `None`: no lines.
    pub region_boundaries: Option<RegionBoundaryOptions>,
}

impl PolyOptions {
    /// True if FRAM is on and the rejected fractures are kept for the output
    pub fn keeps_rejected_fractures(&self) -> bool {
        self.fram
            .as_ref()
            .is_some_and(|fram_options| fram_options.write_rejected_fractures)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain_mesh::test_meshes::l_shaped_mesh;
    use crate::fracture::user_defined_fractures::{insert_user_defined_fractures, UserFracture};

    /// Options with the L-shaped mesh as the domain. The mesh has no cell
    /// where x and y are positive.
    fn poly_options(keep_rejected_fractures: bool) -> PolyOptions {
        PolyOptions {
            min_feature_size: 0.1,
            eps: 1e-9,
            domain_size: Vector3::new(10., 10., 10.),
            domain_mesh: Rc::new(l_shaped_mesh()),
            fram: Some(FramOptions {
                write_rejected_fractures: keep_rejected_fractures,
                ..fram_options()
            }),
            region_boundaries: None,
        }
    }

    /// FRAM on, with the check of two user-defined fractures
    fn fram_options() -> FramOptions {
        FramOptions {
            relaxed: false,
            triple_intersections: false,
            check_user_fractures: true,
            write_rejected_fractures: false,
            rejects_per_fracture: 10,
        }
    }

    /// Inserts the same rectangle two times. FRAM rejects the second one.
    fn insert_two_equal_rectangles(keep_rejected_fractures: bool) -> Dfn {
        let options = poly_options(keep_rejected_fractures);
        let mut dfn = Dfn::new();

        for poly_id in 1..=2 {
            let mut rectangle = Poly::new(4, 1., 1.);
            rectangle.normal = Vector3::new(0., 0., 1.);
            dfn.insert_poly(rectangle, poly_id, DomainCheck::Regions(&[]), &options);
        }

        dfn
    }

    /// Inserts a horizontal square with sides of 2 into the L-shaped domain mesh.
    /// The mesh has no cell where x and y are positive.
    fn insert_square_into_l_shaped_mesh(center: [f64; 3]) -> Dfn {
        insert_square_into_regions(center, &[])
    }

    /// As `insert_square_into_l_shaped_mesh`, for a square which must intersect one of the regions.
    /// Region 1 is where y is negative, and region 2 is where x is negative and y is positive.
    fn insert_square_into_regions(center: [f64; 3], region_ids: &[usize]) -> Dfn {
        let options = poly_options(false);
        let mut dfn = Dfn::new();

        let mut square = Poly::new(4, 2f64.sqrt(), 1.);
        square.normal = Vector3::new(0., 0., 1.);
        square.translate(Vector3::from(center));
        dfn.insert_poly(square, 1, DomainCheck::Regions(region_ids), &options);

        dfn
    }

    #[test]
    fn test_user_fracture_with_center_outside_the_mesh_is_accepted_if_it_intersects_a_cell() {
        let dfn = insert_square_into_l_shaped_mesh([0.5, 0.5, 0.]);

        assert_eq!(dfn.accepted_polys.len(), 1);
    }

    #[test]
    fn test_user_fracture_which_intersects_a_cell_of_its_region_is_accepted() {
        let dfn = insert_square_into_regions([2.5, -2.5, 0.], &[1]);

        assert_eq!(dfn.accepted_polys.len(), 1);
    }

    #[test]
    fn test_user_fracture_which_intersects_only_cells_of_another_region_is_rejected() {
        let dfn = insert_square_into_regions([2.5, -2.5, 0.], &[2]);

        assert_eq!(dfn.stats.rejected_user_fractures, [1]);
    }

    #[test]
    fn test_user_fracture_which_intersects_a_cell_of_one_of_its_regions_is_accepted() {
        let dfn = insert_square_into_regions([2.5, -2.5, 0.], &[2, 1]);

        assert_eq!(dfn.accepted_polys.len(), 1);
    }

    #[test]
    fn test_user_fracture_which_intersects_no_cell_is_rejected() {
        let dfn = insert_square_into_l_shaped_mesh([2.5, 2.5, 0.]);

        assert_eq!(dfn.stats.rejected_user_fractures, [1]);
    }

    #[test]
    fn test_rejected_user_fracture_is_kept_when_the_option_is_on() {
        let dfn = insert_two_equal_rectangles(true);

        assert_eq!(dfn.accepted_polys.len(), 1);
        assert_eq!(dfn.rejected_polys.len(), 1);
        assert_eq!(dfn.rejected_polys[0].reason, RejectReason::CloseToEdge);
    }

    #[test]
    fn test_reject_after_the_first_accepted_fracture_is_counted_for_the_second_attempt() {
        let dfn = insert_two_equal_rectangles(false);

        assert_eq!(dfn.stats.rejects_per_attempt, [0, 1]);
    }

    /// Two horizontal squares with the shared edge x = -2, and a vertical
    /// square which crosses that edge. All are in the L-shaped domain mesh.
    const LEFT_PIECE: [[f64; 3]; 4] = [
        [-3., -3., 0.],
        [-2., -3., 0.],
        [-2., -2., 0.],
        [-3., -2., 0.],
    ];
    const RIGHT_PIECE: [[f64; 3]; 4] = [
        [-2., -3., 0.],
        [-1., -3., 0.],
        [-1., -2., 0.],
        [-2., -2., 0.],
    ];
    /// `RIGHT_PIECE` with a gap of 0.05 (less than min_feature_size) to `LEFT_PIECE`
    const RIGHT_PIECE_AFTER_GAP: [[f64; 3]; 4] = [
        [-1.95, -3., 0.],
        [-0.95, -3., 0.],
        [-0.95, -2., 0.],
        [-1.95, -2., 0.],
    ];
    const CROSSING_SQUARE: [[f64; 3]; 4] = [
        [-2.5, -2.5, -0.5],
        [-1.5, -2.5, -0.5],
        [-1.5, -2.5, 0.5],
        [-2.5, -2.5, 0.5],
    ];

    /// Inserts user-defined polygons, each with its value of the `fracture_id` array.
    fn insert_user_polygons(polygons: &[([[f64; 3]; 4], Option<u64>)]) -> Dfn {
        insert_user_polygons_with_options(polygons, &poly_options(false))
    }

    fn insert_user_polygons_with_options(
        polygons: &[([[f64; 3]; 4], Option<u64>)],
        options: &PolyOptions,
    ) -> Dfn {
        let fractures: Vec<UserFracture> = polygons
            .iter()
            .map(|(vertices, fracture_id)| UserFracture {
                vertices: vertices.to_vec(),
                region_ids: Vec::new(),
                fracture_id: *fracture_id,
            })
            .collect();
        let mut dfn = Dfn::new();

        insert_user_defined_fractures(&fractures, options, &mut dfn);

        dfn
    }

    /// A horizontal square in the bounding box of the L-shaped mesh which
    /// intersects no cell of it.
    const PIECE_OUTSIDE_MESH: [[f64; 3]; 4] =
        [[2., 2., 0.], [3., 2., 0.], [3., 3., 0.], [2., 3., 0.]];

    #[test]
    fn test_piece_outside_the_mesh_is_accepted_with_a_piece_of_its_fracture_in_the_mesh() {
        let dfn = insert_user_polygons(&[(PIECE_OUTSIDE_MESH, Some(5)), (LEFT_PIECE, Some(5))]);

        assert_eq!(dfn.accepted_polys.len(), 2);
    }

    #[test]
    fn test_polygon_outside_the_mesh_is_rejected_with_only_another_fracture_in_the_mesh() {
        let dfn = insert_user_polygons(&[(PIECE_OUTSIDE_MESH, Some(5)), (LEFT_PIECE, Some(6))]);

        assert_eq!(dfn.stats.rejected_user_fractures, [1]);
    }

    #[test]
    fn test_fracture_with_all_pieces_outside_the_mesh_is_rejected() {
        let mut second_piece = PIECE_OUTSIDE_MESH;
        for vertex in &mut second_piece {
            vertex[0] += 1.;
        }

        let dfn = insert_user_polygons(&[(PIECE_OUTSIDE_MESH, Some(5)), (second_piece, Some(5))]);

        assert_eq!(dfn.stats.rejected_user_fractures, [1, 2]);
    }

    /// A vertical square which crosses `LEFT_PIECE` at a distance of 0.05
    /// (less than min_feature_size) from its edge x = -2.
    const SQUARE_NEAR_EDGE: [[f64; 3]; 4] = [
        [-2.05, -2.8, -0.5],
        [-2.05, -2.2, -0.5],
        [-2.05, -2.2, 0.5],
        [-2.05, -2.8, 0.5],
    ];

    /// Inserts `LEFT_PIECE` and `SQUARE_NEAR_EDGE` as two user-defined fractures.
    fn insert_user_fracture_near_edge(check_user_fractures: bool) -> Dfn {
        let options = PolyOptions {
            fram: Some(FramOptions {
                check_user_fractures,
                ..fram_options()
            }),
            ..poly_options(false)
        };

        insert_user_polygons_with_options(&[(LEFT_PIECE, None), (SQUARE_NEAR_EDGE, None)], &options)
    }

    #[test]
    fn test_neighbor_polygon_of_another_fracture_is_accepted_without_the_user_fracture_check() {
        let options = PolyOptions {
            fram: Some(FramOptions {
                check_user_fractures: false,
                ..fram_options()
            }),
            ..poly_options(false)
        };

        let dfn =
            insert_user_polygons_with_options(&[(LEFT_PIECE, None), (RIGHT_PIECE, None)], &options);

        assert_eq!(dfn.accepted_polys.len(), 2);
    }

    #[test]
    fn test_user_fracture_near_the_edge_of_a_user_fracture_is_rejected() {
        let dfn = insert_user_fracture_near_edge(true);

        assert_eq!(dfn.stats.rejected_user_fractures, [2]);
    }

    #[test]
    fn test_user_fracture_near_the_edge_of_a_user_fracture_is_accepted_without_the_user_fracture_check(
    ) {
        let dfn = insert_user_fracture_near_edge(false);

        assert_eq!(dfn.accepted_polys.len(), 2);
        assert_eq!(dfn.intersections.len(), 1);
    }

    #[test]
    fn test_pieces_of_one_fracture_have_no_intersection() {
        let dfn = insert_user_polygons(&[(LEFT_PIECE, Some(0)), (RIGHT_PIECE, Some(0))]);

        assert_eq!(dfn.accepted_polys.len(), 2);
        assert!(dfn.intersections.is_empty());
    }

    #[test]
    fn test_pieces_of_one_fracture_are_in_one_cluster() {
        let dfn = insert_user_polygons(&[(LEFT_PIECE, Some(0)), (RIGHT_PIECE, Some(0))]);

        assert_eq!(dfn.accepted_polys[0].group, dfn.accepted_polys[1].group);
    }

    #[test]
    fn test_neighbor_polygon_of_another_fracture_is_rejected() {
        let dfn = insert_user_polygons(&[(LEFT_PIECE, Some(0)), (RIGHT_PIECE, Some(1))]);

        assert_eq!(dfn.stats.rejected_user_fractures, [2]);
    }

    #[test]
    fn test_fracture_can_cross_the_shared_edge_of_two_pieces() {
        let dfn = insert_user_polygons(&[
            (LEFT_PIECE, Some(0)),
            (RIGHT_PIECE, Some(0)),
            (CROSSING_SQUARE, None),
        ]);

        assert_eq!(dfn.accepted_polys.len(), 3);
        assert_eq!(dfn.intersections.len(), 2);
    }

    #[test]
    fn test_fracture_cannot_cross_two_fractures_with_a_gap_less_than_h() {
        let dfn = insert_user_polygons(&[
            (LEFT_PIECE, Some(0)),
            (RIGHT_PIECE_AFTER_GAP, Some(1)),
            (CROSSING_SQUARE, None),
        ]);

        assert_eq!(dfn.stats.rejected_user_fractures, [3]);
    }

    #[test]
    fn test_rejected_user_fracture_is_not_kept_when_the_option_is_off() {
        let dfn = insert_two_equal_rectangles(false);

        assert_eq!(dfn.accepted_polys.len(), 1);
        assert!(dfn.rejected_polys.is_empty());
    }
}
