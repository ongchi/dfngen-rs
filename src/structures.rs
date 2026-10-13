use std::cell::RefCell;
use std::rc::Rc;

use parry3d_f64::na::{Point3, Vector3};
use rand::Rng;
use rand_mt::Mt64;
use tracing::{error, info, warn};

use crate::cg::intersection_checking;
use crate::distribution::{LogNormalSampling, TruncExp, TruncLogNormal, TruncPowerLaw};
use crate::domain_mesh::DomainMesh;
use crate::error::DfngenError;
use crate::fracture::insert_shape::print_reject_reason;
use crate::fracture::poly::Poly;
use crate::fracture::user_defined_fractures::DomainCheck;

/// A fracture which was not accepted, with the FRAM reject code of its last attempt.
pub struct RejectedFracture {
    pub poly: Poly,
    pub reject_code: i32,
}

/// Intersections structure.
/// This structure contains all data pertaining to
/// one intersection. This includes the IDs for both
/// intersecting fractures, the intersection end points,
/// a list of references to any triple intersection points existing
/// on the intersection, and a flag denotting whether or not
/// the intersection has been shortened by shrinkIntersection()
/// during FRAM.
#[derive(Clone, Default)]
pub struct IntersectionPoints {
    /// Fracture 1, index of fracture in the accpeted polygons list that this
    /// intersection belongs to ('fract1' and 'fract2' are in no particular order).
    pub fract1: isize,

    /// Fracture 2, index of fracture in the accpeted polygons list that this
    /// intersection belongs to ('fract1' and 'fract2' are in no particular order).
    pub fract2: isize,

    /// Intersection endpoint 1
    pub p1: Point3<f64>,
    /// Intersection endpoint 2
    pub p2: Point3<f64>,

    /// Triple intersection points/nodes on intersection.
    pub triple_points_idx: Vec<usize>,

    /// Used to update book keeping for keeping track of overal intersection length
    /// that has been shortened from shrinkIntersection(). Used in intersectionChecking().
    pub intersection_shortened: bool,
}

/// Holds temporary triple point data while FRAM is checking
/// all intersections for a new polygon/fracture.
///
/// Once a fracture is accepted, the temporary triple point
/// 'triplePoint' is moved to its permanent location.
///
/// 'intIdx' is used to update intersection point structures' (IntPoints)
/// rejerences to the new triple intersection points. 'intIdx' contains
/// the index of the triple point in the permanent triple points array
/// IF the fracture is accepted.
///
/// If the fracture is rejected, this data is discraded.
pub struct TriplePtTempData {
    /// Triple intersection point.
    pub triple_point: Point3<f64>,
    /// Index to 'triplePoints' array in main() to where this point would be stored
    /// if the FRAM checks pass.
    pub int_index: Vec<usize>,
}

/// FractureGroups is a structure used to keep track of which fractures are in
/// each cluster. FractureGroups works in conjunction with GroupData.
///
/// FractGroups holds pointers/index numbers to the polygons belonging to  group number 'groupNum'. Unlike
/// GroupData, fractGroups does not stay aligned to cluster group numbers. To keep from copying, deleteing, and
/// re-allocating memory when groups merge together, we simply change the variable 'groupNum' to the new group
/// number and leave the polygon pointer/index list as is.
///
/// Because of this, there will be multiple FractureGroups objects with the same group number when
/// fracture clusters merge together, but with differnt polygons listed. To get all the
/// polygons from a group we must search the fractGroups array for all matching groups and look at each of their
/// polgon lists.
pub struct FractureGroups {
    /// Fracture cluster group number.
    pub group_num: usize,
    /// List of polygon indices in the 'acceptedPoly' array in main() which belong to this group.
    pub poly_list: Vec<usize>,
}

/// GroupData works in conjunction with 'FractureGroups'
/// GroupData keeps track of which group numbers ('groupNum' in above struct)
/// Connect to which boundaries. When a cluster of fractures bridges two
/// groups/clusters, One of the groups 'valid' bool is set to 0 (not valid).
/// The fractures whos valid bool becomes 0 are moved to the other group
/// In the 'GroupData' structure array, FractureGroups.groupNum-1 is the index
/// to the corresponding GroupData struct in the GroupData array.
#[derive(Default)]
pub struct GroupData {
    /// Number of polygons in group.
    pub size: usize,
    /// Valid bit, True if group this structures data is still valid, false otherwise.
    /// Data can become invalid when fracture cluster groups merge together.
    pub valid: bool,
    /// Domain boundary sides/faces that this cluster connects to..
    /// Index Key:
    /// [0]: -x face, [1]: +x face
    /// [2]: -y face, [3]: +y face
    /// [4]: -z face, [5]: +z face
    pub faces: [bool; 6],
}

/// Rejection reason counters.
#[derive(Default)]
pub struct RejectionReasons {
    /// Rejections due to intersection of length less than min_feature_size.
    pub short_intersection: usize,
    /// Rejections due to intersections being too close to close to
    /// polygon verties.
    pub close_to_node: usize,
    /// Rejections due to intersections being too close to polygon edges.
    pub close_to_edge: usize,
    /// Counter no longer in use.
    pub close_point_to_edge: usize,
    /// Rejections due to fractures landing outside of the domain.
    pub outside: usize,
    /// Rejections due to triple intersection problem.
    pub triple: usize,
    /// Rejections due to an intersection landing too close to another
    /// intersection.
    pub inter_close_to_inter: usize,
}

/// TODO: Make singleton
/// Program and DFN statisistics structure. Keeps various statistics, including
/// fracture cluster information, about the DFN being generated.
#[derive(Default)]
pub struct Stats {
    /// Counters for the number of polygons/fractures accepted by each stochastic
    /// family. Elements in this array are in the same order as the stochastic shape
    /// families array 'shapeFamilies' in main(). e.g. The counter for the second family
    /// in the shapeFamilies array is the second element in this array.
    pub accepted_from_fam: Vec<usize>,

    /// Counters for the number of polygons/fractures rejected by each stochastic
    /// family. Elements in this array are in the same order as the stochastic shape
    /// families array 'shapeFamilies' in main(). e.g. The counter for the second family
    /// in the shapeFamilies array is the second element in this array.
    pub rejected_from_fam: Vec<usize>,

    /// Number of fractures estimated for each stochastic family (see dryRun()).
    /// Elements in this array are in the same order as the stochastic shape
    /// families array 'shapeFamilies' in main(). e.g. The estimated number of
    /// fractures for the second family in the shapeFamilies array is the second
    /// element in this array.
    pub expected_from_fam: Vec<usize>,

    /// Total number of accepted polygons/fractures for the DFN. This variable is
    /// updated as the DFN is generated.
    pub accepted_poly_count: usize,

    /// Total number of rejected polygons/fractures for the DFN. This variable is
    /// updated as the DFN is generated.
    pub rejected_poly_count: usize,

    /// Total number of polygon/fracture re-translations. This variable is updated
    /// as the DFN is generated.
    pub retranslated_poly_count: u32,

    /// Total area of fractures before isolated fracture removal. Variable is set in main() after
    /// DFN generation has completed.
    pub area_before_removal: f64,

    /// Total area of fractures after isolated fracture removal. Variable is set in main() after
    /// DFN generation has completed.
    pub area_after_removal: f64,

    /// Used to assign group/cluster numbers. 'nextGroupNum initializes to 1, and is used to
    /// assign a fracture a group number, and is then incremented.
    pub next_group_num: usize,

    /// RejectionReasons structure holds counters for all rejection reasons.
    pub rejection_reasons: RejectionReasons,

    /// Counter for number of intersections that have been shortened by FRAM's
    /// shrinkIntersection() function. Only counts intersections of fractures that have been
    /// accepted into the domain.
    pub intersections_shortened: u32,

    /// Total length of original intersections in the DFN before intersections were shortened.
    /// Length of all intersections if none were shortened
    pub original_length: f64,

    /// Total length of intersections in DFN which were discarded by shrinkIntersection()
    /// inside FRAM.
    /// Final intersection length = originalLength - discardedLength
    pub discarded_length: f64,

    /// Rejects per insertion attempt counter. Each element represents the number of tries
    /// it took to fit a fracture into the DFN. e.g. The number stored in the 10th element
    /// is the number of tries it took before the 10th fracture was accepted. This count
    /// inclues re-translating the same fracture to different locations as well as generating
    /// new fractures.
    pub rejects_per_attempt: Vec<usize>,

    /// Fracture cluster data. See struct FractureGroup.
    pub fract_group: Vec<FractureGroups>,
    /// Fracture cluster data. See struct GroupData.
    pub group_data: Vec<GroupData>,
    /// Rejected User Fractures
    pub rejected_user_fracture: Vec<usize>,
}

pub enum RadiusFunction {
    LogNormal {
        mu: f64,
        sigma: f64,
        sampling: LogNormalSampling,
    },
    TruncatedPowerLaw {
        alpha: f64,
    },
    Exponential {
        lambda: f64,
    },
    Constant(f64),
}

pub struct RadiusDistribution {
    pub min: f64,
    pub max: f64,
    pub function: RadiusFunction,
}

impl RadiusDistribution {
    pub fn new_truncated_power_law(alpha: f64, min: f64, max: f64) -> Self {
        Self {
            function: RadiusFunction::TruncatedPowerLaw { alpha },
            min,
            max,
        }
    }

    pub fn new_log_normal(
        mu: f64,
        sigma: f64,
        min: f64,
        max: f64,
        sampling: LogNormalSampling,
    ) -> Self {
        Self {
            function: RadiusFunction::LogNormal {
                mu,
                sigma,
                sampling,
            },
            min,
            max,
        }
    }

    pub fn new_exponential(lambda: f64, min: f64, max: f64) -> Self {
        Self {
            function: RadiusFunction::Exponential { lambda },
            min,
            max,
        }
    }

    pub fn new_constant(value: f64) -> Self {
        Self {
            function: RadiusFunction::Constant(value),
            min: value,
            max: value,
        }
    }

    /// Sampling from distribution function.
    pub fn sample(&self, generator: Rc<RefCell<Mt64>>) -> Result<f64, DfngenError> {
        let radius = match self.function {
            RadiusFunction::LogNormal {
                mu,
                sigma,
                sampling,
            } => {
                let distr = TruncLogNormal::new(self.min, self.max, mu, sigma, sampling)?;
                generator.borrow_mut().sample(&distr)?
            }
            RadiusFunction::TruncatedPowerLaw { alpha } => {
                let distr = TruncPowerLaw::new(self.min, self.max, alpha);
                generator.borrow_mut().sample(&distr)
            }
            RadiusFunction::Exponential { lambda } => {
                let distr = TruncExp::new(self.min, self.max, lambda)?;
                generator.borrow_mut().sample(&distr)?
            }
            RadiusFunction::Constant(c) => c,
        };

        Ok(radius)
    }
}

impl IntersectionPoints {
    /// Initializes fract1 and fract2 to -1.
    /// Sets intersectionShortened to false;
    pub fn new() -> Self {
        Self {
            fract1: -1,
            fract2: -1,
            ..Default::default()
        }
    }
}

impl FractureGroups {
    /// Initializes polyList vector to reserve enough memory for 100 floats.
    pub fn new(group_num: usize) -> Self {
        Self {
            poly_list: Vec::new(),
            group_num,
        }
    }
}

impl GroupData {
    /// Initializes size to zero, valid to true,
    /// and zeros (set to false) the faces array.
    pub fn new() -> Self {
        Self {
            valid: true,
            ..Default::default()
        }
    }
}

impl Stats {
    /// Initializes all counters to zero. Initializes
    /// nextGroupNum to 1. (See implementation for details)
    pub fn new() -> Self {
        Self {
            next_group_num: 1,
            // Counter of the rejects before the first accepted fracture
            rejects_per_attempt: vec![0],
            ..Default::default()
        }
    }
}

pub struct DFNGen {
    // Statistics structure:
    pub pstats: Stats,

    // Vector to store accepted polygons/fractures
    pub accepted_poly: Vec<Poly>,

    // Rejected fractures. Stays empty unless PolyOptions::keep_rejected_fractures is set.
    pub rejected_poly: Vec<RejectedFracture>,

    // Vector for storing intersections
    pub intpts: Vec<IntersectionPoints>,

    // Vector for storing triple intersection points
    pub triple_points: Vec<Point3<f64>>,
}

impl DFNGen {
    pub fn new() -> Self {
        Self {
            pstats: Stats::new(),
            accepted_poly: Vec::new(),
            rejected_poly: Vec::new(),
            intpts: Vec::new(),
            triple_points: Vec::new(),
        }
    }

    /// Inserts a user-defined fracture, or rejects it.
    ///
    /// # Arguments
    ///
    /// * `new_poly` - The fracture
    /// * `poly_id` - Number of the fracture in its file, from 1
    /// * `domain_check` - Domain check of the fracture
    /// * `opts` - Options of the generation
    pub fn insert_poly(
        &mut self,
        mut new_poly: Poly,
        poly_id: usize,
        domain_check: DomainCheck,
        opts: &PolyOptions,
    ) {
        let is_outside_box_domain =
            new_poly.check_box_domain(opts.min_feature_size, opts.eps, &opts.domain_size);
        let is_rejected = match domain_check {
            // Poly completely outside domain, or it intersects no cell of its region of the domain mesh
            DomainCheck::Regions(region_ids) => {
                is_outside_box_domain
                    || !opts
                        .domain_mesh
                        .intersects_polygon_in_any_region(&new_poly.vertices, region_ids)
            }
            DomainCheck::FractureIsInDomain => false,
        };

        if is_rejected {
            self.pstats.rejection_reasons.outside += 1;
            self.pstats.rejected_poly_count += 1;
            self.pstats.rejected_user_fracture.push(poly_id);
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
        let reject_code = intersection_checking(
            opts,
            &mut new_poly,
            &mut self.accepted_poly,
            &mut self.intpts,
            &mut self.pstats,
            &mut self.triple_points,
        );

        if reject_code == 0 {
            // If intersection is ok (FRAM passed all tests)
            // Incriment counter of accepted polys
            self.pstats.accepted_poly_count += 1;
            // Calculate poly's area
            new_poly.assign_area();
            // Add new rejectsPerAttempt counter
            self.pstats.rejects_per_attempt.push(0);
            self.accepted_poly.push(new_poly); // Save newPoly to accepted polys list
        } else {
            self.pstats.rejects_per_attempt[self.pstats.accepted_poly_count] += 1;
            self.pstats.rejected_poly_count += 1;
            self.pstats.rejected_user_fracture.push(poly_id);

            print_reject_reason(reject_code, &new_poly);

            if opts.keep_rejected_fractures {
                self.rejected_poly.push(RejectedFracture {
                    poly: new_poly,
                    reject_code,
                });
            }
        }
    }

    /// Resets intersection, triple point, and cluster group state.
    /// Used before rebuilding the DFN (e.g. after removing fractures).
    fn reset_state(&mut self) {
        self.pstats.group_data.clear();
        self.pstats.fract_group.clear();
        self.triple_points.clear();
        self.intpts.clear();
        self.pstats.next_group_num = 1;
    }

    /// Remove Fractures Smaller Than Minimum Size
    ///
    /// Function is designed to be used AFTER DFN generation
    /// Originally created to compare the difference in distributions
    /// if small fractures were removed after a DFN was generated
    /// as opposed limiting their insertion during DFN generation.
    ///
    /// The minimum size options in the input files will still be used.
    /// If the user wishes to also removeFractures after DFN generation,
    /// 'minSize' here must be larger than the minimum size fractures in the
    /// input file. Fratrues with radii less than 'minSize' will be removed
    /// after DFN has been created.
    ///
    /// NOTE:
    /// Must be executed before getCluster()
    /// This funciton rebuilds the DFN. Using getCluster() before this
    /// funciton executes causes undefined behavior.
    pub fn remove_small_fractures(&mut self, opts: &PolyOptions, min_size: f64) {
        self.rebuild_with(opts, |poly| poly.xradius >= min_size);
    }

    /// Rebuilds the DFN (intersections, triple points, cluster groups) from
    /// the accepted fractures for which `is_kept` is true.
    fn rebuild_with(&mut self, opts: &PolyOptions, is_kept: impl Fn(&Poly) -> bool) {
        let mut final_poly_list: Vec<Poly> = Vec::new();
        self.reset_state();

        for poly in self.accepted_poly.iter().filter(|poly| is_kept(poly)) {
            let mut new_poly = poly.clone();
            new_poly.group_num = 0; // Reset cluster group number
            new_poly.intersection_index.clear(); // Remove ref to old intersections

            // Find line of intersection and FRAM check
            let reject_code = intersection_checking(
                opts,
                &mut new_poly,
                &mut final_poly_list,
                &mut self.intpts,
                &mut self.pstats,
                &mut self.triple_points,
            );

            if reject_code == 0 {
                // Intersections and triple points are saved within intersection_checking()
                final_poly_list.push(new_poly);
            } else {
                error!("Error rebuilding dfn, previously accepted fracture was rejected during DFN rebuild.");
            }
        }

        info!("Rebuilding DFN complete.");
        self.accepted_poly = final_poly_list;
    }
}

pub struct PolyOptions {
    pub min_feature_size: f64,
    pub eps: f64,
    pub domain_size: Vector3<f64>,
    /// The cells which make the domain.
    pub domain_mesh: Rc<DomainMesh>,
    pub r_fram: bool,
    pub is_fram_enabled: bool,
    pub triple_intersections: bool,
    /// False - FRAM does not check an intersection of two user-defined fractures.
    pub check_user_fractures: bool,
    pub keep_rejected_fractures: bool,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain_mesh::test_meshes::l_shaped_mesh;
    use crate::fracture::user_defined_fractures::{insert_user_defined_fractures, UserFracture};

    fn assert_samples_are_within_limits(distribution: RadiusDistribution) {
        let generator = Rc::new(RefCell::new(Mt64::new(1)));

        for _ in 0..10_000 {
            let radius = distribution.sample(generator.clone()).unwrap();
            assert!(
                (distribution.min..=distribution.max).contains(&radius),
                "radius {radius} is not in [{}, {}]",
                distribution.min,
                distribution.max
            );
        }
    }

    #[test]
    fn test_log_normal_radius_is_within_min_and_max() {
        assert_samples_are_within_limits(RadiusDistribution::new_log_normal(
            1.,
            1.,
            1.,
            4.,
            LogNormalSampling::Rejection,
        ));
    }

    #[test]
    fn test_exponential_radius_is_within_min_and_max() {
        assert_samples_are_within_limits(RadiusDistribution::new_exponential(0.4, 1., 4.));
    }

    #[test]
    fn test_truncated_power_law_radius_is_within_min_and_max() {
        assert_samples_are_within_limits(RadiusDistribution::new_truncated_power_law(2.5, 1., 4.));
    }

    /// Options with the L-shaped mesh as the domain. The mesh has no cell
    /// where x and y are positive.
    fn poly_options(keep_rejected_fractures: bool) -> PolyOptions {
        PolyOptions {
            min_feature_size: 0.1,
            eps: 1e-9,
            domain_size: Vector3::new(10., 10., 10.),
            domain_mesh: Rc::new(l_shaped_mesh()),
            r_fram: false,
            is_fram_enabled: true,
            triple_intersections: false,
            check_user_fractures: true,
            keep_rejected_fractures,
        }
    }

    /// Inserts the same rectangle two times. FRAM rejects the second one.
    fn insert_two_equal_rectangles(keep_rejected_fractures: bool) -> DFNGen {
        let opts = poly_options(keep_rejected_fractures);
        let mut dfngen = DFNGen::new();

        for poly_id in 1..=2 {
            let mut rectangle = Poly::new(4, 1., 1.);
            rectangle.normal = Vector3::new(0., 0., 1.);
            dfngen.insert_poly(rectangle, poly_id, DomainCheck::Regions(&[]), &opts);
        }

        dfngen
    }

    /// Inserts a horizontal square with sides of 2 into the L-shaped domain mesh.
    /// The mesh has no cell where x and y are positive.
    fn insert_square_into_l_shaped_mesh(center: [f64; 3]) -> DFNGen {
        insert_square_into_regions(center, &[])
    }

    /// As `insert_square_into_l_shaped_mesh`, for a square which must intersect one of the regions.
    /// Region 1 is where y is negative, and region 2 is where x is negative and y is positive.
    fn insert_square_into_regions(center: [f64; 3], region_ids: &[usize]) -> DFNGen {
        let opts = poly_options(false);
        let mut dfngen = DFNGen::new();

        let mut square = Poly::new(4, 2f64.sqrt(), 1.);
        square.normal = Vector3::new(0., 0., 1.);
        square.translate(Vector3::from(center));
        dfngen.insert_poly(square, 1, DomainCheck::Regions(region_ids), &opts);

        dfngen
    }

    #[test]
    fn test_user_fracture_with_center_outside_the_mesh_is_accepted_if_it_intersects_a_cell() {
        let dfngen = insert_square_into_l_shaped_mesh([0.5, 0.5, 0.]);

        assert_eq!(dfngen.accepted_poly.len(), 1);
    }

    #[test]
    fn test_user_fracture_which_intersects_a_cell_of_its_region_is_accepted() {
        let dfngen = insert_square_into_regions([2.5, -2.5, 0.], &[1]);

        assert_eq!(dfngen.accepted_poly.len(), 1);
    }

    #[test]
    fn test_user_fracture_which_intersects_only_cells_of_another_region_is_rejected() {
        let dfngen = insert_square_into_regions([2.5, -2.5, 0.], &[2]);

        assert_eq!(dfngen.pstats.rejected_user_fracture, [1]);
    }

    #[test]
    fn test_user_fracture_which_intersects_a_cell_of_one_of_its_regions_is_accepted() {
        let dfngen = insert_square_into_regions([2.5, -2.5, 0.], &[2, 1]);

        assert_eq!(dfngen.accepted_poly.len(), 1);
    }

    #[test]
    fn test_user_fracture_which_intersects_no_cell_is_rejected() {
        let dfngen = insert_square_into_l_shaped_mesh([2.5, 2.5, 0.]);

        assert_eq!(dfngen.pstats.rejected_user_fracture, [1]);
    }

    #[test]
    fn test_rejected_user_fracture_is_kept_when_the_option_is_on() {
        let dfngen = insert_two_equal_rectangles(true);

        assert_eq!(dfngen.accepted_poly.len(), 1);
        assert_eq!(dfngen.rejected_poly.len(), 1);
        assert_ne!(dfngen.rejected_poly[0].reject_code, 0);
    }

    #[test]
    fn test_reject_after_the_first_accepted_fracture_is_counted_for_the_second_attempt() {
        let dfngen = insert_two_equal_rectangles(false);

        assert_eq!(dfngen.pstats.rejects_per_attempt, [0, 1]);
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
    fn insert_user_polygons(polygons: &[([[f64; 3]; 4], Option<u64>)]) -> DFNGen {
        insert_user_polygons_with_options(polygons, &poly_options(false))
    }

    fn insert_user_polygons_with_options(
        polygons: &[([[f64; 3]; 4], Option<u64>)],
        opts: &PolyOptions,
    ) -> DFNGen {
        let fractures: Vec<UserFracture> = polygons
            .iter()
            .map(|(vertices, fracture_id)| UserFracture {
                vertices: vertices.to_vec(),
                region_ids: Vec::new(),
                fracture_id: *fracture_id,
            })
            .collect();
        let mut dfngen = DFNGen::new();

        insert_user_defined_fractures(&fractures, opts, &mut dfngen);

        dfngen
    }

    /// A horizontal square in the bounding box of the L-shaped mesh which
    /// intersects no cell of it.
    const PIECE_OUTSIDE_MESH: [[f64; 3]; 4] =
        [[2., 2., 0.], [3., 2., 0.], [3., 3., 0.], [2., 3., 0.]];

    #[test]
    fn test_piece_outside_the_mesh_is_accepted_with_a_piece_of_its_fracture_in_the_mesh() {
        let dfngen = insert_user_polygons(&[(PIECE_OUTSIDE_MESH, Some(5)), (LEFT_PIECE, Some(5))]);

        assert_eq!(dfngen.accepted_poly.len(), 2);
    }

    #[test]
    fn test_polygon_outside_the_mesh_is_rejected_with_only_another_fracture_in_the_mesh() {
        let dfngen = insert_user_polygons(&[(PIECE_OUTSIDE_MESH, Some(5)), (LEFT_PIECE, Some(6))]);

        assert_eq!(dfngen.pstats.rejected_user_fracture, [1]);
    }

    #[test]
    fn test_fracture_with_all_pieces_outside_the_mesh_is_rejected() {
        let mut second_piece = PIECE_OUTSIDE_MESH;
        for vertex in &mut second_piece {
            vertex[0] += 1.;
        }

        let dfngen =
            insert_user_polygons(&[(PIECE_OUTSIDE_MESH, Some(5)), (second_piece, Some(5))]);

        assert_eq!(dfngen.pstats.rejected_user_fracture, [1, 2]);
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
    fn insert_user_fracture_near_edge(check_user_fractures: bool) -> DFNGen {
        let opts = PolyOptions {
            check_user_fractures,
            ..poly_options(false)
        };

        insert_user_polygons_with_options(&[(LEFT_PIECE, None), (SQUARE_NEAR_EDGE, None)], &opts)
    }

    #[test]
    fn test_neighbor_polygon_of_another_fracture_is_accepted_without_the_user_fracture_check() {
        let opts = PolyOptions {
            check_user_fractures: false,
            ..poly_options(false)
        };

        let dfngen =
            insert_user_polygons_with_options(&[(LEFT_PIECE, None), (RIGHT_PIECE, None)], &opts);

        assert_eq!(dfngen.accepted_poly.len(), 2);
    }

    #[test]
    fn test_user_fracture_near_the_edge_of_a_user_fracture_is_rejected() {
        let dfngen = insert_user_fracture_near_edge(true);

        assert_eq!(dfngen.pstats.rejected_user_fracture, [2]);
    }

    #[test]
    fn test_user_fracture_near_the_edge_of_a_user_fracture_is_accepted_without_the_user_fracture_check(
    ) {
        let dfngen = insert_user_fracture_near_edge(false);

        assert_eq!(dfngen.accepted_poly.len(), 2);
        assert_eq!(dfngen.intpts.len(), 1);
    }

    #[test]
    fn test_pieces_of_one_fracture_have_no_intersection() {
        let dfngen = insert_user_polygons(&[(LEFT_PIECE, Some(0)), (RIGHT_PIECE, Some(0))]);

        assert_eq!(dfngen.accepted_poly.len(), 2);
        assert!(dfngen.intpts.is_empty());
    }

    #[test]
    fn test_pieces_of_one_fracture_are_in_one_cluster() {
        let dfngen = insert_user_polygons(&[(LEFT_PIECE, Some(0)), (RIGHT_PIECE, Some(0))]);

        assert_eq!(
            dfngen.accepted_poly[0].group_num,
            dfngen.accepted_poly[1].group_num
        );
    }

    #[test]
    fn test_neighbor_polygon_of_another_fracture_is_rejected() {
        let dfngen = insert_user_polygons(&[(LEFT_PIECE, Some(0)), (RIGHT_PIECE, Some(1))]);

        assert_eq!(dfngen.pstats.rejected_user_fracture, [2]);
    }

    #[test]
    fn test_fracture_can_cross_the_shared_edge_of_two_pieces() {
        let dfngen = insert_user_polygons(&[
            (LEFT_PIECE, Some(0)),
            (RIGHT_PIECE, Some(0)),
            (CROSSING_SQUARE, None),
        ]);

        assert_eq!(dfngen.accepted_poly.len(), 3);
        assert_eq!(dfngen.intpts.len(), 2);
    }

    #[test]
    fn test_fracture_cannot_cross_two_fractures_with_a_gap_less_than_h() {
        let dfngen = insert_user_polygons(&[
            (LEFT_PIECE, Some(0)),
            (RIGHT_PIECE_AFTER_GAP, Some(1)),
            (CROSSING_SQUARE, None),
        ]);

        assert_eq!(dfngen.pstats.rejected_user_fracture, [3]);
    }

    #[test]
    fn test_rejected_user_fracture_is_not_kept_when_the_option_is_off() {
        let dfngen = insert_two_equal_rectangles(false);

        assert_eq!(dfngen.accepted_poly.len(), 1);
        assert!(dfngen.rejected_poly.is_empty());
    }
}
