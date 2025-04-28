use std::fmt::Formatter;
use std::rc::Rc;
use std::{cell::RefCell, fmt::Display};

use parry3d_f64::na::{Point3, Vector3};
use rand::Rng;
use rand_mt::Mt64;
use tracing::{error, info, warn};

use crate::cg::{intersection_checking, is_in_boundary};
use crate::distribution::{TruncExp, TruncLogNormal, TruncPowerLaw};
use crate::error::DfngenError;
use crate::fracture::insert_shape::print_reject_reason;
use crate::fracture::poly::Poly;

pub struct RejectedUserFracture {
    pub id: usize,
    pub user_fracture_type: i32,
}

impl RejectedUserFracture {
    pub fn new(id: usize, user_fracture_type: i32) -> Self {
        Self {
            id,
            user_fracture_type,
        }
    }
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
    /// Rejections due to intersection of length less than h.
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

    /// Total number of fractures that have been truncated against the domain.
    pub truncated: usize,

    /// Total area of fractures before isolated fracture removal. Variable is set in main() after
    /// DFN generation has completed.
    pub area_before_removal: f64,

    /// Total area of fractures after isolated fracture removal. Variable is set in main() after
    /// DFN generation has completed.
    pub area_after_removal: f64,

    // Total volume of fractures before isolated fracture removal. Variable is set in
    // main() after DFN generation has completed.
    // vol_before_removal: f64,

    // Total volume of fractures after isolated fracture removal. Variable is set in
    // main() after DFN generation has completed.
    // vol_after_removal: f64,
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

    /// Total number of intersection points in DFN after isolated fracture removal.
    /// This variable is used as a counter in writeIntersections() when generating
    /// output files. It is used to determine the number of nodes lagrit will see
    /// as duplicates and remove. Used in error checking.
    pub intersection_node_count: usize,

    /// Total number of triple intersection points in DFN after isolated fracture removal.
    /// This variable is used as a counter in writeIntersections() when generating
    /// output files. It is used to determine the number of nodes lagrit will see
    /// as duplicates and remove. Used in error checking.
    pub triple_node_count: usize,

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
    pub rejected_user_fracture: Vec<RejectedUserFracture>,
}

#[derive(Clone, Copy, Debug)]
pub enum Shape {
    Ellipse(u8),
    Rectangle,
}

impl Default for Shape {
    fn default() -> Self {
        Shape::Ellipse(6)
    }
}

impl Display for Shape {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Shape::Ellipse(_) => write!(f, "Ellipse"),
            Shape::Rectangle => write!(f, "Rectangular"),
        }
    }
}

impl Shape {
    pub fn number_of_nodes(&self) -> u8 {
        match self {
            Shape::Ellipse(n) => *n,
            Shape::Rectangle => 4,
        }
    }
}

pub enum RadiusFunction {
    LogNormal { mu: f64, sigma: f64 },
    TruncatedPowerLaw { alpha: f64 },
    Exponential { lambda: f64 },
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

    pub fn new_log_normal(mu: f64, sigma: f64, min: f64, max: f64) -> Self {
        Self {
            function: RadiusFunction::LogNormal { mu, sigma },
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
            RadiusFunction::LogNormal { mu, sigma } => {
                let distr = TruncLogNormal::new(self.min, f64::INFINITY, mu, sigma)?;
                generator.borrow_mut().sample(&distr)?
            }
            RadiusFunction::TruncatedPowerLaw { alpha } => {
                let distr = TruncPowerLaw::new(self.min, self.max, alpha);
                generator.borrow_mut().sample(&distr)
            }
            RadiusFunction::Exponential { lambda } => {
                let distr = TruncExp::new(self.min, f64::INFINITY, lambda)?;
                generator.borrow_mut().sample(&distr)?
            }
            RadiusFunction::Constant(c) => c,
        };

        Ok(radius)
    }
}

impl Display for RadiusDistribution {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match &self.function {
            RadiusFunction::LogNormal { mu, sigma } => {
                writeln!(f, "Distribution: Lognormal")?;
                writeln!(f, "Mean: {}", mu)?;
                writeln!(f, "Standard Deviation: {}", sigma)?;
                writeln!(f, "Minimum Radius (m): {}", self.min)?;
                writeln!(f, "Maximum Radius (m): {}", self.max)
            }
            RadiusFunction::TruncatedPowerLaw { alpha } => {
                writeln!(f, "Distribution: Truncated Power-Law")?;
                writeln!(f, "Alpha: {}", alpha)?;
                writeln!(f, "Minimum Radius (m): {}", self.min)?;
                writeln!(f, "Maximum Radius (m): {}", self.max)
            }
            RadiusFunction::Exponential { lambda } => {
                writeln!(f, "Distribution: Exponential")?;
                writeln!(f, "Mean: {}", 1. / lambda)?;
                writeln!(f, "Lambda: {}", lambda)?;
                writeln!(f, "Minimum Radius (m): {}", self.min)?;
                writeln!(f, "Maximum Radius (m): {}", self.max)
            }
            RadiusFunction::Constant(value) => {
                writeln!(f, "Distribution: Constant")?;
                writeln!(f, "Radius (m): {}", value)
            }
        }
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
            ..Default::default()
        }
    }
}

pub struct DFNGen {
    // Statistics structure:
    pub pstats: Stats,

    // Vector to store accepted polygons/fractures
    pub accepted_poly: Vec<Poly>,

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
            intpts: Vec::new(),
            triple_points: Vec::new(),
        }
    }

    pub fn insert_poly(&mut self, poly: Poly, poly_id: usize, family_id: i32, opts: &PolyOptions) {
        let mut new_poly = poly;

        if new_poly.domain_truncation(opts.h, opts.eps, &opts.domain_size) {
            // Poly completely outside domain
            self.pstats.rejection_reasons.outside += 1;
            self.pstats.rejected_poly_count += 1;
            self.pstats
                .rejected_user_fracture
                .push(RejectedUserFracture::new(poly_id, family_id));
            warn!(
                "({}, {}): rejected for being outside the defined domain.",
                family_id, poly_id
            );
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
            if new_poly.truncated {
                self.pstats.truncated += 1;
            }

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
            self.pstats
                .rejected_user_fracture
                .push(RejectedUserFracture::new(poly_id, family_id));

            print_reject_reason(reject_code, &new_poly);
        }
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
        let mut final_poly_list = Vec::new();
        // Clear GroupData
        self.pstats.group_data.clear();
        // Clear FractGroup
        self.pstats.fract_group.clear();
        // Clear Triple Points
        self.triple_points.clear();
        // Clear IntPoints
        self.intpts.clear();
        // Re-init nextGroupNum
        self.pstats.next_group_num = 1;

        for poly in self.accepted_poly.iter_mut() {
            if poly.xradius < min_size {
                poly.vertices.clear();
                continue;
            }

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

            // IF POLY ACCEPTED:
            if reject_code == 0 {
                // Intersections are ok
                // SAVING POLYGON (intersection and triple points saved witchin intersectionChecking())
                final_poly_list.push(new_poly); // SAVE newPoly to accepted polys list
            } else {
                // Poly rejected
                warn!("Error rebuilding dfn, previously accepted fracture was rejected during DFN rebuild.");
            }
        }

        info!("Rebuilding DFN complete.");
        self.accepted_poly.clear();
        self.accepted_poly.extend(final_poly_list);
    }

    /// Remove Fractures Outside 2D polygon domain
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
    /// # Arguments
    ///
    /// * `num_of_domain_vertices` - Number of vertices in the domain polygon
    /// * `domain_vertices` -
    pub fn remove_fractures_beyond_domains(
        &mut self,
        opts: &PolyOptions,
        num_of_domain_vertices: usize,
        domain_vertices: &[Point3<f64>],
    ) {
        let mut final_poly_list: Vec<Poly> = Vec::new();

        // Clear GroupData
        self.pstats.group_data.clear();
        // Clear FractGroup
        self.pstats.fract_group.clear();
        // Clear Triple Points
        self.triple_points.clear();
        // Clear IntPoints
        self.intpts.clear();
        // Re-init nextGroupNum
        self.pstats.next_group_num = 1;

        for i in 0..self.accepted_poly.len() {
            let x = self.accepted_poly[i].translation[0];
            let y = self.accepted_poly[i].translation[1];

            // cout << "fracture " << i + 1 << " center " << x << "," << y << endl;
            if !is_in_boundary(num_of_domain_vertices, domain_vertices, x, y) {
                self.accepted_poly.clear();
                continue;
            }

            let new_poly = &mut self.accepted_poly[i];
            new_poly.group_num = 0; // Reset cluster group number
            new_poly.intersection_index.clear(); // Remove ref to old intersections
                                                 // Find line of intersection and FRAM check
            let reject_code = if opts.disable_fram {
                0
            } else {
                intersection_checking(
                    opts,
                    new_poly,
                    &mut final_poly_list,
                    &mut self.intpts,
                    &mut self.pstats,
                    &mut self.triple_points,
                )
            };

            // IF POLY ACCEPTED:
            if reject_code == 0 {
                // Intersections are ok
                // SAVING POLYGON (intersection and triple points saved witchin intersectionChecking())
                final_poly_list.push(new_poly.clone()); // SAVE newPoly to accepted polys list
            } else {
                // Poly rejected
                error!("Error rebuilding dfn, previously accepted fracture was rejected during DFN rebuild.");
            }
        }

        info!("Rebuilding DFN complete.");
        self.accepted_poly.clear();
        self.accepted_poly.extend(final_poly_list);
    }
}

pub struct PolyOptions {
    pub h: f64,
    pub eps: f64,
    pub domain_size: Vector3<f64>,
    pub r_fram: bool,
    pub disable_fram: bool,
    pub triple_intersections: bool,
}
