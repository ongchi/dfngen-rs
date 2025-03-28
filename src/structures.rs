use parry3d_f64::na::{Point3, Vector3};

#[derive(Clone, Default)]
/// The Poly structre is used to create and store fracrures/polygons.
pub struct Poly {
    /// Number of nodes/vertices in the polygon.
    pub number_of_nodes: isize,

    /// Contains the index of the 'shapeFamilies' array in main() from which the fracture was
    /// created from. If the polygon was created from user defined input, it will be marked
    /// (-2 for user-rect, and -1 for user ell.
    /// The stochastic shape family numbers start at 0 for the first family, and increase by
    /// 1 for each addition family. Families are in the same order which they are defined in the
    /// 'famProb' variable in the input file, starting with ellipse families, then rectangular families.
    pub family_num: isize,

    /// Fracture cluster number which the fracture belongs to. This variable is used to keep
    /// track of fracture connectivity. When a fracture first intersects another fracture,
    /// it inherits its cluster group number (groupNum). If a fracture does not intersect
    /// any other fractures, it is given a new and unique cluster group number. When a
    /// fracture bridges two different clusters, all clusters are merged to be have the
    /// cluster group number of the first intersecting fracture.
    pub group_num: usize,

    /// Polygon area (Not calculated until after DFN generation has completed).
    pub area: f64,

    /// X-radius before fracture-domain truncation. In the case of rectangles, radius is
    /// 1/2 the width of the polygon.
    /// X-radius is equal to the value generated from randum distributions or given by the user
    /// in the case of constant distributions and user-defined fractures.
    pub xradius: f64,

    /// Y-radius before fracture-domain truncation. In the case of rectangles, radius is 1/2
    /// the width of the polygon.
    /// Y-radius is equal to x-radius * aspect ratio (yradius = aspectRatio * xradius).
    pub yradius: f64,

    /// Aspect ratio of polygon before fracture-domain truncation. Must be value greater than zero.
    pub aspect_ratio: f64,

    /// Translation of polygon. This variable is set while building the polygon.
    pub translation: Vector3<f64>,

    /// Polygon normal. This variable is set while building the polygon.
    pub normal: Vector3<f64>,

    /// The bounding box of the polygon. Set with createBoundingBox().
    /// Index Key:
    /// [0]: x minimum, [1]: x maximum
    /// [2]: y minimum, [3]: y maximum
    /// [4]: z minimum, [5]: z maximum
    pub bounding_box: [f64; 6],

    /// Double array for which hold the polygon's vertices. Vertices are stored in a 1-D array.
    /// e.g For n number of vertices, array will be: {x1, y1, z1, x2, y2, z2, ... , xn, yn, zn}
    pub vertices: Vec<f64>,

    /// The faces array contains flags (true/false) which denote which sides, if any, the polygon is touching.
    /// True (not zero) - Polygon is touching a domain boundary.
    /// False (0) - Polygon is not touching a domain boundary.
    /// Index Key:
    /// [0]: -x face, [1]: +x face
    /// [2]: -y face, [3]: +y face
    /// [4]: -z face, [5]: +z face
    /// Touching boundary faces of the fractures group, not neccesarily the faces of the fracture
    pub faces: [bool; 6],

    /// When writing intersection points and poly vertices to output, we need them to be entirely on the x-y plane.
    /// XYPlane is used for error checking. During intersection rotations, the polygon is rotated to
    /// the x-y plane and it's vertices are changed within the Poly structure. Errors will occur if
    /// the same polygon was to be rotated again because although the vertices are now on the x-y plane, the
    /// normal is still the normal of the polygon in its 3D space. This variable prevents this from happening. */
    pub xyplane: bool,

    /// True if the polygon has been truncated, false otherwise. This variable is used for re-translating polygons.
    /// If the polygon has been truncated, it must be re-built. Otherwise, it can simply be given a new translation. */
    pub truncated: bool,

    /// List of indices to the permanent intersection array ('intPts' in main()) which belong to this polygon.
    pub intersection_index: Vec<usize>,
}

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

/// Shape is used to hold varibales for all types of stochastic shapes. During getInput(),
/// all stochastic families for both recaangles and ellipses are parsed from the user input
/// and are placed in a Shape structure array.
#[derive(Default)]
pub struct Shape {
    /// 0 = ellipse, 1 = rectangle
    pub shape_family: usize,

    /// 1: Lognormal, 2: truncated power-law, 3: exponential, 4: constant
    pub distribution_type: u8,

    /// Number of vertices used to create the polygon. Ellipse families only.
    pub num_points: usize,

    /// Array of thetas to build poly from, initialized while reading input and building shape structures
    pub theta_list: Vec<f64>,

    /// Current index to the radii list 'radiiList'.
    pub radii_idx: usize,

    /// Initial list of fracture/polygon radii, sorted largest to smallest.
    pub radii_list: Vec<f64>,

    /// Layer the family belongs to. 0 is entire domain, greater than 0 is a layer.
    /// e.g. 2 would be the second layer listed in the input file under "layers:".
    pub layer: usize,

    /// Region the family belongs to. 0 is entire domain, greater than 0 is a region.
    /// e.g. 2 would be the second region listed in the input file under "regions:".
    pub region: usize,

    /// Aspect ratio for family.
    pub aspect_ratio: f64,

    /// Target p32 (fracture intensity) for the family when using p32 program-stopping option.
    pub p32_target: f64,

    /// Current P32 value for this family.
    pub current_p32: f64,

    /// True = degrees, False = radians. This variable is set while readin the user's input file.
    /// After reading in the users input file, any input in degrees
    /// will be  changed to radians.
    pub angle_option: bool,

    /// 'betaOption' is the rotation about the polygon's  normal vector
    /// True - User Specified Rotation
    /// False - Uniform Distribution
    pub beta_distribution: bool,

    /// 'beta' is the rotation, or twist, around z normal before 3d rotation in radians
    /// or degrees depending on 'angleOption'.
    pub beta: f64,

    /// If orientationOption = 0 (Spherical coordinates)
    /// This is the angle the normal vector makes with the z-axis (theta)
    /// If  orientationOption = 1
    /// This is the trend of Rectangle fracture orientation.
    pub angle_one: f64,

    /// If orientationOption = 0 (Spherical coordinates)
    /// This is the angle the normal vector makes with the z-axis (phi)
    /// If  orientationOption = 1
    /// This is the trend of Rectangle fracture orientation.
    pub angle_two: f64,

    /// Parameter for fisher distributions. The
    /// bigger, the more similar (less diverging) are the
    /// rectangular familiy's normal vectors.
    pub kappa: f64,

    /**************** Distribution Variables *********************/
    /*************************************************************/
    /// Value between 0 and 1. Input to distrubution which will generate the user's defined
    /// minimum value from the distribution. Currently used only for exponential
    /// distrubution in the Distributions class. Value is set during Distributions
    /// constructor.
    pub min_dist_input: f64,

    /// Value between 0 and 1. Input to distrubution which will generate the user's defined
    /// maximum value from the distribution. Currently used only for exponential
    /// distrubution in the Distributions class. Value is set during Distributions
    /// constructor.
    pub max_dist_input: f64,

    /// Exponential distribution option. Mean value for exponential distribution.
    pub exp_mean: f64,

    /// Exponential distribution option. Lambda value for exponential distibution. This
    /// value is set by using 1/'expMean' while reading the user's input file.
    pub exp_lambda: f64,

    /// Exponential distribution option. User's chosen minimum value for the distribution.
    /// The distribution will never return a value smaller than this.
    pub exp_min: f64,

    /// Exponential distribution option. User's chosen minimum value for the distribution.
    /// The distribution will never return a value larger than this.
    pub exp_max: f64,

    /// Log-normal distribution option. Mean of underlying normal distribution from which
    /// the log-normal distribution is created.
    pub mean: f64,

    /// Log-normal distribution option. Standard deviation of the underlying
    /// normal distribution from which the log-normal distribution is created.
    pub sd: f64,

    /// Log-normal distribution option. User's chosen minimum value for the distribution.
    /// The distribution will never return a value smaller than this.
    pub log_min: f64,
    /// Log-normal distribution option. User's chosen maximum value for the distribution.
    /// The distribution will never return a value smaller than this.
    pub log_max: f64,

    /// Constant distribution. Constant radii size for all fractures in the family.
    pub const_radi: f64,

    /// Truncated power-law option. Minimum radius for power-law distribution.
    pub min: f64,

    /// Truncated power-law option. Maximum radius for power-law distribution.
    pub max: f64,

    /// Alpha. Used in truncated power-law distribution calculations.
    pub alpha: f64,
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
