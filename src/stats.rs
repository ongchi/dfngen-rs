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

/// Program and DFN statistics structure. Keeps various statistics about the DFN being generated.
pub struct Stats {
    /// Counters for the number of polygons/fractures accepted by each stochastic
    /// family. Elements in this array are in the same order as the stochastic family
    /// families array 'shapeFamilies' in main(). e.g. The counter for the second family
    /// in the shapeFamilies array is the second element in this array.
    pub accepted_per_family: Vec<usize>,

    /// Counters for the number of polygons/fractures rejected by each stochastic
    /// family. Elements in this array are in the same order as the stochastic family
    /// families array 'shapeFamilies' in main(). e.g. The counter for the second family
    /// in the shapeFamilies array is the second element in this array.
    pub rejected_per_family: Vec<usize>,

    /// Number of fractures estimated for each stochastic family (see dryRun()).
    /// Elements in this array are in the same order as the stochastic family
    /// families array 'shapeFamilies' in main(). e.g. The estimated number of
    /// fractures for the second family in the shapeFamilies array is the second
    /// element in this array.
    pub expected_per_family: Vec<usize>,

    /// Total number of accepted polygons/fractures for the DFN. This variable is
    /// updated as the DFN is generated.
    pub accepted_count: usize,

    /// Total number of rejected polygons/fractures for the DFN. This variable is
    /// updated as the DFN is generated.
    pub rejected_count: usize,

    /// Total number of polygon/fracture re-translations. This variable is updated
    /// as the DFN is generated.
    pub retranslated_count: u32,

    /// Total area of fractures before isolated fracture removal. Variable is set in main() after
    /// DFN generation has completed.
    pub area_before_removal: f64,

    /// Total area of fractures after isolated fracture removal. Variable is set in main() after
    /// DFN generation has completed.
    pub area_after_removal: f64,

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

    /// Rejected User Fractures
    pub rejected_user_fractures: Vec<usize>,
}

impl Default for Stats {
    /// All counters are zero.
    fn default() -> Self {
        Self {
            accepted_per_family: Vec::new(),
            rejected_per_family: Vec::new(),
            expected_per_family: Vec::new(),
            accepted_count: 0,
            rejected_count: 0,
            retranslated_count: 0,
            area_before_removal: 0.,
            area_after_removal: 0.,
            rejection_reasons: RejectionReasons::default(),
            intersections_shortened: 0,
            original_length: 0.,
            discarded_length: 0.,
            // Counter of the rejects before the first accepted fracture
            rejects_per_attempt: vec![0],
            rejected_user_fractures: Vec::new(),
        }
    }
}
