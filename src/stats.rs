/// Number of rejects for each reason
#[derive(Default)]
pub struct RejectionReasons {
    /// Rejections due to intersection of length less than min_feature_size.
    pub short_intersection: usize,
    /// Rejections due to an intersection too close to a vertex of a polygon.
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

/// Statistics of the DFN which is generated
pub struct Stats {
    /// Number of accepted fractures of each stochastic family, in the order of the families
    pub accepted_per_family: Vec<usize>,

    /// Number of rejected fractures of each stochastic family, in the order of the families
    pub rejected_per_family: Vec<usize>,

    /// Number of fractures estimated for each stochastic family (the length of its radii list),
    /// in the order of the families
    pub expected_per_family: Vec<usize>,

    /// Number of accepted fractures
    pub accepted_count: usize,

    /// Number of rejected fractures
    pub rejected_count: usize,

    /// Number of times a rejected fracture was moved to a new position
    pub retranslated_count: u32,

    /// Total area of the fractures before the removal of the isolated fractures
    pub area_before_removal: f64,

    /// Total area of the fractures after the removal of the isolated fractures
    pub area_after_removal: f64,

    /// Number of rejects for each reason
    pub rejection_reasons: RejectionReasons,

    /// Number of intersections of accepted fractures which FRAM shortened
    pub intersections_shortened: u32,

    /// Total length of the intersections before FRAM shortened them
    pub original_length: f64,

    /// Total length which FRAM removed from the intersections. The final
    /// intersection length is `original_length - discarded_length`.
    pub discarded_length: f64,

    /// Number of rejects before each accepted fracture. Element 9 is the number of
    /// attempts which were rejected before the 10th fracture was accepted. An
    /// attempt is a new fracture or a new position of a rejected fracture.
    pub rejects_per_attempt: Vec<usize>,

    /// Numbers of the rejected user-defined fractures in their file, from 1
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
