use std::collections::HashSet;
use std::fmt;
use std::ops::{Add, Sub};

use serde::{Deserialize, Serialize};

/// Number of fractures and of their polygons. Polygons with the same
/// `fracture_id` are pieces of one user-defined fracture.
#[derive(Clone, Copy, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct FractureCount {
    pub fractures: usize,
    pub polygons: usize,
}

impl FractureCount {
    /// Counts the polygons with these `fracture_id` values. A polygon with no
    /// `fracture_id` (a stochastic fracture, or a file with no such array) is
    /// one fracture.
    pub fn from_fracture_ids(fracture_ids: impl IntoIterator<Item = Option<u64>>) -> Self {
        let mut ids = HashSet::new();
        let mut count = Self::default();

        for fracture_id in fracture_ids {
            count.polygons += 1;
            if fracture_id.is_none_or(|id| ids.insert(id)) {
                count.fractures += 1;
            }
        }

        count
    }

    /// Count of fractures which have one polygon each
    pub fn of_single_polygons(polygons: usize) -> Self {
        Self {
            fractures: polygons,
            polygons,
        }
    }
}

impl Add for FractureCount {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            fractures: self.fractures + other.fractures,
            polygons: self.polygons + other.polygons,
        }
    }
}

impl Sub for FractureCount {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self {
            fractures: self.fractures - other.fractures,
            polygons: self.polygons - other.polygons,
        }
    }
}

impl fmt::Display for FractureCount {
    fn fmt(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
        write!(formatter, "{} ({} polygons)", self.fractures, self.polygons)
    }
}

/// Number of rejects for each reason
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
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
    /// Rejections due to a boundary line of the region on the fracture:
    /// too short, too near to the edge, or too near to an intersection.
    pub boundary_line: usize,
}

/// Statistics of the DFN which is generated
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
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

    /// Number of rejected user-defined fractures: a fracture with no accepted
    /// polygon
    pub rejected_user_fracture_count: FractureCount,
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
            rejected_user_fracture_count: FractureCount::default(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pieces_with_one_fracture_id_are_one_fracture() {
        let count = FractureCount::from_fracture_ids([Some(7), Some(7), Some(3)]);

        assert_eq!(
            count,
            FractureCount {
                fractures: 2,
                polygons: 3
            }
        );
    }

    #[test]
    fn test_each_polygon_without_fracture_id_is_one_fracture() {
        let count = FractureCount::from_fracture_ids([None, None, Some(7)]);

        assert_eq!(count, FractureCount::of_single_polygons(3));
    }
}
