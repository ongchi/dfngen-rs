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

impl RejectionReasons {
    /// Each counter of `self` minus that of `old` plus that of `new`
    fn replace_part(&self, old: &Self, new: &Self) -> Self {
        let replace = |field: fn(&Self) -> usize| field(self) - field(old) + field(new);
        Self {
            short_intersection: replace(|reasons| reasons.short_intersection),
            close_to_node: replace(|reasons| reasons.close_to_node),
            close_to_edge: replace(|reasons| reasons.close_to_edge),
            close_point_to_edge: replace(|reasons| reasons.close_point_to_edge),
            outside: replace(|reasons| reasons.outside),
            triple: replace(|reasons| reasons.triple),
            inter_close_to_inter: replace(|reasons| reasons.inter_close_to_inter),
            boundary_line: replace(|reasons| reasons.boundary_line),
        }
    }
}

impl Stats {
    /// The statistics with other user-defined fractures: `old_user` and
    /// `new_user` are the statistics after the insertion of the old and the
    /// new user-defined fractures, before the stochastic fractures
    pub fn with_user_part(&self, old_user: &Self, new_user: &Self) -> Self {
        // The counter of the rejects at the end of the user part has also
        // the rejects of the first stochastic attempts
        let old_end = old_user.rejects_per_attempt.len() - 1;
        let new_end = new_user.rejects_per_attempt.len() - 1;
        let mut rejects_per_attempt = new_user.rejects_per_attempt[..new_end].to_vec();
        rejects_per_attempt.push(
            self.rejects_per_attempt[old_end] - old_user.rejects_per_attempt[old_end]
                + new_user.rejects_per_attempt[new_end],
        );
        rejects_per_attempt.extend_from_slice(&self.rejects_per_attempt[old_end + 1..]);

        Self {
            accepted_count: self.accepted_count - old_user.accepted_count + new_user.accepted_count,
            rejected_count: self.rejected_count - old_user.rejected_count + new_user.rejected_count,
            rejection_reasons: self
                .rejection_reasons
                .replace_part(&old_user.rejection_reasons, &new_user.rejection_reasons),
            intersections_shortened: self.intersections_shortened
                - old_user.intersections_shortened
                + new_user.intersections_shortened,
            original_length: self.original_length - old_user.original_length
                + new_user.original_length,
            discarded_length: self.discarded_length - old_user.discarded_length
                + new_user.discarded_length,
            rejects_per_attempt,
            rejected_user_fractures: new_user.rejected_user_fractures.clone(),
            rejected_user_fracture_count: new_user.rejected_user_fracture_count,
            ..self.clone()
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

    /// Statistics after `accepted` user-defined fractures with `rejects`
    /// rejects before each, and `rejected` user-defined fractures
    fn user_stats(accepted: usize, rejects: usize, rejected: &[usize]) -> Stats {
        let mut stats = Stats {
            accepted_count: accepted,
            rejected_count: accepted * rejects + rejected.len(),
            rejects_per_attempt: vec![rejects; accepted + 1],
            rejected_user_fractures: rejected.to_vec(),
            ..Stats::default()
        };
        stats.rejects_per_attempt[accepted] = rejected.len();
        stats.rejection_reasons.outside = stats.rejected_count;
        stats
    }

    #[test]
    fn test_with_user_part_replaces_the_counters_of_the_user_fractures() {
        let old_user = user_stats(2, 1, &[3]);
        // The stochastic part: 3 accepted, with 2, 0 and 4 rejects before them
        let mut total = old_user.clone();
        total.accepted_count += 3;
        total.rejected_count += 6;
        total.rejection_reasons.outside += 6;
        total.rejects_per_attempt[2] += 2;
        total.rejects_per_attempt.extend([0, 4, 0]);
        let new_user = user_stats(1, 2, &[]);

        let stats = total.with_user_part(&old_user, &new_user);

        let expected = Stats {
            accepted_count: 4,
            rejected_count: 8,
            rejects_per_attempt: vec![2, 2, 0, 4, 0],
            rejection_reasons: RejectionReasons {
                outside: 8,
                ..RejectionReasons::default()
            },
            ..Stats::default()
        };
        assert_eq!(stats, expected);
    }
}
