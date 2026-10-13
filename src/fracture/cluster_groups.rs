use tracing::info;

use super::poly::Poly;

/// One cluster group
struct Group {
    /// Number of fractures in the group
    size: usize,
    /// False after the group was merged into another group
    is_valid: bool,
}

/// Fractures which went into one group together
///
/// A merge of two groups changes only `group` of the lists of the merged
/// group. Thus more than one list can have the same group.
struct MemberList {
    group: usize,
    /// Indices of the accepted fractures
    fractures: Vec<usize>,
}

/// Clusters of connected fractures
///
/// Each accepted fracture is in one group. A fracture with no connection gets
/// a new group. A fracture which connects two or more groups merges them. The
/// number of a group is its position in the order in which the groups were made.
#[derive(Default)]
pub struct ClusterGroups {
    /// All groups which were made, by group number. A merged group stays, as not valid.
    groups: Vec<Group>,
    member_lists: Vec<MemberList>,
}

impl ClusterGroups {
    /// Number of groups which were made. The merged groups are included.
    pub fn group_count(&self) -> usize {
        self.groups.len()
    }

    /// Removes all groups.
    pub fn clear(&mut self) {
        self.groups.clear();
        self.member_lists.clear();
    }

    /// Puts a fracture with no connection into a new group.
    ///
    /// # Arguments
    ///
    /// * `fracture_index` - Index of the fracture in the accepted fractures
    ///
    /// # Returns
    ///
    /// The number of the new group
    pub fn add_isolated(&mut self, fracture_index: usize) -> usize {
        let group = self.groups.len();
        self.groups.push(Group {
            size: 1,
            is_valid: true,
        });
        self.member_lists.push(MemberList {
            group,
            fractures: vec![fracture_index],
        });

        group
    }

    /// Puts a fracture into the group of a fracture which it is connected to,
    /// and merges the other groups which it connects into that group.
    ///
    /// # Arguments
    ///
    /// * `group` - Group of the first fracture which the new fracture is connected to
    /// * `fracture_index` - Index which the new fracture gets in the accepted fractures
    /// * `bridged_groups` - Other groups which the new fracture is connected to
    /// * `accepted_poly` - Accepted fractures. The fractures of the merged groups get `group`.
    pub fn add_connected(
        &mut self,
        group: usize,
        fracture_index: usize,
        bridged_groups: &[usize],
        accepted_poly: &mut [Poly],
    ) {
        self.groups[group].size += 1;
        self.member_lists
            .iter_mut()
            .find(|list| list.group == group)
            .expect("the group of an accepted fracture has a member list")
            .fractures
            .push(fracture_index);

        for &bridged in bridged_groups {
            if self.groups[bridged].is_valid {
                self.groups[group].size += self.groups[bridged].size;
                self.groups[bridged].is_valid = false;
            }

            for list in self
                .member_lists
                .iter_mut()
                .filter(|list| list.group == bridged)
            {
                list.group = group;

                for &member in &list.fractures {
                    accepted_poly[member].group = Some(group);
                }
            }
        }
    }

    /// Selects the fractures which stay in the network.
    ///
    /// # Arguments
    ///
    /// * `keep_isolated_fractures` - If true, all accepted fractures stay
    /// * `keep_only_largest_cluster` - If true, only the group with the most fractures stays
    /// * `accepted_count` - Number of accepted fractures
    ///
    /// # Returns
    ///
    /// Indices of the fractures in the groups of 2 or more fractures, in no
    /// special order, or of all accepted fractures.
    pub fn final_fractures(
        &self,
        keep_isolated_fractures: bool,
        keep_only_largest_cluster: bool,
        accepted_count: usize,
    ) -> Vec<usize> {
        info!("In cluster groups");
        info!("Number of fractures: {}", accepted_count);
        info!("Number of groups: {}", self.group_count());

        if keep_isolated_fractures {
            info!("Number of fractures: {}", accepted_count);
            return (0..accepted_count).collect();
        }

        // A group which is valid was not merged into another group, and a
        // group of more than 1 fracture has an intersection
        let mut connected_groups: Vec<usize> = self
            .groups
            .iter()
            .enumerate()
            .filter(|(_, group)| group.is_valid && group.size > 1)
            .map(|(number, _)| number)
            .collect();

        if keep_only_largest_cluster {
            if let Some(largest) = self.largest_group(&connected_groups) {
                connected_groups = vec![largest];
            }
        }

        connected_groups
            .iter()
            .flat_map(|&group| {
                self.member_lists
                    .iter()
                    .filter(move |list| list.group == group)
                    .flat_map(|list| &list.fractures)
            })
            .copied()
            .collect()
    }

    /// Group with the most fractures of the given groups. The first one wins a tie.
    fn largest_group(&self, groups: &[usize]) -> Option<usize> {
        groups.iter().copied().reduce(|largest, group| {
            if self.groups[group].size > self.groups[largest].size {
                group
            } else {
                largest
            }
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Groups 0, 1 and 2 with the given numbers of fractures
    fn groups_with_sizes(sizes: [usize; 3]) -> ClusterGroups {
        ClusterGroups {
            groups: sizes
                .into_iter()
                .map(|size| Group {
                    size,
                    is_valid: true,
                })
                .collect(),
            member_lists: Vec::new(),
        }
    }

    #[test]
    fn test_largest_group_is_selected_by_size_not_by_group_number() {
        let clusters = groups_with_sizes([2, 50, 7]);

        assert_eq!(clusters.largest_group(&[0, 1, 2]), Some(1));
    }

    #[test]
    fn test_largest_group_looks_only_at_the_given_groups() {
        let clusters = groups_with_sizes([2, 50, 7]);

        assert_eq!(clusters.largest_group(&[0, 2]), Some(2));
    }

    #[test]
    fn test_fracture_which_connects_two_groups_merges_them() {
        let mut clusters = ClusterGroups::default();
        let mut fractures = vec![Poly::default(), Poly::default()];
        for (index, fracture) in fractures.iter_mut().enumerate() {
            fracture.group = Some(clusters.add_isolated(index));
        }

        // Fracture 2 is connected to fracture 0 (group 0) and to fracture 1 (group 1)
        clusters.add_connected(0, 2, &[1], &mut fractures);

        let mut final_fractures = clusters.final_fractures(false, false, 3);
        final_fractures.sort_unstable();
        assert_eq!(final_fractures, [0, 1, 2]);
        assert_eq!(fractures[1].group, Some(0));
    }

    #[test]
    fn test_fracture_with_no_connection_is_not_a_final_fracture() {
        let mut clusters = ClusterGroups::default();
        clusters.add_isolated(0);

        assert!(clusters.final_fractures(false, false, 1).is_empty());
    }
}
