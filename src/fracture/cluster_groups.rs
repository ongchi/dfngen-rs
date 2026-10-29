use std::collections::HashSet;

use tracing::{info, warn};

use super::poly::Poly;

/// Which clusters stay in the network (input section `[clusters]`)
#[derive(Clone, Debug)]
pub struct ClusterOptions {
    /// Keep all accepted fractures; the other options then have no effect
    pub keep_isolated_fractures: bool,
    /// Keep only the cluster with the most fractures
    pub keep_only_largest: bool,
    /// Keep only the clusters which touch a boundary surface of each group
    /// (`keep_only_connected_to`; empty = all clusters). Boundary N is
    /// `DomainMesh::boundary_surfaces()[N - 1]`.
    pub connected_surface_groups: Vec<Vec<usize>>,
}

/// The fractures which stay in the network (`ClusterGroups::final_fractures`)
pub struct ClusterSelection {
    /// Indices of the accepted fractures, in no special order
    pub fractures: Vec<usize>,
    /// True if `connected_surface_groups` found clusters but none which
    /// touches a surface of each group: then all clusters stay
    pub is_without_connecting_cluster: bool,
}

/// The boundary surfaces which the accepted fractures touch: the surfaces of
/// the segments of their boundary lines
pub struct SurfaceContacts<'a> {
    pub accepted_polys: &'a [Poly],
}

impl SurfaceContacts<'_> {
    /// True if the fractures together touch a surface of each group
    fn touch_each_group<'b>(
        &self,
        fractures: impl Iterator<Item = &'b usize>,
        groups: &[Vec<usize>],
    ) -> bool {
        let touched: HashSet<usize> = fractures
            .flat_map(|&fracture| &self.accepted_polys[fracture].boundary_lines)
            .flat_map(|line| &line.segments)
            .flat_map(|segment| &segment.side.surface_ids)
            .copied()
            .collect();
        groups
            .iter()
            .all(|group| group.iter().any(|surface_id| touched.contains(surface_id)))
    }
}

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
    /// * `accepted_polys` - Accepted fractures. The fractures of the merged groups get `group`.
    pub fn add_connected(
        &mut self,
        group: usize,
        fracture_index: usize,
        bridged_groups: &[usize],
        accepted_polys: &mut [Poly],
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
                    accepted_polys[member].group = Some(group);
                }
            }
        }
    }

    /// Selects the fractures which stay in the network. Port of `getCluster`
    /// (`clusterGroups.cpp`) of dfnWorks; its `boundaryFaces` are the groups
    /// of boundary surfaces here (each face in its own group gives the rule of
    /// dfnWorks).
    ///
    /// # Arguments
    ///
    /// * `options` - Which clusters stay
    /// * `accepted_count` - Number of accepted fractures
    /// * `contacts` - The boundary surfaces of the accepted fractures, for
    ///   `connected_surface_groups`
    ///
    /// # Returns
    ///
    /// The fractures in the groups of 2 or more fractures, or all accepted
    /// fractures. With `connected_surface_groups` and clusters of which none
    /// touches a surface of each group, all clusters stay (as dfnWorks), with
    /// a warning and `is_without_connecting_cluster`.
    pub fn final_fractures(
        &self,
        options: &ClusterOptions,
        accepted_count: usize,
        contacts: &SurfaceContacts,
    ) -> ClusterSelection {
        info!("In cluster groups");
        info!("Number of fractures: {}", accepted_count);
        info!("Number of groups: {}", self.group_count());

        if options.keep_isolated_fractures {
            info!("Number of fractures: {}", accepted_count);
            return ClusterSelection {
                fractures: (0..accepted_count).collect(),
                is_without_connecting_cluster: false,
            };
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

        let mut is_without_connecting_cluster = false;
        let surface_groups = &options.connected_surface_groups;
        if !surface_groups.is_empty() {
            let connecting_groups: Vec<usize> = connected_groups
                .iter()
                .copied()
                .filter(|&group| contacts.touch_each_group(self.members(group), surface_groups))
                .collect();
            info!(
                "Number of clusters which connect the boundary surface groups: {}",
                connecting_groups.len()
            );
            // No cluster at all: the caller reports that no fractures intersect
            if connecting_groups.is_empty() && !connected_groups.is_empty() {
                warn!("No cluster connects the boundary surface groups: all clusters stay");
                is_without_connecting_cluster = true;
            } else {
                connected_groups = connecting_groups;
            }
        }

        if options.keep_only_largest {
            if let Some(largest) = self.largest_group(&connected_groups) {
                connected_groups = vec![largest];
            }
        }

        ClusterSelection {
            fractures: connected_groups
                .iter()
                .flat_map(|&group| self.members(group))
                .copied()
                .collect(),
            is_without_connecting_cluster,
        }
    }

    /// Indices of the accepted fractures of a group
    fn members(&self, group: usize) -> impl Iterator<Item = &usize> {
        self.member_lists
            .iter()
            .filter(move |list| list.group == group)
            .flat_map(|list| &list.fractures)
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
    use parry3d_f64::na::Point3;

    use super::*;
    use crate::region_boundaries::{BoundaryLine, BoundarySegment, OtherSide, SegmentSide};

    const KEEP_CLUSTERS: ClusterOptions = ClusterOptions {
        keep_isolated_fractures: false,
        keep_only_largest: false,
        connected_surface_groups: Vec::new(),
    };

    /// No boundary surface: for the options which do not look at the surfaces
    const NO_CONTACTS: SurfaceContacts = SurfaceContacts {
        accepted_polys: &[],
    };

    /// Keep the clusters which touch a surface of each group
    fn keep_clusters_connected_to(groups: &[&[usize]]) -> ClusterOptions {
        ClusterOptions {
            connected_surface_groups: groups.iter().map(|group| group.to_vec()).collect(),
            ..KEEP_CLUSTERS
        }
    }

    /// Surface 1 in one group, surface 2 in the other (as in-flow and out-flow)
    fn keep_clusters_connected_to_one_and_two() -> ClusterOptions {
        keep_clusters_connected_to(&[&[1], &[2]])
    }

    /// A fracture with one boundary segment on the given boundary surfaces
    fn fracture_on_surfaces(surface_ids: &[usize]) -> Poly {
        Poly {
            boundary_lines: vec![BoundaryLine {
                segments: vec![BoundarySegment {
                    start: Point3::origin(),
                    end: Point3::new(1., 0., 0.),
                    side: SegmentSide {
                        other_side: OtherSide::MeshSurface,
                        surface_ids: surface_ids.to_vec(),
                    },
                }],
                is_closed: false,
            }],
            ..Poly::default()
        }
    }

    /// Clusters of the given fractures: each cluster is a list of fractures,
    /// each fracture the list of the boundary surfaces which it touches.
    /// Cluster N is group N.
    fn clusters_of(fracture_surfaces: &[&[&[usize]]]) -> (ClusterGroups, Vec<Poly>) {
        let mut clusters = ClusterGroups::default();
        let mut fractures = Vec::new();
        for cluster in fracture_surfaces {
            let group = clusters.add_isolated(fractures.len());
            fractures.push(fracture_on_surfaces(cluster[0]));
            for surfaces in &cluster[1..] {
                clusters.add_connected(group, fractures.len(), &[], &mut fractures);
                fractures.push(fracture_on_surfaces(surfaces));
            }
        }
        (clusters, fractures)
    }

    /// The selection of clusters, with sorted fractures
    fn select(fracture_surfaces: &[&[&[usize]]], options: &ClusterOptions) -> ClusterSelection {
        let (clusters, fractures) = clusters_of(fracture_surfaces);
        let contacts = SurfaceContacts {
            accepted_polys: &fractures,
        };
        let mut selection = clusters.final_fractures(options, fractures.len(), &contacts);
        selection.fractures.sort_unstable();
        selection
    }

    #[test]
    fn test_cluster_which_touches_one_of_two_groups_is_removed() {
        // Cluster 0: fractures 0, 1 on surfaces 1 and 2; cluster 1: fractures 2, 3 on surface 1
        let selection = select(
            &[&[&[1], &[2]], &[&[1], &[]]],
            &keep_clusters_connected_to_one_and_two(),
        );

        assert_eq!(selection.fractures, [0, 1]);
    }

    #[test]
    fn test_surface_in_no_group_is_not_needed() {
        // Cluster 0 touches surfaces 1 and 2, not surface 3
        let selection = select(&[&[&[1], &[2]]], &keep_clusters_connected_to_one_and_two());

        assert_eq!(selection.fractures, [0, 1]);
    }

    #[test]
    fn test_one_surface_of_a_group_is_enough() {
        // Group 1: surfaces 1 and 3; cluster 0 touches surfaces 3 and 2
        let selection = select(
            &[&[&[3], &[2]], &[&[1], &[]]],
            &keep_clusters_connected_to(&[&[1, 3], &[2]]),
        );

        assert_eq!(selection.fractures, [0, 1]);
    }

    #[test]
    fn test_a_connecting_cluster_clears_the_flag() {
        let selection = select(
            &[&[&[1], &[2]], &[&[1], &[]]],
            &keep_clusters_connected_to_one_and_two(),
        );

        assert!(!selection.is_without_connecting_cluster);
    }

    #[test]
    fn test_all_clusters_stay_when_no_cluster_connects_the_groups() {
        // Cluster 0 touches only surface 1, cluster 1 only surface 2
        let selection = select(
            &[&[&[1], &[]], &[&[], &[2]]],
            &keep_clusters_connected_to_one_and_two(),
        );

        assert_eq!(selection.fractures, [0, 1, 2, 3]);
    }

    #[test]
    fn test_no_connecting_cluster_sets_the_flag() {
        let selection = select(
            &[&[&[1], &[]], &[&[], &[2]]],
            &keep_clusters_connected_to_one_and_two(),
        );

        assert!(selection.is_without_connecting_cluster);
    }

    #[test]
    fn test_no_cluster_at_all_does_not_set_the_flag() {
        // One fracture with no intersection: the caller reports NoIntersectingFractures
        let selection = select(&[&[&[1, 2]]], &keep_clusters_connected_to_one_and_two());

        assert!(!selection.is_without_connecting_cluster);
    }

    #[test]
    fn test_keep_only_largest_selects_among_the_connecting_clusters() {
        // Cluster 0 (3 fractures) touches only surface 1; cluster 1 (2 fractures) both
        let options = ClusterOptions {
            keep_only_largest: true,
            ..keep_clusters_connected_to_one_and_two()
        };

        let selection = select(&[&[&[1], &[], &[1]], &[&[1, 2], &[]]], &options);

        assert_eq!(selection.fractures, [3, 4]);
    }

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

        let mut final_fractures = clusters
            .final_fractures(&KEEP_CLUSTERS, 3, &NO_CONTACTS)
            .fractures;
        final_fractures.sort_unstable();
        assert_eq!(final_fractures, [0, 1, 2]);
        assert_eq!(fractures[1].group, Some(0));
    }

    #[test]
    fn test_fracture_with_no_connection_is_not_a_final_fracture() {
        let mut clusters = ClusterGroups::default();
        clusters.add_isolated(0);

        assert!(clusters
            .final_fractures(&KEEP_CLUSTERS, 1, &NO_CONTACTS)
            .fractures
            .is_empty());
    }
}
