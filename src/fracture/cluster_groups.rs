use tracing::info;

use crate::{
    error,
    math_functions::or,
    structures::{FractureGroups, GroupData, Poly, Stats},
};

//     This code goes through all the polygons accepted into the domain and returns the indexes to those polygons
//     which match the users boundary faces option.
//     Isoloated fractures and Polygons with no intersections are removed
//
// groupData structure:
// ***********************
//     groupData keeps information as to how many polygons are in the group or cluster. It also keeps
//     track of the boundary faces which the cluster is in contact with. There is also a valid variable (0 or 1).
//     If a polygon connects two different groups, one of the groups is merged together with the other. The one which
//     nolonger exists has it's valid variable set to 0.
//
//     The groupData structure array, or vector, is aligned with the group numbers. If you want to know about group 3 for
//     example, you will look at groupData[ (3-1) ], minus 1 because the array starts at 0 while groups start at 1.
//     If valid = 1 in the structure, that group still exists and has not been merged
//
// fractGroups structure:
// *************************
//     fractGroups[] holds the actual pointers/index numbers to the polygon array along with the group number. Unlike
//     groupData[] described above, fractGroups does not stay aligned to group numbers. To keep from copying, deleting, and
//     re-allocating memory when groups merge together, we simply change the variable groupNum to the new group number.
//
//     Because of this, there will be multiples of the same group numbers, but with different polygons listed. To get all the
//     polygons from a group we must search the fractGroups array for all matching groups
//

/// Get Matching Fracture Clusters
///
/// Uses boundaryFaces input option to get the wanted fracture cluster before writing output files.
///
/// NOTE: 'boundaryFaces' array is a global variable
///
/// # Returns
///
/// Vector of indices to fractures which remained after isolated and non-matching boundary faces fracture removal.
pub fn get_cluster(
    keep_isolated_fractures: bool,
    keep_only_largest_cluster: bool,
    ignore_boundary_faces: bool,
    boundary_faces: &[bool; 6],
    pstats: &Stats,
) -> Vec<usize> {
    // The max number of groups is pstats.nextGroupNum - 1
    let mut matching_groups = Vec::new();
    let mut final_poly_list = Vec::new();

    info!("In cluster groups");
    info!("Number of fractures: {}", pstats.accepted_poly_count);
    info!("Number of groups: {}", pstats.group_data.len());

    if !keep_isolated_fractures {
        // NOTE: (groupNumber-1) = corresponding groupData structures' index of the arary
        //       similarly, the index of groupData + 1 = groupNumber (due to groupNumber starting at 1, array starting at 0)

        // Find all matching groups:
        // Get matching groups from pstats.groupData[]
        for i in 0..pstats.group_data.len() {
            // If the data is valid, meaning that group still exists (hasn't been merged to a new group) and
            // if the group has more than 1 fracture meaning there are intersections and
            // the cluster matches the requirements of the user's boundaryFaces option
            if !ignore_boundary_faces {
                if pstats.group_data[i].valid
                    && pstats.group_data[i].size > 1
                    && faces_match(boundary_faces, &pstats.group_data[i].faces)
                {
                    matching_groups.push(i + 1); //save matching group number
                }
            } else {
                // Get all cluster groups
                if pstats.group_data[i].valid && pstats.group_data[i].size > 1 {
                    matching_groups.push(i + 1); // Save matching group number
                }
            }
        }

        if keep_only_largest_cluster && matching_groups.len() > 1 {
            // If only keeping the largest cluster, find group with largest size
            // Initialize largestGroup
            let mut largest_group = matching_groups[0];

            for gid in matching_groups.iter() {
                // If largest group < current group, the current group is the new largest group
                if largest_group < pstats.group_data[gid - 1].size {
                    largest_group = *gid;
                }
            }

            matching_groups.clear(); // Clear group numbers
            matching_groups.push(largest_group); // Save only the largest group
        }

        // Gather the final polygon numbers/indecies.
        for gid in matching_groups {
            for k in 0..pstats.fract_group.len() {
                if gid == pstats.fract_group[k].group_num {
                    // If the groupNumbers match
                    // copy all poly indecies of matching group
                    for j in 0..pstats.fract_group[k].poly_list.len() {
                        final_poly_list.push(pstats.fract_group[k].poly_list[j]);
                    }
                }
            }
        }
    } else {
        info!("Number of fractures: {}", pstats.accepted_poly_count);

        for i in 0..pstats.accepted_poly_count {
            final_poly_list.push(i);
        }
    }

    final_poly_list
}

/// Test if Boundary Faces Option Match User's Desired Boundary Faces
///
/// Compares the user's faces option to a fracture cluster
///
/// # Arguments
///
/// * `faces_option` - Boundary faces user input option array ('boundaryFaces' in input file)
/// * `faces` - Fracture cluter's boundary faces array.
///     Similar to 'boundaryFaces' input option, but denotes which faces the fracture cluster connects to.
///
/// # Returns
///
/// * `bool`
///     false - If faces meet user's facesOption requirements
///     true - If faces do not meet requirements
fn faces_match(faces_option: &[bool; 6], faces: &[bool; 6]) -> bool {
    for i in 0..6 {
        // If the user specified the face, check if faces match
        if faces_option[i] & !(faces_option[i] & faces[i]) {
            return false;
        }
    }

    true
}

/// Assign New Polygon to a Cluster
///
/// Assigns a new polygon/fracture to a new cluster group number.
/// Assumes 'newPoly' does not intersect with any other fractures.
///
/// # Arguments
///
/// * `newPoly` - New polygon
/// * `pstats` - Program stats structure
/// * `newPolyIndex` - Index of 'newPoly' in the 'acceptedPoly' array (array of all accepted polys)
pub fn assign_group(new_poly: &mut Poly, pstats: &mut Stats, new_poly_index: usize) {
    new_poly.group_num = pstats.next_group_num;
    let mut new_group_data = GroupData::new(); // Keeps fracture cluster data
                                               // Copy newPoly faces info to groupData
    or(&mut new_group_data.faces, &new_poly.faces);
    // Incriment groupData's poly count
    new_group_data.size += 1;
    pstats.group_data.push(new_group_data); // Save boundary face information to permanent location
    let mut new_group = FractureGroups::new(pstats.next_group_num); // 'newPoly' had no intersections, put it in a new group
    pstats.next_group_num += 1;
    new_group.poly_list.push(new_poly_index); // Save index (number) of newPoly
    pstats.fract_group.push(new_group);
}

/// Update Cluster Groups
///
/// Updates fracture cluster group data for the addition of 'newPoly'.
///
/// If 'newPoly' did not bridge any clusters together, 'newPoly' is added to the cluster of
/// the first polygon it intersected with. The cluster groups  boundary connectivity data is
/// updated, clusters fracture count is incremented, and 'newPoly' is added to
/// the list of polygons.
///
/// If 'newPoly' bridged two or more clusters, the function also merges the multiple
/// cluster groups into a single group. 'newPoly' is first added to group of the first
/// fracture it intersected with. Then, any remaining cluster groups 'newPoly' intersected with
/// will have their group number changed to match that of 'newPolys' group number (see struct FractureGroups).
/// Inside the GroupData structure, all polygons will have their group number updated to the new group number.
/// Once polygons and struct FractureGroups have been updated, the FractureGroups corresponding GroupData
/// structure will have its valid bit set to false. (This is more efficient than deleting the GroupData
/// element from its array, which causes memory re-allocation and copying).
///
/// The GroupData structure contains information about it's corresponding FractureGroups structure (see GroupData).
/// To access the cooresponding GroupData structure from a cluster group number, the index to the GroupData array
/// within the Stats structure (pstats variable) will be the the cluster group number subtracted by 1.
/// e.g. If you need to access the GroupData structure for cluster group 10, it will be the variable:
///     pstats.groupData[10-1]
/// To access the corresponding FractureGroups structure for a cluster group number, you must search the FractureGroups
/// structure array and search for ALL matching group numbers (pstats.fractGroup[i].groupNum).
///
/// # Arguments
///
/// * `newPoly` - Reference to new polygon
/// * `acceptedPoly` - Array of all accepted polygons
/// * `encounteredGroups` - Array of group numbers for any other bridged fracture cluster groups
/// * `pstats` - Program statistics structure (contains fracture cluster data)
/// * `newPolyIndex` - Index of 'newPoly' once placed into the array of all accepted polygons
pub fn update_groups(
    new_poly: &Poly,
    accepted_poly: &mut [Poly],
    encountered_groups: &[usize],
    pstats: &mut Stats,
    new_poly_index: usize,
) {
    if encountered_groups.is_empty() {
        // 'newPoly' didn't encounter more than 1 other group of fractures
        // Save newPoly to the group structure
        let mut i = 0;
        // Update boundary face info
        or(
            &mut pstats.group_data[new_poly.group_num - 1].faces,
            &new_poly.faces,
        );
        // NOTE: Groups start at 1, but array starts at zero hence the groupNum-1 for index
        // Incriment the groups size for the new polygon
        pstats.group_data[new_poly.group_num - 1].size += 1;

        // Search for group newPoly belongs to and add newPoly to it's list, add to first matching group found
        while pstats.fract_group[i].group_num != new_poly.group_num && i != pstats.fract_group.len()
        {
            i += 1;
        }

        // Error check
        if i == pstats.fract_group.len() && pstats.fract_group[i].group_num != new_poly.group_num {
            error!("Group not found");
        }

        // Add newPoly to fracture/cluster group
        pstats.fract_group[i].poly_list.push(new_poly_index);
    } else {
        // Multiple cluter groups were encountered.
        // Merge the groups (change/update the group numbers)
        // Search for encountered groups and make them all have the same group numbers
        // First, add newpoly to correct group
        let mut k = 0;

        while new_poly.group_num != pstats.fract_group[k].group_num {
            k += 1;
        }

        // Add new poly to the group
        pstats.fract_group[k].poly_list.push(new_poly_index);
        // Incriment size of group by one for newPoly
        pstats.group_data[new_poly.group_num - 1].size += 1;
        // Add any boundary face data to group fro newPoly
        or(
            &mut pstats.group_data[new_poly.group_num - 1].faces,
            &new_poly.faces,
        );

        // Need merge groups in encounteredGroups list to newPoly.groupNum
        for gid in encountered_groups {
            // Add merging groups size to new group's size
            if pstats.group_data[gid - 1].valid {
                pstats.group_data[new_poly.group_num - 1].size += pstats.group_data[gid - 1].size;
                // Merge any boundary faces data
                let encountered_faces = pstats.group_data[gid - 1].faces;
                or(
                    &mut pstats.group_data[new_poly.group_num - 1].faces,
                    &encountered_faces,
                );
                // Mark merged groups group data as invalid
                pstats.group_data[gid - 1].valid = false;
            }

            for jj in 0..pstats.fract_group.len() {
                if pstats.fract_group[jj].group_num == *gid {
                    // Change group number to the new group number
                    pstats.fract_group[jj].group_num = new_poly.group_num;

                    // Change polys group number inside the group to the new group num
                    for z in 0..pstats.fract_group[jj].poly_list.len() {
                        accepted_poly[pstats.fract_group[jj].poly_list[z]].group_num =
                            new_poly.group_num;
                    }
                }
            }
        }
    }
}
