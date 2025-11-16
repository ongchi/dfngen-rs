/// Spatial Indexing Module - Octree Implementation
///
/// This module provides efficient spatial indexing for fracture intersection checking.
/// Instead of checking every new fracture against all existing fractures (O(n²)),
/// we use an octree to quickly find candidate fractures that might intersect,
/// reducing average case complexity to O(n log n).

use parry3d_f64::bounding_volume::{Aabb, BoundingVolume};
use parry3d_f64::na::{Point3, Vector3};

/// Helper trait to extend parry3d's Aabb with utility methods
pub trait AabbExt {
    /// Create Aabb from the Poly bounding_box format: [xmin, xmax, ymin, ymax, zmin, zmax]
    fn from_poly_bbox(bbox: &[f64; 6]) -> Self;
}

impl AabbExt for Aabb {
    fn from_poly_bbox(bbox: &[f64; 6]) -> Self {
        let min = Point3::new(bbox[0], bbox[2], bbox[4]);
        let max = Point3::new(bbox[1], bbox[3], bbox[5]);
        Aabb::new(min, max)
    }
}

/// Octree Node
///
/// Stores either leaf data (indices of fractures) or 8 child nodes.
/// A node represents a cubic region of 3D space divided into 8 equal octants.
#[derive(Clone)]
struct OctreeNode {
    /// Bounding box of this node
    bounds: Aabb,
    /// Indices of fractures stored at this node (if leaf)
    fracture_indices: Vec<usize>,
    /// Child nodes (8 octants: x-/x+, y-/y+, z-/z+)
    children: Option<Box<[OctreeNode; 8]>>,
    /// Maximum depth to prevent infinite subdivision
    depth: usize,
}

impl OctreeNode {
    /// Create a new leaf node
    fn new_leaf(bounds: Aabb) -> Self {
        Self {
            bounds,
            fracture_indices: Vec::new(),
            children: None,
            depth: 0,
        }
    }

    /// Create a new internal node with 8 children
    fn subdivide(&mut self) {
        if self.depth >= MAX_DEPTH {
            return;
        }

        let center = self.bounds.center();
        let _half_size = self.bounds.half_extents();

        let mut children: Vec<OctreeNode> = Vec::with_capacity(8);

        // Create 8 child nodes for each octant
        for i in 0..8 {
            let min_x = if (i & 1) == 0 {
                self.bounds.mins.x
            } else {
                center.x
            };
            let max_x = if (i & 1) == 0 {
                center.x
            } else {
                self.bounds.maxs.x
            };

            let min_y = if (i & 2) == 0 {
                self.bounds.mins.y
            } else {
                center.y
            };
            let max_y = if (i & 2) == 0 {
                center.y
            } else {
                self.bounds.maxs.y
            };

            let min_z = if (i & 4) == 0 {
                self.bounds.mins.z
            } else {
                center.z
            };
            let max_z = if (i & 4) == 0 {
                center.z
            } else {
                self.bounds.maxs.z
            };

            let child_bounds = Aabb::new(
                Point3::new(min_x, min_y, min_z),
                Point3::new(max_x, max_y, max_z),
            );

            let mut child = OctreeNode::new_leaf(child_bounds);
            child.depth = self.depth + 1;
            children.push(child);
        }

        self.children = Some(Box::new([
            children[0].clone(),
            children[1].clone(),
            children[2].clone(),
            children[3].clone(),
            children[4].clone(),
            children[5].clone(),
            children[6].clone(),
            children[7].clone(),
        ]));
    }

    /// Find which octant a point belongs to (relative to node center)
    #[inline]
    #[allow(dead_code)]
    fn get_octant_index(&self, point: &Vector3<f64>) -> usize {
        let center = self.bounds.center();
        let mut index = 0;
        if point.x >= center.x {
            index |= 1;
        }
        if point.y >= center.y {
            index |= 2;
        }
        if point.z >= center.z {
            index |= 4;
        }
        index
    }

    /// Insert a fracture index into this node, subdividing if necessary
    fn insert(&mut self, fracture_idx: usize, bbox: &Aabb) {
        // If Aabb doesn't overlap this node, don't insert
        if !self.bounds.intersects(bbox) {
            return;
        }

        // If this is a leaf node
        if self.children.is_none() {
            self.fracture_indices.push(fracture_idx);

            // Subdivide if we have too many fractures at this node
            if self.fracture_indices.len() > MAX_FRACTURES_PER_NODE && self.depth < MAX_DEPTH {
                self.subdivide();

                // Move fractures to children
                let indices = std::mem::take(&mut self.fracture_indices);
                if let Some(children) = &mut self.children {
                    for idx in indices {
                        for child in children.iter_mut() {
                            child.insert(idx, bbox);
                        }
                    }
                }
            }
        } else if let Some(children) = &mut self.children {
            // If internal node, insert into all overlapping children
            for child in children.iter_mut() {
                if child.bounds.intersects(bbox) {
                    child.insert(fracture_idx, bbox);
                }
            }
        }
    }

    /// Query all fractures that might overlap with the given Aabb
    fn query(&self, bbox: &Aabb, results: &mut Vec<usize>) {
        // If this node doesn't overlap the query bbox, skip
        if !self.bounds.intersects(bbox) {
            return;
        }

        // Add fractures from this node
        for &idx in &self.fracture_indices {
            results.push(idx);
        }

        // Recursively query children
        if let Some(children) = &self.children {
            for child in children.iter() {
                child.query(bbox, results);
            }
        }
    }
}

// Parameters for octree tuning
const MAX_DEPTH: usize = 8; // Maximum depth of octree (8 levels = 256³ cells in worst case)
const MAX_FRACTURES_PER_NODE: usize = 16; // Subdivide node if it has more fractures

/// Octree Spatial Index
///
/// Efficiently stores and queries fracture bounding boxes in 3D space.
/// Supports dynamic insertion of fractures and spatial queries.
pub struct Octree {
    root: OctreeNode,
}

#[allow(dead_code)]
impl Octree {
    /// Create a new octree covering the entire domain
    pub fn new(domain_size: &Vector3<f64>) -> Self {
        let bounds = Aabb::new(Point3::origin(), Point3::new(domain_size.x, domain_size.y, domain_size.z));
        Self {
            root: OctreeNode::new_leaf(bounds),
        }
    }

    /// Insert a fracture at the given index with its bounding box
    pub fn insert(&mut self, fracture_idx: usize, bbox: &Aabb) {
        self.root.insert(fracture_idx, bbox);
    }

    /// Query all candidate fractures that might overlap with the given Aabb
    ///
    /// This returns indices of fractures whose bounding boxes overlap the query Aabb.
    /// Further geometric tests are still needed to confirm actual intersections.
    pub fn query_overlapping(&self, bbox: &Aabb) -> Vec<usize> {
        let mut results = Vec::new();
        self.root.query(bbox, &mut results);
        results
    }

    /// Query candidates using the Poly bounding_box format: [xmin, xmax, ymin, ymax, zmin, zmax]
    pub fn query_overlapping_poly_bbox(&self, bbox: &[f64; 6]) -> Vec<usize> {
        self.query_overlapping(&Aabb::from_poly_bbox(bbox))
    }

    /// Get statistics about the octree (for debugging/optimization)
    pub fn stats(&self) -> OctreeStats {
        let mut stats = OctreeStats::default();
        self.root.collect_stats(&mut stats);
        stats
    }

    fn collect_stats(&self, stats: &mut OctreeStats) {
        self.root.collect_stats(stats);
    }
}

/// Statistics about octree structure
#[derive(Default, Debug)]
#[allow(dead_code)]
pub struct OctreeStats {
    pub total_nodes: usize,
    pub leaf_nodes: usize,
    pub total_fractures: usize,
    pub max_fractures_per_node: usize,
    pub max_depth_used: usize,
}

#[allow(dead_code)]
impl OctreeNode {
    fn collect_stats(&self, stats: &mut OctreeStats) {
        stats.total_nodes += 1;
        stats.total_fractures += self.fracture_indices.len();
        stats.max_fractures_per_node =
            stats.max_fractures_per_node.max(self.fracture_indices.len());
        stats.max_depth_used = stats.max_depth_used.max(self.depth);

        if let Some(children) = &self.children {
            for child in children.iter() {
                child.collect_stats(stats);
            }
        } else {
            stats.leaf_nodes += 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_aabb_overlap() {
        let aabb1 = Aabb::new(Point3::new(0.0, 0.0, 0.0), Point3::new(10.0, 10.0, 10.0));
        let aabb2 = Aabb::new(Point3::new(5.0, 5.0, 5.0), Point3::new(15.0, 15.0, 15.0));
        let aabb3 = Aabb::new(Point3::new(20.0, 20.0, 20.0), Point3::new(30.0, 30.0, 30.0));

        assert!(aabb1.intersects(&aabb2));
        assert!(aabb2.intersects(&aabb1));
        assert!(!aabb1.intersects(&aabb3));
    }

    #[test]
    fn test_octree_insertion_and_query() {
        let mut octree = Octree::new(&Vector3::new(100.0, 100.0, 100.0));

        // Insert some fractures
        let bbox1 = Aabb::new(Point3::new(0.0, 0.0, 0.0), Point3::new(10.0, 10.0, 10.0));
        let bbox2 = Aabb::new(Point3::new(50.0, 50.0, 50.0), Point3::new(60.0, 60.0, 60.0));
        let bbox3 = Aabb::new(Point3::new(5.0, 5.0, 5.0), Point3::new(15.0, 15.0, 15.0));

        octree.insert(0, &bbox1);
        octree.insert(1, &bbox2);
        octree.insert(2, &bbox3);

        // Query should return candidates that overlap
        let candidates = octree.query_overlapping(&Aabb::new(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(12.0, 12.0, 12.0),
        ));

        assert!(candidates.contains(&0), "Should find fracture 0");
        assert!(!candidates.contains(&1), "Should not find fracture 1");
        assert!(candidates.contains(&2), "Should find fracture 2");
    }
}
