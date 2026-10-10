use parry3d_f64::na::{center, distance, Point3};

use crate::structures::{DFNGen, PolyOptions};

use super::poly::Poly;

/// Family number of the user-defined fractures.
pub const USER_FRACTURE_FAMILY: isize = -1;

/// Inserts the user-defined fractures in the order of their file.
/// Returns the number of them.
pub fn insert_user_defined_fractures(
    fractures: &[UserFracture],
    poly_opts: &PolyOptions,
    dfngen: &mut DFNGen,
) -> usize {
    for (i, fracture) in fractures.iter().enumerate() {
        dfngen.insert_poly(fracture.to_poly(), i + 1, fracture.region_id, poly_opts);
    }

    fractures.len()
}

/// A fracture which the user defines by its vertices, in the file of the
/// user-defined fractures.
#[derive(Debug)]
pub struct UserFracture {
    pub vertices: Vec<[f64; 3]>,

    /// Region of the domain mesh which the fracture must intersect. 0 is any
    /// region, greater than 0 is the position of the region in
    /// `DomainMesh::regions` plus 1.
    pub region_id: usize,
}

impl UserFracture {
    pub fn to_poly(&self) -> Poly {
        let mut poly = Poly {
            family_num: USER_FRACTURE_FAMILY,
            ..Default::default()
        };
        let nodes = &self.vertices;

        // Get a normal vector
        // Vector from fist node to node accross middle of polygon
        let pt_idx_12 = nodes.len() / 2;
        let p1 = Point3::from_slice(&nodes[0]);
        let p2 = Point3::from_slice(&nodes[1]);
        let p_12 = Point3::from_slice(&nodes[pt_idx_12]);

        let v2 = p2 - p1;

        // For 3 nodes, the node accross the middle is the second node.
        // Then the third node gives the normal vector.
        let p_normal = Point3::from_slice(&nodes[pt_idx_12.max(2)]);
        poly.normal = v2.cross(&(p_normal - p1)).normalize();

        // Estimate translation (middle node)
        // Use midpoint between 1st and and half way around polygon
        // Note: For polygons defined by coordinates, the coordinates
        // themselves provide the translation. We need to estimate the center
        // of the polygon and init. the translation array
        poly.translation = center(&p1, &p_12).coords;

        // Estimate radius
        // across middle if even number of nodes
        poly.xradius = 0.5 * v2.magnitude();

        // Get idx for node 1/4 around polygon
        let pt_idx_14 = nodes.len() / 4;
        // Get idx for node 3/4 around polygon
        // across middle close to perpendicular to xradius magnitude calculation
        let pt_idx_34 = pt_idx_14 + pt_idx_12;

        let p_14 = Point3::from_slice(&nodes[pt_idx_14]);
        let p_34 = Point3::from_slice(&nodes[pt_idx_34]);

        poly.yradius = 0.5 * distance(&p_14, &p_34);

        poly.aspect_ratio = poly.yradius / poly.xradius;
        poly.vertices = nodes.iter().map(|node| Point3::from_slice(node)).collect();

        poly
    }
}

#[cfg(test)]
mod tests {
    use parry3d_f64::na::Vector3;

    use super::*;

    #[test]
    fn test_triangle_by_vertices_has_its_normal() {
        let triangle = UserFracture {
            vertices: vec![[0., 0., 0.], [1., 0., 0.], [1., 1., 0.]],
            region_id: 0,
        };

        assert_eq!(triangle.to_poly().normal, Vector3::new(0., 0., 1.));
    }
}
