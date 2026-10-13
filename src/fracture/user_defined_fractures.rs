use std::collections::HashSet;

use parry3d_f64::na::{center, distance, Point3};

use crate::structures::{DFNGen, PolyOptions};

use super::poly::Poly;

/// Family number of the user-defined fractures.
pub const USER_FRACTURE_FAMILY: isize = -1;

/// Inserts the user-defined fractures in the order of their file.
/// Returns the number of them.
///
/// A fracture of more than one polygon (`fracture_id`) is in the domain if one
/// of its pieces is in the domain. Then no piece of it is rejected as outside.
pub fn insert_user_defined_fractures(
    fractures: &[UserFracture],
    poly_opts: &PolyOptions,
    dfngen: &mut DFNGen,
) -> usize {
    let fracture_ids_in_domain: HashSet<u64> = fractures
        .iter()
        .filter(|fracture| is_piece_in_domain(fracture, poly_opts))
        .filter_map(|fracture| fracture.fracture_id)
        .collect();

    for (i, fracture) in fractures.iter().enumerate() {
        let domain_check = match fracture.fracture_id {
            Some(id) if fracture_ids_in_domain.contains(&id) => DomainCheck::FractureIsInDomain,
            _ => DomainCheck::Regions(&fracture.region_ids),
        };
        dfngen.insert_poly(fracture.to_poly(), i + 1, domain_check, poly_opts);
    }

    fractures.len()
}

/// True if the polygon has a part in the box domain and intersects a cell of
/// one of its regions of the domain mesh.
fn is_piece_in_domain(fracture: &UserFracture, opts: &PolyOptions) -> bool {
    let mut poly = fracture.to_poly();

    !poly.check_box_domain(opts.h, opts.eps, &opts.domain_size)
        && opts
            .domain_mesh
            .intersects_polygon_in_any_region(&poly.vertices, &fracture.region_ids)
}

/// Domain check of a user-defined fracture at its insertion.
#[derive(Clone, Copy)]
pub enum DomainCheck<'a> {
    /// The polygon must intersect a cell of one of these regions of the
    /// domain mesh. No region means a cell of any region.
    Regions(&'a [usize]),

    /// The polygon is a piece of a fracture which has a piece in the domain.
    /// It is not rejected as outside the domain.
    FractureIsInDomain,
}

/// A fracture which the user defines by its vertices, in the file of the
/// user-defined fractures.
#[derive(Debug)]
pub struct UserFracture {
    pub vertices: Vec<[f64; 3]>,

    /// Regions of the domain mesh: the fracture must intersect a cell of one
    /// of them. Each is the position of the region in `DomainMesh::regions`
    /// plus 1. Empty is a cell of any region.
    pub region_ids: Vec<usize>,

    /// Value of the `fracture_id` array of the file. Polygons with the same
    /// value are pieces of one fracture. `None` if the file has no such array.
    pub fracture_id: Option<u64>,
}

impl UserFracture {
    pub fn to_poly(&self) -> Poly {
        let mut poly = Poly {
            family_num: USER_FRACTURE_FAMILY,
            user_fracture_id: self.fracture_id,
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
            region_ids: Vec::new(),
            fracture_id: None,
        };

        assert_eq!(triangle.to_poly().normal, Vector3::new(0., 0., 1.));
    }
}
