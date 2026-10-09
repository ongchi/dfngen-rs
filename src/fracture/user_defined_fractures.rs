use parry3d_f64::na::{center, distance, Point3, Vector3};

use crate::structures::{DFNGen, PolyOptions};

use super::poly::Poly;

/// Family number of the user-defined fractures.
pub const USER_FRACTURE_FAMILY: isize = -1;

/// A fracture which the user defines in the input file.
#[derive(Debug)]
pub enum UserFracture {
    ByParameters(FractureDef),
    ByVertices(FractureNodesDef),
}

impl UserFracture {
    pub fn to_poly(&self) -> Poly {
        let mut poly = match self {
            Self::ByParameters(fracture) => fracture.to_poly(),
            Self::ByVertices(fracture) => fracture.to_poly(),
        };
        poly.family_num = USER_FRACTURE_FAMILY;
        poly
    }
}

/// Inserts the user-defined fractures in the order of the input file.
/// Returns the number of them.
pub fn insert_user_defined_fractures(
    fractures: &[UserFracture],
    poly_opts: &PolyOptions,
    dfngen: &mut DFNGen,
) -> usize {
    for (i, fracture) in fractures.iter().enumerate() {
        dfngen.insert_poly(fracture.to_poly(), i + 1, poly_opts);
    }

    fractures.len()
}

/// User-defined fracture that is defined by its vertices.
#[derive(Debug)]
pub struct FractureNodesDef {
    pub vertices: Vec<[f64; 3]>,
}

impl FractureNodesDef {
    pub fn to_poly(&self) -> Poly {
        let mut poly = Poly::default();
        let nodes = &self.vertices;

        // Get a normal vector
        // Vector from fist node to node accross middle of polygon
        let pt_idx_12 = nodes.len() / 2;
        let p1 = Point3::from_slice(&nodes[0]);
        let p2 = Point3::from_slice(&nodes[1]);
        let p_12 = Point3::from_slice(&nodes[pt_idx_12]);

        let v1 = p_12 - p1;
        let v2 = p2 - p1;

        poly.normal = v2.cross(&v1).normalize();

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

/// User-defined fracture that is defined by a set of properties.
#[derive(Debug)]
pub struct FractureDef {
    pub eps: f64, // TODO: Move this parameter into global
    pub number_of_points: usize,
    pub radius: f64,
    /// Rotation around the normal vector in radians
    pub beta: f64,
    pub aspect_ratio: f64,
    pub translation: [f64; 3],
    pub normal: [f64; 3],
}

impl FractureDef {
    pub fn to_poly(&self) -> Poly {
        let normal = Vector3::from_row_slice(&self.normal);

        let mut new_poly = Poly::new(self.number_of_points, self.radius, self.aspect_ratio);
        new_poly.translation = Vector3::from_row_slice(&self.translation);

        // Apply 2d rotation matrix, twist around origin
        // Assumes polygon on x-y plane
        // Angle must be in rad
        new_poly.rotation_2d(self.beta);

        // Rotate vertices to uenormal[index] (new normal)
        new_poly.rotation_3d(&normal, self.eps);

        // Save newPoly's new normal vector
        new_poly.normal = normal;

        // Translate newPoly to uetranslation
        new_poly.translate(Vector3::from_row_slice(&self.translation));

        new_poly
    }
}
