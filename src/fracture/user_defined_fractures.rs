use parry3d_f64::na::{center, distance, Point3, Vector3};

use crate::structures::{DFNGen, PolyOptions};

use super::poly::Poly;

/// Family number of user-defined ellipses.
pub const USER_ELLIPSE_FAMILY: isize = -1;
/// Family number of user-defined rectangles.
pub const USER_RECTANGLE_FAMILY: isize = -2;
/// Family number of user-defined polygons.
pub const USER_POLYGON_FAMILY: isize = -3;

/// A fracture which the user defines in the input file.
#[derive(Debug)]
pub enum UserFracture {
    ByParameters(FractureDef),
    ByVertices(FractureNodesDef),
}

impl UserFracture {
    pub fn to_poly(&self) -> Poly {
        match self {
            Self::ByParameters(fracture) => fracture.to_poly(),
            Self::ByVertices(fracture) => fracture.to_poly(),
        }
    }

    pub fn family_num(&self) -> isize {
        match self {
            Self::ByParameters(fracture) => fracture.family_num(),
            Self::ByVertices(fracture) => fracture.family_num(),
        }
    }
}

/// Inserts the user-defined fractures. Returns the number of them.
///
/// Polygons are inserted first, then rectangles and ellipses in the order
/// that `rect_first` gives. Each group keeps the order of the input file.
pub fn insert_user_defined_fractures(
    rect_first: bool,
    fractures: &[UserFracture],
    poly_opts: &PolyOptions,
    dfngen: &mut DFNGen,
) -> usize {
    let insertion_order = if rect_first {
        [
            USER_POLYGON_FAMILY,
            USER_RECTANGLE_FAMILY,
            USER_ELLIPSE_FAMILY,
        ]
    } else {
        [
            USER_POLYGON_FAMILY,
            USER_ELLIPSE_FAMILY,
            USER_RECTANGLE_FAMILY,
        ]
    };

    for family_id in insertion_order {
        let family_fractures = fractures
            .iter()
            .filter(|fracture| fracture.family_num() == family_id);

        for (i, fracture) in family_fractures.enumerate() {
            dfngen.insert_poly(fracture.to_poly(), i + 1, family_id as i32, poly_opts);
        }
    }

    fractures.len()
}

/// User-defined fractures that are defined by a set of nodes.
#[derive(Debug)]
pub enum FractureNodesDef {
    Ellipse(Vec<[f64; 3]>),
    Rectangle(Vec<[f64; 3]>),
    Polygon(Vec<[f64; 3]>),
}

impl FractureNodesDef {
    pub fn family_num(&self) -> isize {
        match self {
            Self::Ellipse(_) => USER_ELLIPSE_FAMILY,
            Self::Rectangle(_) => USER_RECTANGLE_FAMILY,
            Self::Polygon(_) => USER_POLYGON_FAMILY,
        }
    }

    pub fn to_poly(&self) -> Poly {
        let mut poly = Poly {
            family_num: self.family_num(),
            ..Default::default()
        };

        let (Self::Ellipse(nodes) | Self::Rectangle(nodes) | Self::Polygon(nodes)) = self;

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
        match self {
            Self::Ellipse(_) | Self::Polygon(_) => {
                // Get idx for node 1/4 around polygon
                let pt_idx_14 = nodes.len() / 4;
                // Get idx for node 3/4 around polygon
                // across middle close to perpendicular to xradius magnitude calculation
                let pt_idx_34 = pt_idx_14 + pt_idx_12;

                let p_14 = Point3::from_slice(&nodes[pt_idx_14]);
                let p_34 = Point3::from_slice(&nodes[pt_idx_34]);

                poly.yradius = 0.5 * distance(&p_14, &p_34);
            }
            Self::Rectangle(_) => {
                let p4 = Point3::from_slice(&nodes[3]);
                let v3 = p4 - p1;

                // Set radius (x and y radii might be switched based on order of users coordinates)
                poly.yradius = 0.5 * v3.magnitude();
            }
        }

        poly.aspect_ratio = poly.yradius / poly.xradius;
        poly.vertices = nodes.iter().map(|node| Point3::from_slice(node)).collect();

        poly
    }
}

/// User-defined fractures that are defined by a set of properties.
#[derive(Debug)]
pub struct FractureDef {
    pub eps: f64, // TODO: Move this parameter into global
    /// Number of nodes of an ellipse. 0 for a rectangle.
    pub num_nodes: usize,
    pub radius: f64,
    /// Rotation around the normal vector in radians
    pub beta: f64,
    pub aspect_ratio: f64,
    pub translation: [f64; 3],
    pub normal: [f64; 3],
}

impl FractureDef {
    fn is_rect(&self) -> bool {
        self.num_nodes == 0
    }

    pub fn family_num(&self) -> isize {
        if self.is_rect() {
            USER_RECTANGLE_FAMILY
        } else {
            USER_ELLIPSE_FAMILY
        }
    }

    pub fn to_poly(&self) -> Poly {
        let is_rect = self.is_rect();
        let normal = Vector3::from_row_slice(&self.normal);

        let mut new_poly = if is_rect {
            Poly::new_rect(self.radius, self.aspect_ratio)
        } else {
            Poly::new_ell(self.num_nodes, self.radius, self.aspect_ratio)
        };
        new_poly.family_num = self.family_num();
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
