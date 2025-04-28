use std::{fs::File, io::Read};

use itertools::zip_eq;
use parry3d_f64::na::{distance, Point3, Vector3};
use text_io::read;

use crate::error::DfngenError;
use crate::io::{
    input::ExternalFractureFiles,
    read_input_functions::{search_var, ReadFromTextFile, UserFractureReader},
};
use crate::structures::{DFNGen, PolyOptions};

use super::poly::Poly;

pub fn insert_user_defined_fractures(
    rect_first: bool,
    files: &ExternalFractureFiles,
    poly_opts: &PolyOptions,
    dfngen: &mut DFNGen,
) -> usize {
    let mut num_user_defined_fractures = 0;

    // User Polygons are always inserted first
    if let Some(ref path) = files.user_poly_by_coord_file {
        let poly_def = FractureNodesDef::from_poly_file(path);
        for (i, frac) in poly_def.into_iter().enumerate() {
            dfngen.insert_poly(frac.try_into().unwrap(), i + 1, -3, poly_opts);
            num_user_defined_fractures += 1;
        }
    }

    let mut rect_polys = Vec::new();
    let mut ell_polys = Vec::new();

    if let Some(ref path) = files.user_rect_file {
        let rect_def = FractureDef::from_file(path, poly_opts.eps, false);
        for frac in rect_def.into_iter() {
            rect_polys.push(frac.try_into().unwrap());
            num_user_defined_fractures += 1;
        }
    }

    if let Some(ref path) = files.user_rect_by_coord_file {
        let rect_def = FractureNodesDef::from_rect_file(path);
        for frac in rect_def.into_iter() {
            rect_polys.push(frac.try_into().unwrap());
            num_user_defined_fractures += 1;
        }
    }

    if let Some(ref path) = files.user_ell_file {
        let ell_def = FractureDef::from_file(path, poly_opts.eps, true);
        for frac in ell_def.into_iter() {
            ell_polys.push(frac.try_into().unwrap());
            num_user_defined_fractures += 1;
        }
    }

    if let Some(ref path) = files.user_ell_by_coord_file {
        let ell_def = FractureNodesDef::from_ell_file(path);
        for frac in ell_def.into_iter() {
            ell_polys.push(frac.try_into().unwrap());
            num_user_defined_fractures += 1;
        }
    }

    if rect_first {
        for (i, poly) in rect_polys.into_iter().enumerate() {
            dfngen.insert_poly(poly, i + 1, -2, poly_opts);
        }
        for (i, poly) in ell_polys.into_iter().enumerate() {
            dfngen.insert_poly(poly, i + 1, -1, poly_opts);
        }
    } else {
        for (i, poly) in ell_polys.into_iter().enumerate() {
            dfngen.insert_poly(poly, i + 1, -1, poly_opts);
        }
        for (i, poly) in rect_polys.into_iter().enumerate() {
            dfngen.insert_poly(poly, i + 1, -2, poly_opts);
        }
    }

    num_user_defined_fractures
}

/// User-defined fractures that are defined by a set of nodes.
#[derive(Debug)]
pub enum FractureNodesDef {
    Ellipse(Vec<[f64; 3]>),
    Rectangle(Vec<[f64; 3]>),
    Polygon(Vec<[f64; 3]>),
}

impl FractureNodesDef {
    pub fn new_ell(nodes: Vec<[f64; 3]>) -> Result<Self, DfngenError> {
        Ok(Self::Ellipse(nodes))
    }

    pub fn new_rect(nodes: Vec<[f64; 3]>) -> Result<Self, DfngenError> {
        if nodes.len() != 4 {
            Err(DfngenError::ValueError(
                "User defined rectangle must have 4 nodes".to_string(),
            ))
        } else {
            Ok(Self::Rectangle(nodes))
        }
    }

    pub fn new_poly(nodes: Vec<[f64; 3]>) -> Result<Self, DfngenError> {
        Ok(Self::Polygon(nodes))
    }

    /// Reads a file containing user-defined fractures defined by ellipse vertices.
    pub fn from_ell_file(path: &str) -> Vec<Self> {
        let mut file = File::open(path).unwrap();

        search_var(&mut file, "nEllipses:");
        let mut n_frac: usize = 0;
        n_frac.read_from_text(&mut file);

        search_var(&mut file, "nNodes:");
        let mut num_points: usize = 0;
        num_points.read_from_text(&mut file);

        search_var(&mut file, "Coordinates:");
        let mut bytes = file.bytes().map(|ch| ch.unwrap());

        let mut fractures = Vec::new();

        for _ in 0..n_frac {
            let mut nodes = Vec::new();
            for _ in 0..num_points {
                let mut tmp = [0., 0., 0.];
                for val in &mut tmp {
                    *val = read!("{}", bytes);
                }
                nodes.push(tmp);
            }
            fractures.push(FractureNodesDef::new_ell(nodes).unwrap());
        }

        fractures
    }

    /// Reads a file containing user-defined fractures defined by rectangle vertices.
    pub fn from_rect_file(path: &str) -> Vec<Self> {
        let mut file = File::open(path).unwrap();

        search_var(&mut file, "nRectangles:");
        let mut n_frac: usize = 0;
        n_frac.read_from_text(&mut file);

        search_var(&mut file, "Coordinates:");

        let mut bytes = file.bytes().map(|ch| ch.unwrap());

        let mut fractures = Vec::new();

        for _ in 0..n_frac {
            let mut nodes = Vec::new();
            for _ in 0..4 {
                let mut tmp = [0., 0., 0.];
                for val in &mut tmp {
                    *val = read!("{}", bytes);
                }
                nodes.push(tmp);
            }
            fractures.push(FractureNodesDef::new_rect(nodes).unwrap());
        }

        fractures
    }

    /// Reads a file containing user-defined fractures defined by polygon vertices.
    pub fn from_poly_file(path: &str) -> Vec<Self> {
        let mut file = File::open(path).unwrap();

        search_var(&mut file, "nPolygons:");
        let mut n_frac: usize = 0;
        n_frac.read_from_text(&mut file);

        let mut num_points: usize = 0;
        num_points.read_from_text(&mut file);

        let mut fractures = Vec::new();

        for _ in 0..n_frac {
            let mut nodes = Vec::new();
            for _ in 0..num_points {
                let mut tmp = [0., 0., 0.];
                tmp.read_from_text(&mut file);
                nodes.push(tmp);
            }
            fractures.push(FractureNodesDef::new_poly(nodes).unwrap());
        }

        fractures
    }
}

impl TryInto<Poly> for FractureNodesDef {
    type Error = DfngenError;

    fn try_into(self) -> Result<Poly, Self::Error> {
        let mut poly = Poly::default();

        let nodes = match self {
            Self::Ellipse(ref nodes) => {
                poly.family_num = -1;
                nodes
            }
            Self::Rectangle(ref nodes) => {
                poly.family_num = -2;
                nodes
            }
            Self::Polygon(ref nodes) => {
                poly.family_num = -3;
                nodes
            }
        };

        poly.number_of_nodes = nodes.len() as isize;

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
        poly.translation = 0.5 * Vector3::new(p1.x + p_12.x, p1.y + p_12.y, p1.z + p_12.z);

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
        poly.vertices = nodes.iter().flat_map(|&node| node.to_vec()).collect();

        Ok(poly)
    }
}

/// User-defined fractures that are defined by a set of properties.
#[derive(Debug)]
pub struct FractureDef {
    eps: f64, // TODO: Move this parameter into global
    num_nodes: usize,
    radius: f64,
    beta: f64,
    aspect_ratio: f64,
    translation: [f64; 3],
    normal: [f64; 3],
}

impl FractureDef {
    pub fn new(
        eps: f64,
        num_nodes: usize,
        radius: f64,
        beta: f64,
        aspect_ratio: f64,
        translation: [f64; 3],
        normal: [f64; 3],
    ) -> Self {
        Self {
            eps,
            num_nodes,
            radius,
            beta,
            aspect_ratio,
            translation,
            normal,
        }
    }

    pub fn from_file(path: &str, eps: f64, is_ell: bool) -> Vec<Self> {
        let mut frac_reader = if is_ell {
            UserFractureReader::new(path, "nUserEll:")
        } else {
            UserFractureReader::new(path, "nUserRect:")
        };

        let mut radii: Vec<f64> = Vec::new();
        let mut beta: Vec<f64> = Vec::new();
        let mut aspect: Vec<f64> = Vec::new();
        let mut translation: Vec<[f64; 3]> = Vec::new();

        let mut orientation_option: u8 = 0;

        let mut normal: Vec<[f64; 3]> = Vec::new();
        let mut trend_plunge: Vec<[f64; 2]> = Vec::new();
        let mut dip_strike: Vec<[f64; 2]> = Vec::new();

        let mut num_points: Vec<usize> = Vec::new();

        frac_reader.read_vec("Radii:", &mut radii);

        let mut angle_option = String::new();
        frac_reader.read_value("AngleOption:", &mut angle_option);
        let angle_convertion_factor = match angle_option.as_str() {
            "degree" => std::f64::consts::PI / 180.,
            "radian" => 1.,
            _ => panic!("Invalid angle option"),
        };

        frac_reader.read_vec("Beta:", &mut beta);
        beta.iter_mut().for_each(|v| *v *= angle_convertion_factor);

        frac_reader.read_vec("Aspect_Ratio:", &mut aspect);
        frac_reader.read_vec_arr3("Translation:", &mut translation);

        frac_reader.read_value("userOrientationOption:", &mut orientation_option);

        if orientation_option == 0 {
            frac_reader.read_vec_arr3("Normal:", &mut normal);
        } else if orientation_option == 1 {
            frac_reader.read_vec_arr2("Trend_Plunge:", &mut trend_plunge);

            // Convert Trend and Plunge into normal vectors
            for [trend, plunge] in &trend_plunge {
                let trend = trend * angle_convertion_factor;
                let plunge = plunge * angle_convertion_factor;
                normal.push([
                    trend.cos() * plunge.cos(),
                    trend.sin() * plunge.cos(),
                    plunge.sin(),
                ]);
            }
        } else if orientation_option == 2 {
            frac_reader.read_vec_arr2("Dip_Strike:", &mut dip_strike);

            // Convert Dip and Strike into normal vectors
            for [dip, strike] in &dip_strike {
                let dip = dip * angle_convertion_factor;
                let strike = strike * angle_convertion_factor;
                normal.push([
                    dip.sin() * strike.sin(),
                    -dip.sin() * strike.cos(),
                    dip.cos(),
                ]);
            }
        }

        if is_ell {
            frac_reader.read_vec("Number_of_Vertices:", &mut num_points);
        } else {
            num_points = vec![0; frac_reader.n_frac];
        }

        let mut fractures = Vec::new();

        for (((((num_nodes, radius), beta), aspect_ratio), translation), normal) in zip_eq(
            zip_eq(
                zip_eq(zip_eq(zip_eq(num_points, radii), beta), aspect),
                translation,
            ),
            normal,
        ) {
            fractures.push(FractureDef::new(
                eps,
                if is_ell { num_nodes } else { 0 },
                radius,
                beta,
                aspect_ratio,
                translation,
                normal,
            ))
        }

        fractures
    }
}

impl TryInto<Poly> for FractureDef {
    type Error = DfngenError;

    fn try_into(self) -> Result<Poly, Self::Error> {
        let is_rect = self.num_nodes == 0;

        let family_num = if is_rect { -2 } else { -1 };
        let normal = Vector3::from_row_slice(&self.normal);

        let mut new_poly = if is_rect {
            Poly::new_rect(self.radius, self.aspect_ratio)
        } else {
            Poly::new_ell(self.num_nodes, self.radius, self.aspect_ratio)
        };
        new_poly.family_num = family_num;
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

        Ok(new_poly)
    }
}
