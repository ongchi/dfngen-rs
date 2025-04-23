use std::{fs::File, io::Read};

use parry3d_f64::na::{distance, Point3, Vector3};
use text_io::read;

use crate::{
    computational_geometry::{apply_rotation2_d, apply_rotation3_d, translate},
    distribution::generating_points::generate_theta,
    io::{
        input::ExternalFractureFiles,
        read_input_functions::{search_var, ReadFromTextFile, UserFractureReader},
    },
    structures::{DFNGen, Poly, PolyOptions},
};

use super::insert_shape::{initialize_ell_vertices, initialize_rect_vertices};

pub fn insert_user_defined_fractures(
    rect_first: bool,
    files: &ExternalFractureFiles,
    poly_opts: &PolyOptions,
    dfngen: &mut DFNGen,
) -> usize {
    let mut num_user_defined_fractures = 0;

    // User Polygons are always inserted first
    if let Some(ref path) = files.user_poly_by_coord_file {
        let poly_def = UserDefinedPolygonByCoord::from_file(path);
        for (i, poly) in poly_def.create_polys().into_iter().enumerate() {
            dfngen.insert_poly(poly, i + 1, -3, poly_opts);
        }
        num_user_defined_fractures += poly_def.n_frac;
    }

    let mut rect_polys = Vec::new();
    let mut ell_polys = Vec::new();

    if let Some(ref path) = files.user_rect_file {
        let rect_def = UserDefinedFractures::from_file(path, false);
        rect_polys.extend(rect_def.create_polys(poly_opts.eps));
        num_user_defined_fractures += rect_def.n_frac;
    }

    if let Some(ref path) = files.user_rect_by_coord_file {
        let rect_def = UserDefinedRectByCoord::from_file(path);
        rect_polys.extend(rect_def.create_polys());
        num_user_defined_fractures += rect_def.n_frac;
    }

    if let Some(ref path) = files.user_ell_file {
        let ell_def = UserDefinedFractures::from_file(path, true);
        ell_polys.extend(ell_def.create_polys(poly_opts.eps));
        num_user_defined_fractures += ell_def.n_frac;
    }

    if let Some(ref path) = files.user_ell_by_coord_file {
        let ell_def = UserDefinedEllByCoord::from_file(path);
        ell_polys.extend(ell_def.create_polys());
        num_user_defined_fractures += ell_def.n_frac;
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

#[derive(Debug)]
pub struct UserDefinedFractures {
    pub n_frac: usize,
    pub radii: Vec<f64>,
    pub beta: Vec<f64>,
    pub aspect: Vec<f64>,
    pub translation: Vec<[f64; 3]>,
    pub normal: Vec<Vector3<f64>>,
    pub num_points: Vec<usize>,
}

impl UserDefinedFractures {
    pub fn from_file(path: &str, is_ell: bool) -> Self {
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
        }

        UserDefinedFractures {
            n_frac: frac_reader.n_frac,
            radii,
            beta,
            aspect,
            translation,
            normal: normal
                .into_iter()
                .map(|norm| Vector3::from_row_slice(&norm).normalize())
                .collect(),
            num_points,
        }
    }

    pub fn create_polys(&self, eps: f64) -> Vec<Poly> {
        let is_rect = self.num_points.is_empty();

        let mut polys = Vec::new();

        for idx in 0..self.n_frac {
            let number_of_nodes = if is_rect {
                4
            } else {
                self.num_points[idx] as isize
            };
            let family_num = if is_rect { -2 } else { -1 };

            let mut new_poly = Poly {
                // Set number of nodes. Needed for rotations.
                number_of_nodes,
                family_num,
                group_num: 0,
                // Initialize normal to {0,0,1}. need initialized for 3D rotation
                normal: Vector3::new(0., 0., 1.),
                translation: Vector3::from_row_slice(&self.translation[idx]),
                vertices: Vec::with_capacity(number_of_nodes as usize * 3),
                ..Default::default()
            };

            if is_rect {
                // initialize_rect_vertices() sets newpoly.xradius, newpoly.yradius, newpoly.aperture
                initialize_rect_vertices(&mut new_poly, self.radii[idx], self.aspect[idx]);
            } else {
                // Generate theta array used to place vertices
                let theta_ary = generate_theta(self.aspect[idx], self.num_points[idx]);
                // Initialize vertices on x-y plane
                initialize_ell_vertices(
                    &mut new_poly,
                    self.radii[idx],
                    self.aspect[idx],
                    &theta_ary,
                    self.num_points[idx],
                );
            }

            // Apply 2d rotation matrix, twist around origin
            // Assumes polygon on x-y plane
            // Angle must be in rad
            apply_rotation2_d(&mut new_poly, self.beta[idx]);

            // Rotate vertices to uenormal[index] (new normal)
            apply_rotation3_d(&mut new_poly, &self.normal[idx], eps);

            // Save newPoly's new normal vector
            new_poly.normal = self.normal[idx];

            // Translate newPoly to uetranslation
            translate(
                &mut new_poly,
                Vector3::from_row_slice(&self.translation[idx]),
            );

            polys.push(new_poly)
        }

        polys
    }
}

pub struct UserDefinedPolygonByCoord {
    pub n_frac: usize,
    pub num_points: usize,
    pub vertices: Vec<f64>,
}

impl UserDefinedPolygonByCoord {
    pub fn from_file(path: &str) -> Self {
        let mut file = File::open(path).unwrap();

        search_var(&mut file, "nPolygons:");
        let mut n_frac: usize = 0;
        n_frac.read_from_text(&mut file);

        let mut num_points: usize = 0;
        num_points.read_from_text(&mut file);

        let mut vertices = Vec::with_capacity(num_points * 3);

        for _ in 0..num_points {
            let mut tmp = [0., 0., 0.];
            tmp.read_from_text(&mut file);
            vertices.extend(tmp);
        }

        Self {
            n_frac,
            num_points,
            vertices,
        }
    }

    pub fn create_polys(&self) -> Vec<Poly> {
        let mut new_polys = Vec::with_capacity(self.n_frac);

        let n_poly_nodes = self.num_points;

        for _ in 0..self.n_frac {
            let mut new_poly = Poly {
                family_num: -3,
                // Set number of nodes. Needed for rotations.
                number_of_nodes: n_poly_nodes as isize,
                vertices: self.vertices.clone(),
                ..Default::default()
            };

            // Get a normal vector
            // Vector from fist node to node across middle of polygon
            let mut pt_idx_12 = 3 * (n_poly_nodes / 2);

            if n_poly_nodes == 3 {
                pt_idx_12 = 8;
            }

            let p1 = Point3::from_slice(&new_poly.vertices[0..3]);
            let p2 = Point3::from_slice(&new_poly.vertices[3..6]);
            let p_12 = Point3::from_slice(&new_poly.vertices[pt_idx_12..pt_idx_12 + 3]);

            let v1 = p_12 - p1;
            let v2 = p2 - p1;

            new_poly.normal = v2.cross(&v1).normalize();
            // Estimate radius
            // across middle if even number of nodes
            // across middle close to perpendicular to xradius magnitude calculation
            new_poly.xradius = 0.5 * v2.magnitude();

            // Get idx for node 1/4 around polygon
            let pt_idx_14 = 3 * (n_poly_nodes / 4);
            // Get idx for node 3/4 around polygon
            let pt_idx_34 = 3 * (3 * n_poly_nodes / 4);

            let p_14 = Point3::from_slice(&new_poly.vertices[pt_idx_14..pt_idx_14 + 3]);
            let p_34 = Point3::from_slice(&new_poly.vertices[pt_idx_34..pt_idx_34 + 3]);

            new_poly.yradius = 0.5 * distance(&p_14, &p_34);
            new_poly.aspect_ratio = new_poly.yradius / new_poly.xradius;

            // Estimate translation (middle of poly)
            // Use midpoint between 1st and and half way around polygon
            // Note: For polygons defined by coordinates, the coordinates
            // themselves provide the translation. We need to estimate the center
            // of the polygon and init. the translation array
            new_poly.translation = 0.5 * Vector3::new(p1.x + p_12.x, p1.y + p_12.y, p1.z + p_12.z);

            new_polys.push(new_poly)
        }

        new_polys
    }
}

pub struct UserDefinedEllByCoord {
    pub n_frac: usize,
    pub num_points: usize,
    pub vertices: Vec<f64>,
}

impl UserDefinedEllByCoord {
    pub fn from_file(path: &str) -> Self {
        let mut file = File::open(path).unwrap();

        search_var(&mut file, "nEllipses:");
        let mut n_frac = 0;
        n_frac.read_from_text(&mut file);

        search_var(&mut file, "nNodes:");
        let mut num_points = 0;
        num_points.read_from_text(&mut file);

        search_var(&mut file, "Coordinates:");

        let mut bytes = file.bytes().map(|ch| ch.unwrap());
        let size = n_frac * num_points * 3;

        let mut vertices = Vec::with_capacity(size * 3);

        for _ in 0..size {
            let val = read!("{}", bytes);
            vertices.push(val)
        }

        Self {
            n_frac,
            num_points,
            vertices,
        }
    }

    pub fn create_polys(&self) -> Vec<Poly> {
        let mut new_polys = Vec::with_capacity(self.n_frac);

        for idx in 0..self.n_frac {
            let mut new_poly = Poly {
                family_num: -1,
                // Set number of nodes. Needed for rotations.
                number_of_nodes: self.num_points as isize,
                // Initialize normal to {0,0,1}. need initialized for 3D rotation
                normal: Vector3::new(0., 0., 1.),
                ..Default::default()
            };

            new_poly.vertices.reserve(self.num_points * 3);

            let poly_vert_idx = idx * 3 * self.num_points; // Each polygon has nEllNodes * 3 vertices

            // Initialize vertices
            for j in 0..self.num_points {
                let v_idx = j * 3;
                new_poly.vertices[v_idx] = self.vertices[poly_vert_idx + v_idx];
                new_poly.vertices[v_idx + 1] = self.vertices[poly_vert_idx + 1 + v_idx];
                new_poly.vertices[v_idx + 2] = self.vertices[poly_vert_idx + 2 + v_idx];
            }

            // Get a normal vector
            // Vector from fist node to node accross middle of polygon
            let pt_idx_12 = 3 * (self.num_points / 2);
            let p1 = Point3::from_slice(&new_poly.vertices[0..3]);
            let p2 = Point3::from_slice(&new_poly.vertices[3..6]);
            let p_12 = Point3::from_slice(&new_poly.vertices[pt_idx_12..pt_idx_12 + 3]);

            let v1 = p_12 - p1;
            let v2 = p2 - p1;

            new_poly.normal = v2.cross(&v1).normalize();
            // Estimate radius
            // across middle if even number of nodes
            new_poly.xradius = 0.5 * v2.magnitude();

            // Get idx for node 1/4 around polygon
            let pt_idx_14 = 3 * (pt_idx_12 / 2);
            // Get idx for node 3/4 around polygon
            // across middle close to perpendicular to xradius magnitude calculation
            let pt_idx_34 = 3 * (pt_idx_14 + pt_idx_12);

            let p_14 = Point3::from_slice(&new_poly.vertices[pt_idx_14..pt_idx_14 + 3]);
            let p_34 = Point3::from_slice(&new_poly.vertices[pt_idx_34..pt_idx_34 + 3]);

            new_poly.yradius = 0.5 * distance(&p_14, &p_34);
            new_poly.aspect_ratio = new_poly.yradius / new_poly.xradius;

            // Estimate translation (middle of poly)
            // Use midpoint between 1st and and half way around polygon
            // Note: For polygons defined by coordinates, the coordinates
            // themselves provide the translation. We need to estimate the center
            // of the polygon and init. the translation array
            new_poly.translation = 0.5 * Vector3::new(p1.x + p_12.x, p1.y + p_12.y, p1.z + p_12.z);

            new_polys.push(new_poly)
        }

        new_polys
    }
}

pub struct UserDefinedRectByCoord {
    pub n_frac: usize,
    pub vertices: Vec<f64>,
}

impl UserDefinedRectByCoord {
    pub fn from_file(path: &str) -> Self {
        let mut file = File::open(path).unwrap();

        search_var(&mut file, "nRectangles:");
        let mut n_frac = 0;
        n_frac.read_from_text(&mut file);

        search_var(&mut file, "Coordinates:");

        let mut bytes = file.bytes().map(|ch| ch.unwrap());

        let mut vertices = Vec::with_capacity(n_frac * 12);

        for _ in 0..n_frac * 12 {
            let val = read!("{}", bytes);
            vertices.push(val)
        }

        Self { n_frac, vertices }
    }

    pub fn create_polys(&self) -> Vec<Poly> {
        let mut new_polys = Vec::with_capacity(self.n_frac);

        for idx in 0..self.n_frac {
            let mut new_poly = Poly {
                family_num: -2,
                // Set number of nodes. Needed for rotations.
                number_of_nodes: 4,
                ..Default::default()
            };

            new_poly.vertices.reserve(12); // 4 * {x,y,z}
            let poly_vert_idx = idx * 12; // Each polygon has 4 vertices (12 elements, 4*{x,y,z}))

            // Initialize vertices
            for j in 0..4 {
                let v_idx = j * 3;
                new_poly.vertices[v_idx] = self.vertices[poly_vert_idx + v_idx];
                new_poly.vertices[v_idx + 1] = self.vertices[poly_vert_idx + 1 + v_idx];
                new_poly.vertices[v_idx + 2] = self.vertices[poly_vert_idx + 2 + v_idx];
            }

            // Check that rectangle lays one a single plane:
            // let xProd1 = cross Product vector (1st node to 2nd node) with vector(1st node to 3rd node)
            // and xProd2 = cross product vector (1st node to 3th node) with vector (1st node to 4th node)
            // Then, cross product xProd1 and xProd2, if this produces zero vector, all coords are on the same plane
            // v1 is vector from first vertice to third vertice
            // Vector from fist node to 3rd node (vector through middle of sqare)
            let p1 = Point3::from_slice(&new_poly.vertices[0..3]);
            let p2 = Point3::from_slice(&new_poly.vertices[3..6]);
            let p3 = Point3::from_slice(&new_poly.vertices[6..9]);
            let p4 = Point3::from_slice(&new_poly.vertices[9..12]);

            let v1 = p3 - p1;
            let v2 = p2 - p1;
            let v3 = p4 - p1;
            let x_prod1 = v2.cross(&v1).normalize();
            // let x_prod2 = crossProduct(&v3, &v1);
            // let x_prod3 = crossProduct(&x_prod1, &x_prod2);
            //will be zero vector if all vertices are on the same plane
            //TODO: Error check below is too sensitive. Adjust it.
            // Error check for points not on the same plane
            //        if (std::abs(magnitude(xProd3[0],xProd3[1],xProd3[2])) > eps) { //points do not lay on the same plane. reject poly else meshing will fail
            // if (!(std::abs(xProd3[0]) < eps && std::abs(xProd3[1]) < eps && std::abs(xProd3[2]) < eps)) {
            //     delete[] newPoly.vertices;
            //     pstats.rejectedPolyCount++;
            //     std::cout << "\nUser Rectangle (defined by coordinates) " << i+1 << " was rejected. The defined vertices are not co-planar.\n";
            //     std::cout << "Please check user defined coordinates for rectanle " << i+1 << " in input file\n";
            //     delete[] xProd1;
            //     delete[] xProd2;
            //     delete[] xProd3;
            //     continue; //go to next poly
            // }

            // Set normal vector
            new_poly.normal = x_prod1;

            // Set radius (x and y radii might be switched based on order of users coordinates)
            new_poly.xradius = 0.5 * v2.magnitude();
            new_poly.yradius = 0.5 * v3.magnitude();
            new_poly.aspect_ratio = new_poly.yradius / new_poly.xradius;

            // Estimate translation
            // Use midpoint between 1st and 3rd vertices
            // Note: For polygons defined by coordinates, the coordinates
            // themselves provide the translation. We are just filling the
            // translation array for completeness even though the translation
            // array might not be used
            new_poly.translation = 0.5 * Vector3::new(p1.x + p3.x, p1.y + p3.y, p1.z + p3.z);

            new_polys.push(new_poly)
        }

        new_polys
    }
}
