mod domain_truncation;

use parry3d_f64::na::Vector3;

use crate::{
    cg::{is_parallel, rotation_matrix},
    distribution::generating_points::generate_theta,
};

#[derive(Clone, Default)]
/// The Poly structre is used to create and store fracrures/polygons.
pub struct Poly {
    /// Number of nodes/vertices in the polygon.
    pub number_of_nodes: isize,

    /// Contains the index of the 'shapeFamilies' array in main() from which the fracture was
    /// created from. If the polygon was created from user defined input, it will be marked
    /// (-2 for user-rect, and -1 for user ell.
    /// The stochastic shape family numbers start at 0 for the first family, and increase by
    /// 1 for each addition family. Families are in the same order which they are defined in the
    /// 'famProb' variable in the input file, starting with ellipse families, then rectangular families.
    pub family_num: isize,

    /// Fracture cluster number which the fracture belongs to. This variable is used to keep
    /// track of fracture connectivity. When a fracture first intersects another fracture,
    /// it inherits its cluster group number (groupNum). If a fracture does not intersect
    /// any other fractures, it is given a new and unique cluster group number. When a
    /// fracture bridges two different clusters, all clusters are merged to be have the
    /// cluster group number of the first intersecting fracture.
    pub group_num: usize,

    /// Polygon area (Not calculated until after DFN generation has completed).
    pub area: f64,

    /// X-radius before fracture-domain truncation. In the case of rectangles, radius is
    /// 1/2 the width of the polygon.
    /// X-radius is equal to the value generated from randum distributions or given by the user
    /// in the case of constant distributions and user-defined fractures.
    pub xradius: f64,

    /// Y-radius before fracture-domain truncation. In the case of rectangles, radius is 1/2
    /// the width of the polygon.
    /// Y-radius is equal to x-radius * aspect ratio (yradius = aspectRatio * xradius).
    pub yradius: f64,

    /// Aspect ratio of polygon before fracture-domain truncation. Must be value greater than zero.
    pub aspect_ratio: f64,

    /// Translation of polygon. This variable is set while building the polygon.
    pub translation: Vector3<f64>,

    /// Polygon normal. This variable is set while building the polygon.
    pub normal: Vector3<f64>,

    /// The bounding box of the polygon. Set with assign_bounding_box().
    /// Index Key:
    /// [0]: x minimum, [1]: x maximum
    /// [2]: y minimum, [3]: y maximum
    /// [4]: z minimum, [5]: z maximum
    pub bounding_box: [f64; 6],

    /// Double array for which hold the polygon's vertices. Vertices are stored in a 1-D array.
    /// e.g For n number of vertices, array will be: {x1, y1, z1, x2, y2, z2, ... , xn, yn, zn}
    pub vertices: Vec<f64>,

    /// The faces array contains flags (true/false) which denote which sides, if any, the polygon is touching.
    /// True (not zero) - Polygon is touching a domain boundary.
    /// False (0) - Polygon is not touching a domain boundary.
    /// Index Key:
    /// [0]: -x face, [1]: +x face
    /// [2]: -y face, [3]: +y face
    /// [4]: -z face, [5]: +z face
    /// Touching boundary faces of the fractures group, not neccesarily the faces of the fracture
    pub faces: [bool; 6],

    /// When writing intersection points and poly vertices to output, we need them to be entirely on the x-y plane.
    /// XYPlane is used for error checking. During intersection rotations, the polygon is rotated to
    /// the x-y plane and it's vertices are changed within the Poly structure. Errors will occur if
    /// the same polygon was to be rotated again because although the vertices are now on the x-y plane, the
    /// normal is still the normal of the polygon in its 3D space. This variable prevents this from happening. */
    pub xyplane: bool,

    /// True if the polygon has been truncated, false otherwise. This variable is used for re-translating polygons.
    /// If the polygon has been truncated, it must be re-built. Otherwise, it can simply be given a new translation. */
    pub truncated: bool,

    /// List of indices to the permanent intersection array ('intPts' in main()) which belong to this polygon.
    pub intersection_index: Vec<usize>,
}

impl Poly {
    /// Create a Rectangular
    ///
    /// # Arguments
    ///
    /// * `radius` - Radius (1/2 x dimension length)
    /// * `aspect_ratio` - Aspect ratio
    pub fn new_rect(radius: f64, aspect_ratio: f64) -> Self {
        let mut new_poly = Self {
            number_of_nodes: 4,
            normal: Vector3::new(0., 0., 1.),
            vertices: Vec::with_capacity(3 * 4),
            ..Default::default()
        };

        // Initializes vertices for rectangular poly using radius (1/2 x length)
        // and aspcet ratio. (xradius = radius, yradius = radius * aspectRatio)
        // Poly will be on x-y plane
        let x = radius;
        let y = radius * aspect_ratio;
        new_poly.xradius = x;
        new_poly.yradius = y;
        new_poly.aspect_ratio = aspect_ratio;
        // Initialize vertices
        new_poly.vertices.push(x);
        new_poly.vertices.push(y);
        new_poly.vertices.push(0.);
        new_poly.vertices.push(-x);
        new_poly.vertices.push(y);
        new_poly.vertices.push(0.);
        new_poly.vertices.push(-x);
        new_poly.vertices.push(-y);
        new_poly.vertices.push(0.);
        new_poly.vertices.push(x);
        new_poly.vertices.push(-y);
        new_poly.vertices.push(0.);

        new_poly
    }

    /// Create a Ellipse
    ///
    /// # Arguments
    ///
    /// * `radius` - Radius (xradius = radius. yradius = radius * aspectRatio)
    /// * `aspect_ratio` - Aspect ratio
    pub fn new_ell(num_nodes: usize, radius: f64, aspect_ratio: f64) -> Self {
        let mut new_poly = Self {
            number_of_nodes: num_nodes as isize,
            normal: Vector3::new(0., 0., 1.),
            vertices: Vec::with_capacity(3 * num_nodes),
            ..Default::default()
        };

        new_poly.xradius = radius;
        new_poly.yradius = radius * aspect_ratio;
        new_poly.aspect_ratio = aspect_ratio;

        // Generate theta array used to place vertices
        let theta_list = generate_theta(aspect_ratio, num_nodes);

        // Initializes ellipse vertices on x-y plane
        for theta in theta_list {
            new_poly.vertices.push(radius * theta.cos());
            new_poly.vertices.push(radius * aspect_ratio * theta.sin());
            new_poly.vertices.push(0.);
        }

        new_poly
    }

    /// Assign Poly's Area
    ///
    /// Calculate exact area of polygon (after truncation)
    ///
    /// Summary of algorithm:
    /// 1: Creates a point on the inside of the polygon (insidePt)
    /// 2: Uses 'insidePt' to create vectors to all outside vertices, breaking the polygon into triangles
    /// 3: Uses .5 * (magnitude of cross product) for each triangle for area calculation
    /// 4: Sums the areas for each trianlge for total area of polygon
    pub fn assign_area(&mut self) {
        let area = if self.number_of_nodes == 3 {
            //area = 1/2 mag of xProd
            let v1 = Vector3::new(
                self.vertices[3] - self.vertices[0],
                self.vertices[4] - self.vertices[1],
                self.vertices[5] - self.vertices[2],
            );
            let v2 = Vector3::new(
                self.vertices[6] - self.vertices[0],
                self.vertices[7] - self.vertices[1],
                self.vertices[8] - self.vertices[2],
            );
            let x_prod = v1.cross(&v2);

            0.5 * x_prod.magnitude()
        } else {
            // More than 3 vertices
            let mut poly_area = 0.; // For summing area over trianlges of polygon
                                    // Get coordinate within polygon
            let idx_across = (self.number_of_nodes / 2 * 3) as usize;
            let inside_pt = [
                self.vertices[0] + (0.5 * (self.vertices[idx_across] - self.vertices[0])),
                self.vertices[1] + (0.5 * (self.vertices[idx_across + 1] - self.vertices[1])),
                self.vertices[2] + (0.5 * (self.vertices[idx_across + 2] - self.vertices[2])),
            ];

            for i in 0..self.number_of_nodes - 1 {
                let idx = (i * 3) as usize;
                let v1 = Vector3::new(
                    self.vertices[idx] - inside_pt[0],
                    self.vertices[idx + 1] - inside_pt[1],
                    self.vertices[idx + 2] - inside_pt[2],
                );
                let v2 = Vector3::new(
                    self.vertices[idx + 3] - inside_pt[0],
                    self.vertices[idx + 4] - inside_pt[1],
                    self.vertices[idx + 5] - inside_pt[2],
                );
                let x_prod = v1.cross(&v2);
                let area = 0.5 * x_prod.magnitude();
                poly_area += area; // Accumulate area
            }

            // Last portion of polygon, insidePt to first vertice and insidePt to last vertice
            let last = (3 * (self.number_of_nodes - 1)) as usize;
            let v1 = Vector3::new(
                self.vertices[0] - inside_pt[0],
                self.vertices[1] - inside_pt[1],
                self.vertices[2] - inside_pt[2],
            );
            let v2 = Vector3::new(
                self.vertices[last] - inside_pt[0],
                self.vertices[last + 1] - inside_pt[1],
                self.vertices[last + 2] - inside_pt[2],
            );
            let x_prod = v1.cross(&v2);
            let area = 0.5 * x_prod.magnitude();
            poly_area += area; // Accumulate area

            poly_area
        };

        self.area = area;
    }

    /// Assign Bounding box
    ///
    /// Creates bounding box for polygon/fracture
    /// Sets bounding box in poly struct
    ///
    /// # Arguments
    pub fn assign_bounding_box(&mut self) {
        // Initialize mins and maxs
        let mut min_x = self.vertices[0]; // x1
        let mut max_x = min_x;
        let mut min_y = self.vertices[1]; // y1
        let mut max_y = min_y;
        let mut min_z = self.vertices[2]; // z1
        let mut max_z = min_z;

        for i in 1..self.number_of_nodes {
            let idx = (i * 3) as usize;
            max_x = f64::max(max_x, self.vertices[idx]);
            min_x = f64::min(min_x, self.vertices[idx]);
            max_y = f64::max(max_y, self.vertices[idx + 1]);
            min_y = f64::min(min_y, self.vertices[idx + 1]);
            max_z = f64::max(max_z, self.vertices[idx + 2]);
            min_z = f64::min(min_z, self.vertices[idx + 2]);
        }

        self.bounding_box = [min_x, max_x, min_y, max_y, min_z, max_z];
    }

    /// 2D rotation matrix
    ///
    /// Rotates poly around its normal vecotor on x-y plane
    /// Assumes poly is on x-y plane
    /// Assumes poly.number_of_nodes is set
    /// Angle must be in radians
    ///
    /// # Arguments
    ///
    /// * `angle` - Angle to rotate to
    pub fn rotation_2d(&mut self, angle: f64) {
        let sin_calc = angle.sin();
        let cos_calc = angle.cos();

        // Rotates polygon on x-y plane counter-clockwise
        for i in 0..self.number_of_nodes {
            let idx = (i * 3) as usize;
            let x = self.vertices[idx];
            let y = self.vertices[idx + 1];
            self.vertices[idx] = (x * cos_calc) + (y * sin_calc); // x
            self.vertices[idx + 1] = (x * -sin_calc) + (y * cos_calc); // y
            self.vertices[idx + 2] = 0.; // z
        }
    }

    /// Applies a Rotation Matrix to poly vertices
    ///
    /// RotMatrix = I + V + V^2((1-cos)/sin^2)) rotate to new normal
    ///
    /// Rotates 'newPoly' from newPoly's current normal
    /// so that newPoly's new normal will be 'normalB'
    /// Assumes poly.numberOfPoints and newPoly.normal are initialized and normalized
    ///
    /// # Arguments
    ///
    /// * `new_poly` - Poly to be rotated
    /// * `normal_b` - Normal vector to rotate to
    /// * `eps` - Epsilon value for floating point comparisons
    pub fn rotation_3d(&mut self, normal_b: &Vector3<f64>, eps: f64) {
        // Normals should already be normalized by this point!!!
        let normal_a = self.normal;

        if !is_parallel(&normal_a, normal_b, eps) {
            // NOTE: rotationMatrix() requires normals to be normalized
            let r = rotation_matrix(&normal_a, normal_b, eps);

            // Apply rotation to all vertices
            for i in 0..self.number_of_nodes {
                let idx = (i * 3) as usize;

                let x = self.vertices[idx] * r[0]
                    + self.vertices[idx + 1] * r[1]
                    + self.vertices[idx + 2] * r[2];
                let y = self.vertices[idx] * r[3]
                    + self.vertices[idx + 1] * r[4]
                    + self.vertices[idx + 2] * r[5];
                let z = self.vertices[idx] * r[6]
                    + self.vertices[idx + 1] * r[7]
                    + self.vertices[idx + 2] * r[8];

                self.vertices[idx] = x;
                self.vertices[idx + 1] = y;
                self.vertices[idx + 2] = z;
            }
        }
    }

    /// Translate the poly by vector
    ///
    /// Assumes newPoly.number_of_nodes is initialized
    ///
    /// # Arguments
    ///
    /// * `translation` - Translation vector
    pub fn translate(&mut self, translation: Vector3<f64>) {
        self.translation = translation;

        for i in 0..self.number_of_nodes {
            let idx = (i * 3) as usize;
            self.vertices[idx] += translation.x;
            self.vertices[idx + 1] += translation.y;
            self.vertices[idx + 2] += translation.z;
        }
    }
}
