mod domain_truncation;

pub use domain_truncation::DomainBoundaryMode;

use parry3d_f64::bounding_volume::Aabb;
use parry3d_f64::na::{center, Point3, Rotation3, Vector3};
use parry3d_f64::shape::Triangle;

use crate::cg::{is_parallel, rotation_between};

#[derive(Clone)]
/// The Poly structre is used to create and store fracrures/polygons.
pub struct Poly {
    /// Contains the index of the 'shapeFamilies' array in main() from which the fracture was
    /// created from. If the polygon was created from user defined input, it will be marked
    /// with `USER_FRACTURE_FAMILY` (-1).
    /// The stochastic family numbers start at 0 for the first family, and increase by
    /// 1 for each addition family. Families are in the order of the input file.
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

    /// X-radius before fracture-domain truncation.
    /// X-radius is equal to the value generated from randum distributions or given by the user
    /// in the case of constant distributions and user-defined fractures.
    pub xradius: f64,

    /// Y-radius before fracture-domain truncation.
    /// Y-radius is equal to x-radius * aspect ratio (yradius = aspectRatio * xradius).
    pub yradius: f64,

    /// Aspect ratio of polygon before fracture-domain truncation. Must be value greater than zero.
    pub aspect_ratio: f64,

    /// Translation of polygon. This variable is set while building the polygon.
    pub translation: Vector3<f64>,

    /// Polygon normal. This variable is set while building the polygon.
    pub normal: Vector3<f64>,

    /// The bounding box of the polygon. Set with assign_bounding_box().
    pub bounding_box: Aabb,

    /// Vertices of the polygon, in order around its boundary.
    pub vertices: Vec<Point3<f64>>,

    /// The faces array contains flags (true/false) which denote which sides, if any, the polygon is touching.
    /// True (not zero) - Polygon is touching a domain boundary.
    /// False (0) - Polygon is not touching a domain boundary.
    /// Index Key:
    /// [0]: -x face, [1]: +x face
    /// [2]: -y face, [3]: +y face
    /// [4]: -z face, [5]: +z face
    /// Touching boundary faces of the fractures group, not neccesarily the faces of the fracture
    pub faces: [bool; 6],

    /// True if the polygon has been truncated, false otherwise. This variable is used for re-translating polygons.
    /// If the polygon has been truncated, it must be re-built. Otherwise, it can simply be given a new translation. */
    pub truncated: bool,

    /// List of indices to the permanent intersection array ('intPts' in main()) which belong to this polygon.
    pub intersection_index: Vec<usize>,
}

impl Default for Poly {
    fn default() -> Self {
        Self {
            family_num: 0,
            group_num: 0,
            area: 0.,
            xradius: 0.,
            yradius: 0.,
            aspect_ratio: 0.,
            translation: Vector3::zeros(),
            normal: Vector3::zeros(),
            // Intersects nothing until assign_bounding_box() sets it
            bounding_box: Aabb::new_invalid(),
            vertices: Vec::new(),
            faces: [false; 6],
            truncated: false,
            intersection_index: Vec::new(),
        }
    }
}

impl Poly {
    /// Creates a polygon on the x-y plane, with its center at the origin.
    ///
    /// The polygon is a regular polygon with its vertices on a circle of the
    /// given radius, stretched in y by the aspect ratio. 4 points give a
    /// rectangle with its sides parallel to the axes. Many points give an ellipse.
    ///
    /// # Arguments
    ///
    /// * `number_of_points` - Number of vertices
    /// * `radius` - Radius (xradius = radius. yradius = radius * aspectRatio)
    /// * `aspect_ratio` - Aspect ratio
    pub fn new(number_of_points: usize, radius: f64, aspect_ratio: f64) -> Self {
        let angle_step = 2. * std::f64::consts::PI / number_of_points as f64;

        Self {
            normal: Vector3::new(0., 0., 1.),
            xradius: radius,
            yradius: radius * aspect_ratio,
            aspect_ratio,
            // The first vertex is half a step from the x axis, so that
            // the sides of a polygon with 4 points are parallel to the axes
            vertices: (0..number_of_points)
                .map(|k| {
                    let theta = angle_step * (k as f64 + 0.5);
                    Point3::new(
                        radius * theta.cos(),
                        radius * aspect_ratio * theta.sin(),
                        0.,
                    )
                })
                .collect(),
            ..Default::default()
        }
    }

    /// Assign Poly's Area
    ///
    /// Calculate exact area of polygon (after truncation)
    ///
    /// Splits the polygon into triangles around a point inside of it,
    /// and sums the areas of the triangles.
    pub fn assign_area(&mut self) {
        let vertices = &self.vertices;
        let number_of_nodes = vertices.len();

        self.area = if number_of_nodes == 3 {
            Triangle::new(vertices[0], vertices[1], vertices[2]).area()
        } else {
            let inside_pt = center(&vertices[0], &vertices[number_of_nodes / 2]);

            (0..number_of_nodes)
                .map(|i| {
                    Triangle::new(inside_pt, vertices[i], vertices[(i + 1) % number_of_nodes])
                        .area()
                })
                .sum()
        };
    }

    /// Assign Bounding box
    ///
    /// Creates bounding box for polygon/fracture
    /// Sets bounding box in poly struct
    ///
    /// # Arguments
    pub fn assign_bounding_box(&mut self) {
        self.bounding_box = Aabb::from_points(&self.vertices);
    }

    /// 2D rotation matrix
    ///
    /// Rotates poly around its normal vecotor on x-y plane
    /// Assumes poly is on x-y plane
    /// Angle must be in radians
    ///
    /// # Arguments
    ///
    /// * `angle` - Angle to rotate to
    pub fn rotation_2d(&mut self, angle: f64) {
        let sin_calc = angle.sin();
        let cos_calc = angle.cos();

        // Rotates polygon on x-y plane counter-clockwise
        for vertex in &mut self.vertices {
            *vertex = Point3::new(
                (vertex.x * cos_calc) + (vertex.y * sin_calc),
                (vertex.x * -sin_calc) + (vertex.y * cos_calc),
                0.,
            );
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
        if !is_parallel(&self.normal, normal_b, eps) {
            self.rotate(&rotation_between(&self.normal, normal_b));
        }
    }

    /// Applies the rotation to all vertices.
    pub fn rotate(&mut self, rotation: &Rotation3<f64>) {
        for vertex in &mut self.vertices {
            *vertex = rotation * *vertex;
        }
    }

    /// Translate the poly by vector
    ///
    /// # Arguments
    ///
    /// * `translation` - Translation vector
    pub fn translate(&mut self, translation: Vector3<f64>) {
        self.translation = translation;

        for vertex in &mut self.vertices {
            *vertex += translation;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const TOLERANCE: f64 = 1e-12;

    #[test]
    fn test_four_points_give_a_rectangle_with_sides_parallel_to_the_axes() {
        let poly = Poly::new(4, std::f64::consts::SQRT_2, 2.);

        let corners = [[1., 2.], [-1., 2.], [-1., -2.], [1., -2.]];
        for (vertex, corner) in poly.vertices.iter().zip(corners) {
            let offset = vertex - Point3::new(corner[0], corner[1], 0.);
            assert!(offset.norm() < TOLERANCE, "{} is not {:?}", vertex, corner);
        }
    }

    #[test]
    fn test_all_vertices_are_on_the_ellipse() {
        let (radius, aspect_ratio) = (3., 0.5);

        let poly = Poly::new(7, radius, aspect_ratio);

        for vertex in &poly.vertices {
            let ellipse_equation =
                (vertex.x / radius).powi(2) + (vertex.y / (radius * aspect_ratio)).powi(2);
            assert!((ellipse_equation - 1.).abs() < TOLERANCE, "{}", vertex);
        }
    }
}
