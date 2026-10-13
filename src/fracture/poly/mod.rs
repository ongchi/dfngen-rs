mod box_domain;

use parry3d_f64::bounding_volume::Aabb;
use parry3d_f64::na::{center, Point3, Rotation3, Vector3};
use parry3d_f64::shape::Triangle;

use crate::cg::{is_parallel, rotation_between};

/// Where a fracture comes from
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FractureOrigin {
    /// Stochastic family, by its index in the families of the input file (starts at 0)
    Family(usize),
    /// User-defined fracture file
    UserDefined,
}

impl FractureOrigin {
    /// Index of the stochastic family, `None` for a user-defined fracture
    pub fn family_index(self) -> Option<usize> {
        match self {
            Self::Family(index) => Some(index),
            Self::UserDefined => None,
        }
    }

    /// Family number of the output files: 0 for a user-defined fracture,
    /// and the family index plus 1 for a stochastic family
    pub fn output_number(self) -> usize {
        self.family_index().map_or(0, |index| index + 1)
    }
}

#[derive(Clone)]
/// The Poly structre is used to create and store fracrures/polygons.
pub struct Poly {
    /// Stochastic family or user-defined fracture file which the fracture comes from.
    pub origin: FractureOrigin,

    /// Value of the `fracture_id` array of the user-defined fracture file.
    /// Polygons with the same value are pieces of one fracture. `None` for a
    /// stochastic fracture and for a file without the array.
    pub user_fracture_id: Option<u64>,

    /// Cluster group of the fracture (see `ClusterGroups`). It keeps track of the fracture
    /// connectivity. A fracture gets the group of the first fracture which it intersects.
    /// A fracture which intersects no fracture gets a new group. When a fracture bridges
    /// two groups, the groups are merged into the group of the first intersecting fracture.
    /// `None` until the fracture is checked against the accepted fractures.
    pub group: Option<usize>,

    /// Polygon area (Not calculated until after DFN generation has completed).
    pub area: f64,

    /// X-radius.
    /// X-radius is equal to the value generated from randum distributions or given by the user
    /// in the case of constant distributions and user-defined fractures.
    pub x_radius: f64,

    /// Y-radius.
    /// Y-radius is equal to x-radius * aspect ratio (y_radius = aspectRatio * x_radius).
    pub y_radius: f64,

    /// Aspect ratio of polygon. Must be value greater than zero.
    pub aspect_ratio: f64,

    /// Translation of polygon. This variable is set while building the polygon.
    pub translation: Vector3<f64>,

    /// Polygon normal. This variable is set while building the polygon.
    pub normal: Vector3<f64>,

    /// The bounding box of the polygon. Set with assign_bounding_box().
    pub bounding_box: Aabb,

    /// Vertices of the polygon, in order around its boundary.
    pub vertices: Vec<Point3<f64>>,

    /// List of indices to the permanent intersection array ('intPts' in main()) which belong to this polygon.
    pub intersection_indices: Vec<usize>,
}

impl Default for Poly {
    fn default() -> Self {
        Self {
            origin: FractureOrigin::Family(0),
            user_fracture_id: None,
            group: None,
            area: 0.,
            x_radius: 0.,
            y_radius: 0.,
            aspect_ratio: 0.,
            translation: Vector3::zeros(),
            normal: Vector3::zeros(),
            // Intersects nothing until assign_bounding_box() sets it
            bounding_box: Aabb::new_invalid(),
            vertices: Vec::new(),
            intersection_indices: Vec::new(),
        }
    }
}

impl Poly {
    /// True if the two polygons are pieces of the same user-defined fracture.
    pub fn is_piece_of_same_fracture(&self, other: &Poly) -> bool {
        self.user_fracture_id.is_some() && self.user_fracture_id == other.user_fracture_id
    }

    /// Creates a polygon on the x-y plane, with its center at the origin.
    ///
    /// The polygon is a regular polygon with its vertices on a circle of the
    /// given radius, stretched in y by the aspect ratio. 4 points give a
    /// rectangle with its sides parallel to the axes. Many points give an ellipse.
    ///
    /// # Arguments
    ///
    /// * `number_of_points` - Number of vertices
    /// * `radius` - Radius (x_radius = radius. y_radius = radius * aspectRatio)
    /// * `aspect_ratio` - Aspect ratio
    pub fn new(number_of_points: usize, radius: f64, aspect_ratio: f64) -> Self {
        let angle_step = 2. * std::f64::consts::PI / number_of_points as f64;

        Self {
            normal: Vector3::new(0., 0., 1.),
            x_radius: radius,
            y_radius: radius * aspect_ratio,
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
    /// Calculate exact area of polygon
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
