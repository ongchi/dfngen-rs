use parry3d_f64::na::{Point3, Vector3};
use parry3d_f64::query::details::clip_halfspace_polygon;

use crate::fracture::poly::Poly;

/// What to do with the part of a fracture which is outside the domain
#[derive(Default, Debug, Clone, Copy, PartialEq, Eq)]
pub enum DomainBoundaryMode {
    /// Cut the fracture at the domain boundary.
    #[default]
    Truncate,
    /// Keep the fracture whole.
    KeepWhole,
}

/// Cuts a polygon with the 6 sides of the domain, then removes a vertex of each
/// edge which is shorter than 2h
///
/// # Arguments
///
/// * `vertices` - Vertices of the polygon
/// * `h` - Minimum feature size
/// * `eps` - Epsilon value for floating point comparisons
/// * `half_domain` - Half of the domain size in x, y, z
///
/// # Returns
///
/// Vertices of the part in the domain. Fewer than 3 vertices: the polygon has no
/// part in the domain.
fn clip_to_domain(
    vertices: &[Point3<f64>],
    h: f64,
    eps: f64,
    half_domain: &Vector3<f64>,
) -> Vec<Point3<f64>> {
    // Outward normal of each side of the domain, and the center point of the side
    let domain_sides = [
        (Vector3::new(0., 0., 1.), Point3::new(0., 0., half_domain.z)),
        (
            Vector3::new(0., 0., -1.),
            Point3::new(0., 0., -half_domain.z),
        ),
        (Vector3::new(0., 1., 0.), Point3::new(0., half_domain.y, 0.)),
        (
            Vector3::new(0., -1., 0.),
            Point3::new(0., -half_domain.y, 0.),
        ),
        (Vector3::new(1., 0., 0.), Point3::new(half_domain.x, 0., 0.)),
        (
            Vector3::new(-1., 0., 0.),
            Point3::new(-half_domain.x, 0., 0.),
        ),
    ];

    // Cut the polygon with each side of the domain
    let mut vertices = vertices.to_vec();
    let mut clipped = Vec::with_capacity(vertices.len() + 6);

    for (side_normal, side_center) in domain_sides {
        let Some(last) = vertices.last() else {
            return vertices;
        };
        let is_last_inside = (last - side_center).dot(&side_normal) <= 0.;

        clip_halfspace_polygon(&side_center, &side_normal, &vertices, &mut clipped);

        if is_last_inside {
            // The result starts with the last vertex. Put it back at the end,
            // so that the first vertex of the polygon stays the first one.
            clipped.rotate_left(1);
        }

        std::mem::swap(&mut vertices, &mut clipped);
    }

    // Remove a vertex of each edge which is shorter than 2h
    let mut i = 0;

    while i < vertices.len() {
        let vertex = vertices[i];
        let next = (i + 1) % vertices.len();

        if (vertex - vertices[next]).magnitude() < (2. * h) {
            // If point is NOT on a boundary, delete current indexed point, ELSE delete next point
            let is_on_boundary =
                (0..3).any(|axis| (vertex[axis].abs() - half_domain[axis]).abs() <= eps);

            vertices.remove(if is_on_boundary { next } else { i });
        } else {
            i += 1;
        }
    }

    vertices
}

/// Finds the domain faces (+x, -x, +y, -y, +z, -z) on which a polygon has a vertex
fn find_domain_faces(vertices: &[Point3<f64>], eps: f64, half_domain: &Vector3<f64>) -> [bool; 6] {
    let mut faces = [false; 6];

    for vertex in vertices {
        for axis in 0..3 {
            if vertex[axis] >= half_domain[axis] - eps {
                faces[2 * axis] = true;
            } else if vertex[axis] <= -half_domain[axis] + eps {
                faces[2 * axis + 1] = true;
            }
        }
    }

    faces
}

impl Poly {
    /// Domain Truncation
    ///
    /// Finds the part of the polygon which is in the domain, and the domain faces
    /// which the polygon touches. With `DomainBoundaryMode::Truncate`, the polygon
    /// becomes that part. With `DomainBoundaryMode::KeepWhole`, the vertices of
    /// the polygon do not change.
    ///
    /// # Arguments
    ///
    /// * `h` - Minimum feature size
    /// * `eps` - Epsilon value for floating point comparisons
    /// * `domain_size` - Domain size in x, y, z
    /// * `mode` - What to do with the part outside the domain
    ///
    /// # Returns
    ///
    /// false - If the polygon has a part with 3 vertices or more in the domain
    /// true - If rejected, because the part in the domain has less than 3 vertices
    pub fn domain_truncation(
        &mut self,
        h: f64,
        eps: f64,
        domain_size: &Vector3<f64>,
        mode: DomainBoundaryMode,
    ) -> bool {
        let half_domain = domain_size * 0.5;
        let is_inside_domain = |vertex: &Point3<f64>| {
            (0..3)
                .all(|axis| vertex[axis] <= half_domain[axis] && vertex[axis] >= -half_domain[axis])
        };

        if mode == DomainBoundaryMode::KeepWhole {
            // A whole polygon is not rebuilt for a new translation. Thus the
            // faces of the last translation must be removed here.
            self.faces = [false; 6];
        }

        if self.vertices.iter().all(is_inside_domain) {
            // Truncation is not necessary
            return false;
        }

        let part_in_domain = clip_to_domain(&self.vertices, h, eps, &half_domain);
        let is_rejected = part_in_domain.len() < 3;

        if !is_rejected {
            self.faces = find_domain_faces(&part_in_domain, eps, &half_domain);
        }

        if mode == DomainBoundaryMode::Truncate {
            self.truncated = true;
            self.vertices = part_in_domain;
        }

        is_rejected
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const H: f64 = 0.1;
    const EPS: f64 = 1e-9;

    /// A square with sides of 2 in the plane z = 0 and with its center on the x axis
    fn square_in_domain_of_size_10(center_x: f64, mode: DomainBoundaryMode) -> (Poly, bool) {
        let mut square = Poly::new(4, std::f64::consts::SQRT_2, 1.);
        square.translate(Vector3::new(center_x, 0., 0.));

        let is_rejected = square.domain_truncation(H, EPS, &Vector3::new(10., 10., 10.), mode);

        (square, is_rejected)
    }

    fn max_x(poly: &Poly) -> f64 {
        poly.vertices
            .iter()
            .map(|vertex| vertex.x)
            .fold(f64::MIN, f64::max)
    }

    #[test]
    fn test_truncated_fracture_ends_at_the_face_which_it_crosses() {
        let (square, _) = square_in_domain_of_size_10(4.5, DomainBoundaryMode::Truncate);

        assert!((max_x(&square) - 5.).abs() < EPS);
    }

    #[test]
    fn test_truncated_fracture_is_marked_as_truncated() {
        let (square, _) = square_in_domain_of_size_10(4.5, DomainBoundaryMode::Truncate);

        assert!(square.truncated);
    }

    #[test]
    fn test_whole_fracture_keeps_its_vertices_when_it_crosses_a_face() {
        let (square, _) = square_in_domain_of_size_10(4.5, DomainBoundaryMode::KeepWhole);

        assert!((max_x(&square) - 5.5).abs() < EPS);
    }

    #[test]
    fn test_whole_fracture_is_connected_to_the_face_which_it_crosses() {
        let (square, _) = square_in_domain_of_size_10(4.5, DomainBoundaryMode::KeepWhole);

        // Only the +x face
        assert_eq!(square.faces, [true, false, false, false, false, false]);
    }

    #[test]
    fn test_whole_fracture_is_not_marked_as_truncated() {
        let (square, _) = square_in_domain_of_size_10(4.5, DomainBoundaryMode::KeepWhole);

        assert!(!square.truncated);
    }

    #[test]
    fn test_whole_fracture_with_a_part_in_the_domain_is_accepted() {
        let (_, is_rejected) = square_in_domain_of_size_10(4.5, DomainBoundaryMode::KeepWhole);

        assert!(!is_rejected);
    }

    #[test]
    fn test_whole_fracture_outside_the_domain_is_rejected() {
        let (_, is_rejected) = square_in_domain_of_size_10(7., DomainBoundaryMode::KeepWhole);

        assert!(is_rejected);
    }

    #[test]
    fn test_whole_fracture_loses_the_faces_of_its_last_position() {
        let (mut square, _) = square_in_domain_of_size_10(4.5, DomainBoundaryMode::KeepWhole);
        square.translate(Vector3::new(-4.5, 0., 0.));

        square.domain_truncation(
            H,
            EPS,
            &Vector3::new(10., 10., 10.),
            DomainBoundaryMode::KeepWhole,
        );

        assert_eq!(square.faces, [false; 6]);
    }
}
