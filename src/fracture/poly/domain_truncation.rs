use parry3d_f64::na::{Point3, Vector3};
use parry3d_f64::query::details::clip_halfspace_polygon;

use crate::fracture::poly::Poly;

/// Domain Truncation
///
/// NOTE:
///     domainTruncation() may benefit from rewriting in a more
///     efficient way which does not reallocate memory as often
///     Truncates polygons along the defined domain ('domainSize' in input file)
///
/// # Arguments
///
/// * `h` - Minimum feature size
/// * `eps` - Epsilon value for floating point comparisons
/// * `new_poly` - Polygon being truncated (if truncation is necessary)
/// * `domain_size` - Point to domain size array, 3 doubles: {x, y, z}
///
/// # Returns
///
/// 0 - If Poly is inside domain and was truncated to more than 2 vertices,
///     of poly truncation was not needed
/// 1 - If rejected due to being outside the domain or was truncated to
///     less than 3 vertices
impl Poly {
    pub fn domain_truncation(&mut self, h: f64, eps: f64, domain_size: &Vector3<f64>) -> bool {
        let half_domain = domain_size * 0.5;
        let is_inside_domain = |vertex: &Point3<f64>| {
            (0..3)
                .all(|axis| vertex[axis] <= half_domain[axis] && vertex[axis] >= -half_domain[axis])
        };

        if self.vertices.iter().all(is_inside_domain) {
            // Truncation is not necessary
            return false;
        }

        self.truncated = true;

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
        let mut clipped = Vec::with_capacity(self.vertices.len() + 6);

        for (side_normal, side_center) in domain_sides {
            let Some(last) = self.vertices.last() else {
                return true; // Reject
            };
            let is_last_inside = (last - side_center).dot(&side_normal) <= 0.;

            clip_halfspace_polygon(&side_center, &side_normal, &self.vertices, &mut clipped);

            if is_last_inside {
                // The result starts with the last vertex. Put it back at the end,
                // so that the first vertex of the polygon stays the first one.
                clipped.rotate_left(1);
            }

            std::mem::swap(&mut self.vertices, &mut clipped);
        }

        // Remove a vertex of each edge which is shorter than 2h
        let mut i = 0;

        while i < self.vertices.len() {
            let vertex = self.vertices[i];
            let next = (i + 1) % self.vertices.len();

            if (vertex - self.vertices[next]).magnitude() < (2. * h) {
                // If point is NOT on a boundary, delete current indexed point, ELSE delete next point
                let is_on_boundary =
                    (0..3).any(|axis| (vertex[axis].abs() - half_domain[axis]).abs() <= eps);

                self.vertices.remove(if is_on_boundary { next } else { i });
            } else {
                i += 1;
            }
        }

        if self.vertices.len() < 3 {
            return true; // Reject
        }

        // Update which boundaries the polygon touches
        for vertex in &self.vertices {
            for axis in 0..3 {
                if vertex[axis] >= half_domain[axis] - eps {
                    self.faces[2 * axis] = true;
                } else if vertex[axis] <= -half_domain[axis] + eps {
                    self.faces[2 * axis + 1] = true;
                }
            }
        }

        false
    }
}
