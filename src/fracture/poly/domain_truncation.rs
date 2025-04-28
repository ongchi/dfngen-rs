use parry3d_f64::na::Vector3;

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
        let mut points = Vec::with_capacity(18);
        let mut n_nodes = 0;
        let domain_x = domain_size[0] * 0.5;
        let domain_y = domain_size[1] * 0.5;
        let domain_z = domain_size[2] * 0.5;
        let stop = self.number_of_nodes;

        // Check if truncation is necessay:
        for k in 0..stop {
            let idx = (k * 3) as usize;

            if self.vertices[idx] > domain_x || self.vertices[idx] < -domain_x {
                break;
            }

            if self.vertices[idx + 1] > domain_y || self.vertices[idx + 1] < -domain_y {
                break;
            }

            if self.vertices[idx + 2] > domain_z || self.vertices[idx + 2] < -domain_z {
                break;
            }

            // If finished checking last node and code still hasn't broke
            // from the for loop, then all nodes are inside domain
            if k + 1 == stop {
                return false;
            }
        }

        // If code does not return above, truncation is needed/
        self.truncated = true; // Mark poly as truncated
        let mut n_vertices; // Same as newPoly.numberOfNodes;
        let mut ntmp = Vector3::new(0., 0., 0.); // Normal of domain side
        let mut pttmp = Vector3::new(0., 0., 0.); // Center point on domain side

        // Check against the all the walls of the domain
        for j in 0..6 {
            match j {
                0 => {
                    ntmp[0] = 0.;
                    ntmp[1] = 0.;
                    ntmp[2] = 1.;
                    pttmp[0] = 0.;
                    pttmp[1] = 0.;
                    pttmp[2] = domain_z;
                }

                1 => {
                    ntmp[0] = 0.;
                    ntmp[1] = 0.;
                    ntmp[2] = -1.;
                    pttmp[0] = 0.;
                    pttmp[1] = 0.;
                    pttmp[2] = -domain_z;
                }

                2 => {
                    ntmp[0] = 0.;
                    ntmp[1] = 1.;
                    ntmp[2] = 0.;
                    pttmp[0] = 0.;
                    pttmp[1] = domain_y;
                    pttmp[2] = 0.;
                }

                3 => {
                    ntmp[0] = 0.;
                    ntmp[1] = -1.;
                    ntmp[2] = 0.;
                    pttmp[0] = 0.;
                    pttmp[1] = -domain_y;
                    pttmp[2] = 0.;
                }

                4 => {
                    ntmp[0] = 1.;
                    ntmp[1] = 0.;
                    ntmp[2] = 0.;
                    pttmp[0] = domain_x;
                    pttmp[1] = 0.;
                    pttmp[2] = 0.;
                }

                5 => {
                    ntmp[0] = -1.;
                    ntmp[1] = 0.;
                    ntmp[2] = 0.;
                    pttmp[0] = -domain_x;
                    pttmp[1] = 0.;
                    pttmp[2] = 0.;
                }
                _ => unreachable!(),
            }

            n_vertices = self.number_of_nodes;

            if n_vertices <= 0 {
                return true; // Reject
            }

            n_nodes = 0; // Counter for final number of nodes
            points.clear(); // Clear any points from last iteration
            let last = ((n_vertices - 1) * 3) as usize; // Index to last node starting position in array
            let mut temp = Vector3::new(
                self.vertices[last] - pttmp[0],
                self.vertices[last + 1] - pttmp[1],
                self.vertices[last + 2] - pttmp[2],
            );
            // Previous distance - the dot product of the domain side normal
            // and the distance between vertex number nVertices and the temp point on the
            // domain side.
            let mut prevdist = temp.dot(&ntmp);

            for i in 0..n_vertices {
                let index = (i * 3) as usize;
                temp.x = self.vertices[index] - pttmp.x;
                temp.y = self.vertices[index + 1] - pttmp.y;
                temp.z = self.vertices[index + 2] - pttmp.z;
                // Current distance, the dot product of domain side normal and
                // the distance between the ii'th vertex and the temporary point
                let currdist = temp.dot(&ntmp);

                if currdist <= 0. {
                    // if vertex is towards the domain relative to the domain side
                    // Save vertex
                    points.push(self.vertices[index]); // x
                    points.push(self.vertices[index + 1]); // y
                    points.push(self.vertices[index + 2]); // z
                    n_nodes += 1;
                }

                if currdist * prevdist < 0. {
                    //if crosses boundary
                    n_nodes += 1;

                    if i == 0 {
                        // Store point on boundary
                        // 'last' is index to last vertice
                        temp[0] = self.vertices[last]
                            + (self.vertices[0] - self.vertices[last]) * prevdist.abs()
                                / (currdist.abs() + prevdist.abs());
                        temp[1] = self.vertices[last + 1]
                            + (self.vertices[1] - self.vertices[last + 1]) * prevdist.abs()
                                / (currdist.abs() + prevdist.abs());
                        temp[2] = self.vertices[last + 2]
                            + (self.vertices[2] - self.vertices[last + 2]) * prevdist.abs()
                                / (currdist.abs() + prevdist.abs());
                        points.push(temp[0]);
                        points.push(temp[1]);
                        points.push(temp[2]);
                    } else {
                        // If from outside to inside
                        temp[0] = self.vertices[index - 3]
                            + (self.vertices[index] - self.vertices[index - 3]) * prevdist.abs()
                                / (currdist.abs() + prevdist.abs());
                        temp[1] = self.vertices[index - 2]
                            + (self.vertices[index + 1] - self.vertices[index - 2])
                                * prevdist.abs()
                                / (currdist.abs() + prevdist.abs());
                        temp[2] = self.vertices[index - 1]
                            + (self.vertices[index + 2] - self.vertices[index - 1])
                                * prevdist.abs()
                                / (currdist.abs() + prevdist.abs());
                        points.reserve(points.len() + 3); // Resize vector if needed
                        points.push(temp[0]);
                        points.push(temp[1]);
                        points.push(temp[2]);
                    }

                    if currdist < 0. {
                        // If from outside to inside - swap order
                        let tmp_index = (n_nodes - 1) * 3; // Index to last saved point
                        points[tmp_index] = points[tmp_index - 3];
                        points[tmp_index + 1] = points[tmp_index - 2];
                        points[tmp_index + 2] = points[tmp_index - 1];
                        points[tmp_index - 3] = temp[0];
                        points[tmp_index - 2] = temp[1];
                        points[tmp_index - 1] = temp[2];
                    }
                }

                prevdist = currdist;
            } // End vertice loops

            self.number_of_nodes = n_nodes as isize;

            if n_nodes > 0 {
                self.vertices = Vec::with_capacity(3 * n_nodes);

                // Copy new nodes back to newPoly
                for k in 0..n_nodes {
                    let idxx = k * 3;
                    self.vertices.push(points[idxx]);
                    self.vertices.push(points[idxx + 1]);
                    self.vertices.push(points[idxx + 2]);
                }
            }
        } // End main loop

        let mut i = 0;

        while i < n_nodes {
            let idx = i * 3;

            let next = if i == n_nodes - 1 {
                // If last node
                0
            } else {
                (i + 1) * 3
            };

            let temp = Vector3::new(
                self.vertices[idx] - self.vertices[next],
                self.vertices[idx + 1] - self.vertices[next + 1],
                self.vertices[idx + 2] - self.vertices[next + 2],
            );

            if temp.magnitude() < (2. * h) {
                // If distance between current and next vertex < h
                // If point is NOT on a boundary, delete current indexed point, ELSE delete next point
                if (self.vertices[idx].abs() - domain_x).abs() > eps
                    && (self.vertices[idx + 1].abs() - domain_y).abs() > eps
                    && (self.vertices[idx + 2].abs() - domain_z).abs() > eps
                {
                    // Loop deletes a vertice by shifting elements to the right of the element, to the left
                    for j in i..n_nodes - 1 {
                        let idx = j * 3;
                        self.vertices[idx] = self.vertices[idx + 3]; // Shift x coordinate left
                        self.vertices[idx + 1] = self.vertices[idx + 4]; // Shift y coordinate left
                        self.vertices[idx + 2] = self.vertices[idx + 5];
                        // Shift z coordinate left
                    }
                } else {
                    // Loop deletes a vertice by shifting elements to the right of the element to the left
                    let end = (n_nodes - 1) * 3;

                    for j in (next..end).step_by(3) {
                        self.vertices[j] = self.vertices[j + 3]; // Shift x coordinate left
                        self.vertices[j + 1] = self.vertices[j + 4]; // Shift y coordinate left
                        self.vertices[j + 2] = self.vertices[j + 5]; // Shift z coordinate left
                    }
                }

                n_nodes -= 1; // Update node counter
            } else {
                i += 1;
            }
        } // End while loop

        if n_nodes < 3 {
            return true; // Reject
        }

        self.number_of_nodes = n_nodes as isize;
        let temp = self.number_of_nodes;

        for k in 0..temp {
            // Check boundary faces
            // Update which boundaries newPoly touches
            let idx = (k * 3) as usize;

            if self.vertices[idx] >= domain_x - eps {
                self.faces[0] = true;
            } else if self.vertices[idx] <= -domain_x + eps {
                self.faces[1] = true;
            }

            if self.vertices[idx + 1] >= domain_y - eps {
                self.faces[2] = true;
            } else if self.vertices[idx + 1] <= -domain_y + eps {
                self.faces[3] = true;
            }

            if self.vertices[idx + 2] >= domain_z - eps {
                self.faces[4] = true;
            } else if self.vertices[idx + 2] <= -domain_z + eps {
                self.faces[5] = true;
            }
        }

        false
    }
}
