//! The points which two or more fractures share: the points of the
//! intersections, and the points on the common edges and corners of the
//! pieces of a user-defined fracture
//!
//! Each shared line is divided once, in space. The fractures on both sides
//! use the same points, thus the mesh is connected across the line.

use std::collections::HashMap;

use parry3d_f64::na::{distance, Point3};

use super::plane::{segment_crossing, FracturePlane};
use super::MERGE_DISTANCE_PART;
use crate::region_boundaries::{distance_to_segment, position_on_segment};

/// Points in space which two or more fractures share. A point with a
/// "group" can merge with another point of the same group; the group of a
/// fracture is its own position, or the first piece of its user-defined
/// fracture.
pub struct SharedPoints {
    pub points: Vec<Point3<f64>>,
    /// Groups of each point (empty for a point which does not merge)
    groups: Vec<Vec<usize>>,
    tolerance: f64,
    /// Points which can merge, by the grid cell of size `tolerance` which has them
    cells: HashMap<[i64; 3], Vec<usize>>,
}

impl SharedPoints {
    fn new(tolerance: f64) -> Self {
        Self {
            points: Vec::new(),
            groups: Vec::new(),
            tolerance,
            cells: HashMap::new(),
        }
    }

    fn cell_of(&self, point: &Point3<f64>) -> [i64; 3] {
        [point.x, point.y, point.z].map(|coordinate| (coordinate / self.tolerance).floor() as i64)
    }

    /// The point nearer than the tolerance which has one of the groups, or a
    /// new point. The point gets the groups.
    fn find_or_insert(&mut self, point: Point3<f64>, groups: &[usize]) -> usize {
        let [x, y, z] = self.cell_of(&point);
        let nearest = (x - 1..=x + 1)
            .flat_map(|cell_x| (y - 1..=y + 1).map(move |cell_y| (cell_x, cell_y)))
            .flat_map(|(cell_x, cell_y)| {
                (z - 1..=z + 1).map(move |cell_z| [cell_x, cell_y, cell_z])
            })
            .filter_map(|cell| self.cells.get(&cell))
            .flatten()
            .copied()
            .filter(|&id| distance(&self.points[id], &point) < self.tolerance)
            .filter(|&id| self.groups[id].iter().any(|group| groups.contains(group)))
            .min_by(|&a, &b| {
                distance(&self.points[a], &point).total_cmp(&distance(&self.points[b], &point))
            });

        if let Some(id) = nearest {
            for &group in groups {
                if !self.groups[id].contains(&group) {
                    self.groups[id].push(group);
                }
            }
            return id;
        }

        let id = self.add(point);
        self.groups[id] = groups.to_vec();
        self.cells.entry(self.cell_of(&point)).or_default().push(id);
        id
    }

    /// A new point which does not merge with other points
    fn add(&mut self, point: Point3<f64>) -> usize {
        self.points.push(point);
        self.groups.push(Vec::new());
        self.points.len() - 1
    }
}
/// A line which two fractures share, divided into points
pub struct SharedLine {
    /// Positions of the two fractures in the final fractures
    pub fractures: [usize; 2],
    /// Ids of the points in `SharedPoints`, in order along the line
    pub point_ids: Vec<usize>,
}

/// A point in space which the fractures share: a corner of a piece of a
/// user-defined fracture, or the point where two fractures only touch
pub struct PointOnFractures {
    pub fractures: Vec<usize>,
    pub point: Point3<f64>,
}

/// A line in space which two fractures share, before it is divided
pub struct LineOnFractures {
    pub fractures: [usize; 2],
    pub ends: [Point3<f64>; 2],
}

/// The shared points and lines of the network
pub struct SharedNetwork {
    pub shared_points: SharedPoints,
    pub lines: Vec<SharedLine>,
    /// Ids of the shared points on each final fracture
    pub points_of_fracture: Vec<Vec<usize>>,
}

/// Divides the intersections and the common edges of the pieces, and joins
/// the common corners of the pieces and the points where two fractures touch.
///
/// Lines which are on one straight line and overlap (for example an
/// intersection along the common edge of two pieces) are divided together,
/// thus they have the same points.
///
/// # Arguments
///
/// * `planes` - The final fractures in their planes
/// * `groups` - Group of each final fracture: its own position, or the
///   position of the first piece of its user-defined fracture
/// * `intersections` - Intersections between final fractures
/// * `piece_edges` - Common edges of two pieces of one user-defined fracture
/// * `points` - The corners of the pieces, and the points where two fractures touch
pub fn build_shared_network(
    planes: &[FracturePlane],
    groups: &[usize],
    intersections: &[LineOnFractures],
    piece_edges: &[LineOnFractures],
    points: &[PointOnFractures],
) -> SharedNetwork {
    let min_size = planes.first().map_or(1., |plane| plane.min_size);
    let tolerance = MERGE_DISTANCE_PART * min_size;
    let mut network = SharedNetwork {
        shared_points: SharedPoints::new(tolerance),
        lines: Vec::new(),
        points_of_fracture: vec![Vec::new(); planes.len()],
    };

    for point in points {
        let point_groups: Vec<usize> = point
            .fractures
            .iter()
            .map(|&fracture| groups[fracture])
            .collect();
        let id = network
            .shared_points
            .find_or_insert(point.point, &point_groups);
        for &fracture in &point.fractures {
            network.points_of_fracture[fracture].push(id);
        }
    }

    // The intersections first, then the piece edges
    let lines: Vec<&LineOnFractures> = intersections.iter().chain(piece_edges).collect();
    let is_piece_edge = |line: usize| line >= intersections.len();
    let line_size = |line: usize, point: &Point3<f64>| {
        if !is_piece_edge(line) {
            return min_size;
        }
        lines[line]
            .fractures
            .iter()
            .map(|&fracture| {
                let plane = &planes[fracture];
                plane.size_at(plane.frame.to_plane(point))
            })
            .fold(f64::INFINITY, f64::min)
    };

    let mut point_ids_of_line: Vec<Vec<usize>> = vec![Vec::new(); lines.len()];

    for cluster in overlapping_lines(&lines, planes.len(), tolerance) {
        let ids = divide_cluster(
            &cluster,
            &lines,
            planes,
            groups,
            &line_size,
            &mut network.shared_points,
        );
        for (line, point_ids) in cluster.into_iter().zip(ids) {
            point_ids_of_line[line] = point_ids;
        }
    }

    for (line, point_ids) in lines.iter().zip(point_ids_of_line) {
        for fracture in line.fractures {
            network.points_of_fracture[fracture].extend(&point_ids);
        }
        network.lines.push(SharedLine {
            fractures: line.fractures,
            point_ids,
        });
    }

    for points in &mut network.points_of_fracture {
        points.sort_unstable();
        points.dedup();
    }

    network
}

/// Groups of lines which are on one straight line and overlap, on a common
/// fracture. A line which overlaps no other line is a group of one.
fn overlapping_lines(
    lines: &[&LineOnFractures],
    fracture_count: usize,
    tolerance: f64,
) -> Vec<Vec<usize>> {
    let mut lines_of_fracture: Vec<Vec<usize>> = vec![Vec::new(); fracture_count];
    for (index, line) in lines.iter().enumerate() {
        for fracture in line.fractures {
            lines_of_fracture[fracture].push(index);
        }
    }

    let mut parent: Vec<usize> = (0..lines.len()).collect();
    fn root(parent: &mut [usize], mut index: usize) -> usize {
        while parent[index] != index {
            parent[index] = parent[parent[index]];
            index = parent[index];
        }
        index
    }

    for fracture_lines in &lines_of_fracture {
        for (i, &first) in fracture_lines.iter().enumerate() {
            for &second in &fracture_lines[i + 1..] {
                if are_overlapping(lines[first], lines[second], tolerance) {
                    let (first_root, second_root) =
                        (root(&mut parent, first), root(&mut parent, second));
                    parent[first_root] = second_root;
                }
            }
        }
    }

    let mut clusters: HashMap<usize, Vec<usize>> = HashMap::new();
    for index in 0..lines.len() {
        clusters
            .entry(root(&mut parent, index))
            .or_default()
            .push(index);
    }

    let mut clusters: Vec<Vec<usize>> = clusters.into_values().collect();
    clusters.sort_unstable();
    clusters
}

/// True if the two lines are on one straight line and overlap by more than
/// the tolerance
fn are_overlapping(first: &LineOnFractures, second: &LineOnFractures, tolerance: f64) -> bool {
    let [start, end] = first.ends;
    let length = distance(&start, &end);
    if length <= tolerance {
        return false;
    }
    let direction = (end - start) / length;

    let distance_to_line = |point: &Point3<f64>| {
        let offset = point - start;
        (offset - direction * offset.dot(&direction)).norm()
    };
    if second
        .ends
        .iter()
        .any(|point| distance_to_line(point) > tolerance)
    {
        return false;
    }

    let positions = second.ends.map(|point| (point - start).dot(&direction));
    let overlap =
        positions[0].max(positions[1]).min(length) - positions[0].min(positions[1]).max(0.);
    overlap > tolerance
}

/// A point on a line: its position along the line, and the point
type LinePoint = (f64, Point3<f64>);

/// Registers the special points of a group of overlapping lines, and
/// divides the lines between them. Returns the point ids of each line.
fn divide_cluster(
    cluster: &[usize],
    lines: &[&LineOnFractures],
    planes: &[FracturePlane],
    groups: &[usize],
    line_size: &impl Fn(usize, &Point3<f64>) -> f64,
    shared_points: &mut SharedPoints,
) -> Vec<Vec<usize>> {
    let tolerance = shared_points.tolerance;
    let [origin, end] = lines[cluster[0]].ends;
    let direction = (end - origin).normalize();
    let position_of = |point: &Point3<f64>| (point - origin).dot(&direction);

    let mut cluster_groups: Vec<usize> = cluster
        .iter()
        .flat_map(|&line| lines[line].fractures.map(|fracture| groups[fracture]))
        .collect();
    cluster_groups.sort_unstable();
    cluster_groups.dedup();

    let mut special: Vec<LinePoint> = cluster
        .iter()
        .flat_map(|&line| {
            let mut points = special_points(lines[line], planes, lines);
            points.extend(lines[line].ends);
            points
        })
        .map(|point| (position_of(&point), point))
        .collect();
    special.sort_by(|a, b| a.0.total_cmp(&b.0));
    special.dedup_by(|later, earlier| later.0 - earlier.0 < tolerance);

    // The ranges of the lines along the cluster
    let ranges: Vec<[f64; 2]> = cluster
        .iter()
        .map(|&line| {
            let [first, second] = lines[line].ends.map(|point| position_of(&point));
            [first.min(second), first.max(second)]
        })
        .collect();

    let mut points: Vec<(f64, usize)> = Vec::new();
    for (position, point) in special {
        let id = shared_points.find_or_insert(point, &cluster_groups);

        if let Some(&(previous_position, previous)) = points.last() {
            if previous == id {
                continue;
            }
            let middle = (previous_position + position) / 2.;
            let covering: Vec<usize> = cluster
                .iter()
                .zip(&ranges)
                .filter(|(_, [low, high])| (*low..=*high).contains(&middle))
                .map(|(&line, _)| line)
                .collect();
            let size = |point: &Point3<f64>| {
                covering
                    .iter()
                    .map(|&line| line_size(line, point))
                    .fold(f64::INFINITY, f64::min)
            };
            let start = shared_points.points[previous];
            for between in divide(&start, &point, &size) {
                points.push((position_of(&between), shared_points.add(between)));
            }
        }

        points.push((position, id));
    }

    cluster
        .iter()
        .zip(&ranges)
        .map(|(&line, [low, high])| {
            let mut ids: Vec<usize> = points
                .iter()
                .filter(|(position, _)| {
                    *low - tolerance <= *position && *position <= *high + tolerance
                })
                .map(|&(_, id)| id)
                .collect();
            let [first, second] = lines[line].ends.map(|point| position_of(&point));
            if first > second {
                ids.reverse();
            }
            ids
        })
        .collect()
}

/// Points between the two points (not the two points), so that no part is
/// longer than the size at its middle (by halving). Without a size (no line
/// covers the part), no points.
pub fn divide(
    start: &Point3<f64>,
    end: &Point3<f64>,
    size: &impl Fn(&Point3<f64>) -> f64,
) -> Vec<Point3<f64>> {
    let middle = Point3::from((start.coords + end.coords) / 2.);

    if distance(start, end) <= size(&middle) {
        return Vec::new();
    }

    let mut points = divide(start, &middle, size);
    points.push(middle);
    points.extend(divide(&middle, end, size));
    points
}

/// The points of a line where something on one of its two fractures
/// crosses or touches it: a boundary line or another shared line (a triple
/// point, or the end of an intersection on a common edge of two pieces), or
/// a corner of the fracture on it
fn special_points(
    line: &LineOnFractures,
    planes: &[FracturePlane],
    lines: &[&LineOnFractures],
) -> Vec<Point3<f64>> {
    let mut positions = Vec::new();

    for fracture in line.fractures {
        let plane = &planes[fracture];
        let line_2d = line.ends.map(|end| plane.frame.to_plane(&end));
        let tolerance = MERGE_DISTANCE_PART * plane.min_size;

        let other_lines = lines
            .iter()
            .filter(|other| other.fractures.contains(&fracture))
            .filter(|other| !std::ptr::eq(**other, line))
            .map(|other| other.ends.map(|end| plane.frame.to_plane(&end)));
        let boundary_lines = plane
            .boundary_segments
            .iter()
            .map(|boundary| boundary.segment);

        for other in other_lines.chain(boundary_lines) {
            if let Some((position, _)) = segment_crossing(line_2d, other) {
                positions.push(position);
            }
            for end in other {
                if distance_to_segment(end, line_2d) < tolerance {
                    positions.push(position_on_segment(end, line_2d));
                }
            }
        }

        for corner in &plane.outline {
            if distance_to_segment(*corner, line_2d) < tolerance {
                positions.push(position_on_segment(*corner, line_2d));
            }
        }
    }

    let [start, end] = line.ends;
    positions
        .into_iter()
        .map(|position| start + (end - start) * position)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn line(start: [f64; 3], end: [f64; 3]) -> LineOnFractures {
        LineOnFractures {
            fractures: [0, 1],
            ends: [Point3::from(start), Point3::from(end)],
        }
    }

    #[test]
    fn test_collinear_lines_which_overlap_are_overlapping() {
        let first = line([0., 0., 0.], [4., 0., 0.]);
        let second = line([6., 0., 0.], [2., 0., 0.]);

        assert!(are_overlapping(&first, &second, 1e-9));
    }

    #[test]
    fn test_collinear_lines_which_only_touch_are_not_overlapping() {
        let first = line([0., 0., 0.], [4., 0., 0.]);
        let second = line([4., 0., 0.], [8., 0., 0.]);

        assert!(!are_overlapping(&first, &second, 1e-9));
    }

    #[test]
    fn test_parallel_lines_are_not_overlapping() {
        let first = line([0., 0., 0.], [4., 0., 0.]);
        let second = line([0., 1., 0.], [4., 1., 0.]);

        assert!(!are_overlapping(&first, &second, 1e-9));
    }

    #[test]
    fn test_divide_gives_parts_not_longer_than_the_size() {
        let start = Point3::new(0., 0., 0.);
        let end = Point3::new(10., 0., 0.);

        let points = divide(&start, &end, &|_| 3.);

        // Halving: 10 → 5 → 2.5
        assert_eq!(points.len(), 3);
        assert_eq!(points[1], Point3::new(5., 0., 0.));
    }

    #[test]
    fn test_points_of_one_group_merge() {
        let mut shared_points = SharedPoints::new(0.1);

        let first = shared_points.find_or_insert(Point3::new(1., 1., 1.), &[0, 1]);
        let second = shared_points.find_or_insert(Point3::new(1.05, 1., 1.), &[1, 2]);

        assert_eq!(first, second);
    }

    #[test]
    fn test_points_without_a_common_group_do_not_merge() {
        let mut shared_points = SharedPoints::new(0.1);

        let first = shared_points.find_or_insert(Point3::new(1., 1., 1.), &[0, 1]);
        let second = shared_points.find_or_insert(Point3::new(1.05, 1., 1.), &[2, 3]);

        assert_ne!(first, second);
    }
}
