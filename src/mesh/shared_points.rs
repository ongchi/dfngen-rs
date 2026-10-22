//! The points which two or more fractures share: the points of the
//! intersections, and the points on the common edges and corners of the
//! pieces of a user-defined fracture
//!
//! Each shared line is divided once, in space. The fractures on both sides
//! use the same points, thus the mesh is connected across the line.

use std::collections::HashMap;

use parry3d_f64::na::{center, distance, Point3};

use super::plane::{segment_crossing, FracturePlane, Point2};
use super::{JOIN_DISTANCE_PART, MERGE_DISTANCE_PART};
use crate::region_boundaries::{distance as distance_2d, distance_to_segment, position_on_segment};

/// How a point merges with the other points. Points of user-defined
/// fractures which FRAM did not check can be nearer to each other than
/// `min_feature_size`: their tiny features are joined.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Join {
    /// Merges only within the rounding tolerance, and other points do not
    /// join it
    None,
    /// Does not move, but a `Joinable` point can join it (a corner of a piece
    /// of a user-defined fracture)
    Target,
    /// Joins a `Target` or `Joinable` point within the join distance
    Joinable,
}

/// Points in space which two or more fractures share. A point with a
/// "group" can merge with another point of the same group; the group of a
/// fracture is its own position, or the first piece of its user-defined
/// fracture.
pub struct SharedPoints {
    pub points: Vec<Point3<f64>>,
    /// Groups of each point (empty for a point which does not merge)
    groups: Vec<Vec<usize>>,
    tolerance: f64,
    join_distance: f64,
    /// Points which can merge, by the grid cell of size `tolerance` which has them
    cells: HashMap<[i64; 3], Vec<usize>>,
    /// The join targets, by the grid cell of size `join_distance` which has them
    join_cells: HashMap<[i64; 3], Vec<usize>>,
}

impl SharedPoints {
    fn new(tolerance: f64, join_distance: f64) -> Self {
        Self {
            points: Vec::new(),
            groups: Vec::new(),
            tolerance,
            join_distance,
            cells: HashMap::new(),
            join_cells: HashMap::new(),
        }
    }

    fn cell_of(point: &Point3<f64>, cell_size: f64) -> [i64; 3] {
        [point.x, point.y, point.z].map(|coordinate| (coordinate / cell_size).floor() as i64)
    }

    /// The nearest point within the distance (`cell_size` of the grid) which
    /// has one of the groups
    fn nearest(
        &self,
        point: &Point3<f64>,
        groups: &[usize],
        cells: &HashMap<[i64; 3], Vec<usize>>,
        max_distance: f64,
    ) -> Option<usize> {
        let [x, y, z] = Self::cell_of(point, max_distance);
        (x - 1..=x + 1)
            .flat_map(|cell_x| (y - 1..=y + 1).map(move |cell_y| (cell_x, cell_y)))
            .flat_map(|(cell_x, cell_y)| {
                (z - 1..=z + 1).map(move |cell_z| [cell_x, cell_y, cell_z])
            })
            .filter_map(|cell| cells.get(&cell))
            .flatten()
            .copied()
            .filter(|&id| distance(&self.points[id], point) < max_distance)
            .filter(|&id| self.groups[id].iter().any(|group| groups.contains(group)))
            .min_by(|&a, &b| {
                distance(&self.points[a], point).total_cmp(&distance(&self.points[b], point))
            })
    }

    /// The point nearer than the tolerance (or the join distance, for a
    /// joinable point) which has one of the groups, or a new point. The point
    /// gets the groups.
    fn find_or_insert(&mut self, point: Point3<f64>, groups: &[usize], join: Join) -> usize {
        let found = self
            .nearest(&point, groups, &self.cells, self.tolerance)
            .or_else(|| {
                (join == Join::Joinable)
                    .then(|| self.nearest(&point, groups, &self.join_cells, self.join_distance))
                    .flatten()
            });

        if let Some(id) = found {
            for &group in groups {
                if !self.groups[id].contains(&group) {
                    self.groups[id].push(group);
                }
            }
            return id;
        }

        let id = self.add(point);
        self.groups[id] = groups.to_vec();
        self.cells
            .entry(Self::cell_of(&point, self.tolerance))
            .or_default()
            .push(id);
        if join != Join::None {
            self.join_cells
                .entry(Self::cell_of(&point, self.join_distance))
                .or_default()
                .push(id);
        }
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
    /// `Joinable`: an outline edge near the line can use its points
    pub join: Join,
}

/// A point in space which the fractures share: a corner of a piece of a
/// user-defined fracture, or the point where two fractures only touch
pub struct PointOnFractures {
    pub fractures: Vec<usize>,
    pub point: Point3<f64>,
    pub join: Join,
}

/// A line in space which two fractures share, before it is divided
pub struct LineOnFractures {
    pub fractures: [usize; 2],
    pub ends: [Point3<f64>; 2],
    /// `Joinable` for two user-defined fractures which FRAM did not check,
    /// else `None`
    pub join: Join,
}

/// The shared points and lines of the network
pub struct SharedNetwork {
    pub shared_points: SharedPoints,
    pub lines: Vec<SharedLine>,
    /// Ids of the shared points on each final fracture
    pub points_of_fracture: Vec<Vec<usize>>,
}

impl SharedNetwork {
    /// Makes each shared point of a fracture which is on a shared line of
    /// that fracture (in its plane, nearer than the merge distance to an
    /// edge of the line, not at an end of the edge) a point of the line.
    /// Else a triangle of the fracture has its three points on one line.
    /// The point then is a point of the other fracture of the line too, thus
    /// this repeats until no point is inserted. Returns the fractures of the
    /// changed lines.
    pub fn insert_points_on_lines(&mut self, planes: &[FracturePlane]) -> Vec<usize> {
        let mut changed_fractures = Vec::new();

        loop {
            for points in &mut self.points_of_fracture {
                points.sort_unstable();
                points.dedup();
            }

            let insertions = self.points_on_lines(planes);
            if insertions.is_empty() {
                break;
            }

            for (line_index, mut points) in insertions {
                let line = &mut self.lines[line_index];
                // By the edge, then along the edge
                points.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.total_cmp(&b.1)));
                points.dedup_by_key(|(_, _, id)| *id);

                let mut point_ids = Vec::with_capacity(line.point_ids.len() + points.len());
                let mut next = points.iter().peekable();
                for (edge, &id) in line.point_ids.iter().enumerate() {
                    point_ids.push(id);
                    while let Some(&(_, _, inserted)) =
                        next.next_if(|(point_edge, _, _)| *point_edge == edge)
                    {
                        point_ids.push(inserted);
                    }
                }
                line.point_ids = point_ids;

                for fracture in line.fractures {
                    self.points_of_fracture[fracture].extend(points.iter().map(|&(_, _, id)| id));
                    changed_fractures.push(fracture);
                }
            }
        }

        changed_fractures.sort_unstable();
        changed_fractures.dedup();
        changed_fractures
    }

    /// The shared points on the edges of each line (see
    /// `insert_points_on_lines`): for each line, (edge, position on the edge,
    /// point id); edge k goes from point k to point k + 1 of the line
    fn points_on_lines(
        &self,
        planes: &[FracturePlane],
    ) -> HashMap<usize, Vec<(usize, f64, usize)>> {
        let tolerance = self.shared_points.tolerance;
        let mut lines_of_fracture: Vec<Vec<usize>> = vec![Vec::new(); planes.len()];
        for (index, line) in self.lines.iter().enumerate() {
            for fracture in line.fractures {
                lines_of_fracture[fracture].push(index);
            }
        }

        let mut insertions: HashMap<usize, Vec<(usize, f64, usize)>> = HashMap::new();
        for (fracture, plane) in planes.iter().enumerate() {
            if lines_of_fracture[fracture].is_empty() {
                continue;
            }

            // The shared points of the fracture in its plane, by grid cell
            let cell_size = plane.min_size;
            let cell_of = |[x, y]: Point2| {
                (
                    (x / cell_size).floor() as i64,
                    (y / cell_size).floor() as i64,
                )
            };
            let mut cells: HashMap<(i64, i64), Vec<(usize, Point2)>> = HashMap::new();
            for &id in &self.points_of_fracture[fracture] {
                let position = plane.frame.to_plane(&self.shared_points.points[id]);
                cells
                    .entry(cell_of(position))
                    .or_default()
                    .push((id, position));
            }

            for &line_index in &lines_of_fracture[fracture] {
                let point_ids = &self.lines[line_index].point_ids;
                for (edge, pair) in point_ids.windows(2).enumerate() {
                    let segment = [pair[0], pair[1]]
                        .map(|id| plane.frame.to_plane(&self.shared_points.points[id]));
                    let (low, high) = (
                        cell_of([
                            segment[0][0].min(segment[1][0]) - tolerance,
                            segment[0][1].min(segment[1][1]) - tolerance,
                        ]),
                        cell_of([
                            segment[0][0].max(segment[1][0]) + tolerance,
                            segment[0][1].max(segment[1][1]) + tolerance,
                        ]),
                    );

                    for cell_x in low.0..=high.0 {
                        for cell_y in low.1..=high.1 {
                            let Some(points) = cells.get(&(cell_x, cell_y)) else {
                                continue;
                            };
                            for &(id, position) in points {
                                let is_on_edge = !pair.contains(&id)
                                    && distance_to_segment(position, segment) < tolerance
                                    && segment
                                        .iter()
                                        .all(|&end| distance_2d(end, position) > tolerance)
                                    && !point_ids.contains(&id);
                                if is_on_edge {
                                    insertions.entry(line_index).or_default().push((
                                        edge,
                                        position_on_segment(position, segment),
                                        id,
                                    ));
                                }
                            }
                        }
                    }
                }
            }
        }

        insertions
    }

    /// Splits each edge (two consecutive points of a shared line) at its
    /// middle, in each line which has it. Returns the fractures of these lines.
    pub fn split_edges(&mut self, edges: &[[usize; 2]]) -> Vec<usize> {
        // The middle point of each edge, by the edge with the smaller id first
        let mut middles: HashMap<(usize, usize), usize> = HashMap::new();
        for &[first, second] in edges {
            let points = &self.shared_points.points;
            let middle = center(&points[first], &points[second]);
            let id = self.shared_points.add(middle);
            middles.insert((first.min(second), first.max(second)), id);
        }

        let mut changed_fractures = Vec::new();
        for line in &mut self.lines {
            let mut point_ids = Vec::with_capacity(line.point_ids.len());
            let mut is_changed = false;

            for (index, &id) in line.point_ids.iter().enumerate() {
                point_ids.push(id);
                let Some(&next) = line.point_ids.get(index + 1) else {
                    continue;
                };
                if let Some(&middle) = middles.get(&(id.min(next), id.max(next))) {
                    point_ids.push(middle);
                    for fracture in line.fractures {
                        self.points_of_fracture[fracture].push(middle);
                    }
                    is_changed = true;
                }
            }

            if is_changed {
                line.point_ids = point_ids;
                changed_fractures.extend(line.fractures);
            }
        }

        changed_fractures.sort_unstable();
        changed_fractures.dedup();
        changed_fractures
    }
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
        shared_points: SharedPoints::new(tolerance, JOIN_DISTANCE_PART * min_size),
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
            .find_or_insert(point.point, &point_groups, point.join);
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

    let join_distance = network.shared_points.join_distance;
    for cluster in overlapping_lines(&lines, planes.len(), [tolerance, join_distance]) {
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
            join: line.join,
        });
    }

    network.insert_points_on_lines(planes);

    network
}

/// Groups of lines which are on one straight line and overlap, on a common
/// fracture. A line which overlaps no other line is a group of one.
///
/// Two `Joinable` lines are one group also if they are only near each other
/// (within the join distance), almost parallel and overlap: they then get
/// the same points.
fn overlapping_lines(
    lines: &[&LineOnFractures],
    fracture_count: usize,
    [tolerance, join_distance]: [f64; 2],
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
                let is_joinable =
                    lines[first].join == Join::Joinable && lines[second].join == Join::Joinable;
                let distance = if is_joinable {
                    join_distance
                } else {
                    tolerance
                };
                if are_overlapping(lines[first], lines[second], distance) {
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
    // Tiny features of user-defined fractures which FRAM did not check are
    // joined: their points within the join distance are one point
    let join = if cluster
        .iter()
        .all(|&line| lines[line].join == Join::Joinable)
    {
        Join::Joinable
    } else {
        Join::None
    };
    let merge_distance = match join {
        Join::Joinable => shared_points.join_distance,
        _ => tolerance,
    };
    special.sort_by(|a, b| a.0.total_cmp(&b.0));
    special.dedup_by(|later, earlier| later.0 - earlier.0 < merge_distance);

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
        let id = shared_points.find_or_insert(point, &cluster_groups, join);

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
            join: Join::None,
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
        let mut shared_points = SharedPoints::new(0.1, 1.);

        let first = shared_points.find_or_insert(Point3::new(1., 1., 1.), &[0, 1], Join::None);
        let second = shared_points.find_or_insert(Point3::new(1.05, 1., 1.), &[1, 2], Join::None);

        assert_eq!(first, second);
    }

    #[test]
    fn test_joinable_point_joins_a_target_within_the_join_distance() {
        let mut shared_points = SharedPoints::new(0.1, 1.);

        let corner = shared_points.find_or_insert(Point3::new(1., 1., 1.), &[0], Join::Target);
        let end = shared_points.find_or_insert(Point3::new(1.5, 1., 1.), &[0, 1], Join::Joinable);

        assert_eq!(corner, end);
        assert_eq!(shared_points.points[corner], Point3::new(1., 1., 1.));
    }

    #[test]
    fn test_point_which_is_not_joinable_does_not_join_within_the_join_distance() {
        let mut shared_points = SharedPoints::new(0.1, 1.);

        let first = shared_points.find_or_insert(Point3::new(1., 1., 1.), &[0], Join::Joinable);
        let second = shared_points.find_or_insert(Point3::new(1.5, 1., 1.), &[0], Join::None);

        assert_ne!(first, second);
    }

    #[test]
    fn test_joinable_point_does_not_join_a_point_which_is_no_target() {
        let mut shared_points = SharedPoints::new(0.1, 1.);

        let first = shared_points.find_or_insert(Point3::new(1., 1., 1.), &[0], Join::None);
        let second = shared_points.find_or_insert(Point3::new(1.5, 1., 1.), &[0], Join::Joinable);

        assert_ne!(first, second);
    }

    #[test]
    fn test_points_without_a_common_group_do_not_merge() {
        let mut shared_points = SharedPoints::new(0.1, 1.);

        let first = shared_points.find_or_insert(Point3::new(1., 1., 1.), &[0, 1], Join::None);
        let second = shared_points.find_or_insert(Point3::new(1.05, 1., 1.), &[2, 3], Join::None);

        assert_ne!(first, second);
    }
}
