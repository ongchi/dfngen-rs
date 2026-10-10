use parry3d_f64::bounding_volume::{Aabb, BoundingVolume};
use parry3d_f64::na::{Point3, Vector3};
use parry3d_f64::query::details::clip_halfspace_polygon;
use rand::Rng;

/// Tolerance of the point test, as a part of the diagonal of the mesh bounding box.
/// A point on the face between two cells must be in one of them.
const RELATIVE_TOLERANCE: f64 = 1e-9;

/// One tetrahedron of a mesh cell, in the coordinates of the mesh file
pub struct MeshTetrahedron {
    pub vertices: [Point3<f64>; 4],
    /// Value of the `region` cell property. 0: the cell has no region.
    pub region: u64,
}

/// The cells of the mesh which have the same region value above 0
#[derive(Debug, Clone, PartialEq)]
pub struct MeshRegion {
    /// Value of the `region` cell property
    pub value: u64,
    /// Bounding box of the cells, with the mesh center at the origin
    pub bounds: Aabb,
    /// Sum of the cell volumes
    pub volume: f64,
}

#[derive(Debug)]
struct Tetrahedron {
    vertices: [Point3<f64>; 4],
    /// A point on each face, and the outward unit normal of the face
    faces: [(Point3<f64>, Vector3<f64>); 4],
    /// 0: the cell has no region. N > 0: region N of `DomainMesh::regions`.
    region_id: usize,
}

impl Tetrahedron {
    fn new(vertices: &[Point3<f64>; 4], region_id: usize) -> Self {
        let faces = std::array::from_fn(|i| {
            let opposite_vertex = vertices[i];
            let [p, q, r] = std::array::from_fn(|k| vertices[(i + 1 + k) % 4]);
            let normal = (q - p).cross(&(r - p)).normalize();

            if normal.dot(&(opposite_vertex - p)) > 0. {
                (p, -normal)
            } else {
                (p, normal)
            }
        });

        Self {
            vertices: *vertices,
            faces,
            region_id,
        }
    }

    /// Point of the tetrahedron for 3 numbers in [0, 1). Uniform random
    /// numbers give a uniform random point.
    fn point_at(&self, mut s: f64, mut t: f64, mut u: f64) -> Point3<f64> {
        // Fold the unit cube into the tetrahedron s + t + u <= 1
        // (Rocchini and Cignoni, Generating random points in a tetrahedron, 2000)
        if s + t > 1. {
            s = 1. - s;
            t = 1. - t;
        }
        if t + u > 1. {
            let old_u = u;
            u = 1. - s - t;
            t = 1. - old_u;
        } else if s + t + u > 1. {
            let old_u = u;
            u = s + t + u - 1.;
            s = 1. - t - old_u;
        }

        let [a, b, c, d] = self.vertices;
        a + (b - a) * s + (c - a) * t + (d - a) * u
    }

    fn contains_point(&self, point: &Point3<f64>, tolerance: f64) -> bool {
        self.faces
            .iter()
            .all(|(face_point, normal)| (point - face_point).dot(normal) <= tolerance)
    }

    fn intersects_polygon(&self, vertices: &[Point3<f64>]) -> bool {
        let mut part_in_tetrahedron = vertices.to_vec();
        let mut clipped = Vec::with_capacity(vertices.len() + 4);

        for (face_point, normal) in &self.faces {
            clip_halfspace_polygon(face_point, normal, &part_in_tetrahedron, &mut clipped);
            std::mem::swap(&mut part_in_tetrahedron, &mut clipped);
        }

        !part_in_tetrahedron.is_empty()
    }
}

fn tetrahedron_volume(vertices: &[Point3<f64>; 4]) -> f64 {
    let [a, b, c, d] = vertices;
    ((b - a).cross(&(c - a)).dot(&(d - a)) / 6.).abs()
}

/// The tetrahedra of one region with their cumulative volumes, to select a
/// tetrahedron with a probability in proportion to its volume
#[derive(Debug, Default, Clone)]
struct VolumeTable {
    tetrahedra: Vec<usize>,
    cumulative_volumes: Vec<f64>,
}

impl VolumeTable {
    fn push(&mut self, tetrahedron: usize, volume: f64) {
        let total_volume = self.cumulative_volumes.last().copied().unwrap_or(0.);
        self.tetrahedra.push(tetrahedron);
        self.cumulative_volumes.push(total_volume + volume);
    }

    /// The tetrahedron at a fraction in [0, 1) of the total volume
    fn select(&self, fraction: f64) -> usize {
        let total_volume = self.cumulative_volumes[self.cumulative_volumes.len() - 1];
        let position = self
            .cumulative_volumes
            .partition_point(|&volume| volume <= fraction * total_volume);

        self.tetrahedra[position.min(self.tetrahedra.len() - 1)]
    }
}

/// Uniform grid over the bounding box of the mesh. Each grid cell has the
/// indices of the tetrahedra which can have a point in that grid cell.
#[derive(Debug)]
struct Grid {
    half_size: Vector3<f64>,
    cell_size: Vector3<f64>,
    cell_counts: [usize; 3],
    tetrahedra_of_cell: Vec<Vec<usize>>,
}

impl Grid {
    /// Makes an empty grid with about one grid cell for each tetrahedron.
    fn new(size: &Vector3<f64>, tetrahedron_count: usize) -> Self {
        let cell_edge = (size.x * size.y * size.z / tetrahedron_count as f64).cbrt();
        let cell_counts =
            std::array::from_fn(|axis| ((size[axis] / cell_edge).ceil() as usize).max(1));
        let cell_size = Vector3::from_fn(|axis, _| size[axis] / cell_counts[axis] as f64);

        Self {
            half_size: size * 0.5,
            cell_size,
            cell_counts,
            tetrahedra_of_cell: vec![Vec::new(); cell_counts.iter().product()],
        }
    }

    /// Grid cell of a point. A point outside the grid gives the nearest grid cell.
    fn cell_of(&self, point: &Point3<f64>) -> [usize; 3] {
        std::array::from_fn(|axis| {
            // A negative number becomes 0 in the conversion
            let cell = ((point[axis] + self.half_size[axis]) / self.cell_size[axis]) as usize;
            cell.min(self.cell_counts[axis] - 1)
        })
    }

    fn index_of(&self, cell: [usize; 3]) -> usize {
        (cell[2] * self.cell_counts[1] + cell[1]) * self.cell_counts[0] + cell[0]
    }

    /// Indices in `tetrahedra_of_cell` of the grid cells which a box touches
    fn cells_in(&self, bounds: &Aabb) -> impl Iterator<Item = usize> + '_ {
        let first = self.cell_of(&bounds.mins);
        let last = self.cell_of(&bounds.maxs);

        (first[2]..=last[2]).flat_map(move |z| {
            (first[1]..=last[1])
                .flat_map(move |y| (first[0]..=last[0]).map(move |x| self.index_of([x, y, z])))
        })
    }

    fn insert(&mut self, tetrahedron: usize, bounds: &Aabb) {
        let cells: Vec<usize> = self.cells_in(bounds).collect();

        for cell in cells {
            self.tetrahedra_of_cell[cell].push(tetrahedron);
        }
    }
}

/// Domain which is given by the cells of a mesh
///
/// The bounding box of the mesh is the box domain of the generation. The
/// generation uses coordinates with the center of this box at the origin, and
/// so do the functions of this struct. `center` is the center in the
/// coordinates of the mesh file.
#[derive(Debug)]
pub struct DomainMesh {
    center: Vector3<f64>,
    size: Vector3<f64>,
    volume: f64,
    regions: Vec<MeshRegion>,
    tetrahedra: Vec<Tetrahedron>,
    /// Table 0 has all tetrahedra, table N > 0 those of region N
    volume_tables: Vec<VolumeTable>,
    grid: Grid,
    tolerance: f64,
}

impl DomainMesh {
    /// Makes the domain from the tetrahedra of the mesh cells.
    /// Tetrahedra with no volume are ignored.
    pub fn new(mesh_tetrahedra: &[MeshTetrahedron]) -> Result<Self, String> {
        let mesh_tetrahedra: Vec<&MeshTetrahedron> = mesh_tetrahedra
            .iter()
            .filter(|tetrahedron| tetrahedron_volume(&tetrahedron.vertices) > 0.)
            .collect();

        if mesh_tetrahedra.is_empty() {
            return Err("the mesh has no cell with a volume".to_string());
        }

        let mesh_bounds = Aabb::from_points(
            mesh_tetrahedra
                .iter()
                .flat_map(|tetrahedron| &tetrahedron.vertices),
        );
        let center = mesh_bounds.center().coords;
        let size = mesh_bounds.extents();
        let tolerance = RELATIVE_TOLERANCE * size.norm();

        let mut region_values: Vec<u64> = mesh_tetrahedra
            .iter()
            .map(|tetrahedron| tetrahedron.region)
            .filter(|&region| region > 0)
            .collect();
        region_values.sort();
        region_values.dedup();

        let mut regions: Vec<MeshRegion> = region_values
            .iter()
            .map(|&value| MeshRegion {
                value,
                bounds: Aabb::new_invalid(),
                volume: 0.,
            })
            .collect();

        let mut volume = 0.;
        let mut tetrahedra = Vec::with_capacity(mesh_tetrahedra.len());
        let mut volume_tables = vec![VolumeTable::default(); regions.len() + 1];
        let mut grid = Grid::new(&size, mesh_tetrahedra.len());

        for (i, mesh_tetrahedron) in mesh_tetrahedra.iter().enumerate() {
            let vertices = mesh_tetrahedron.vertices.map(|vertex| vertex - center);
            let bounds = Aabb::from_points(&vertices);
            let tetrahedron_volume = tetrahedron_volume(&vertices);
            let region_id = region_values
                .binary_search(&mesh_tetrahedron.region)
                .map_or(0, |position| position + 1);

            volume += tetrahedron_volume;
            volume_tables[0].push(i, tetrahedron_volume);

            if region_id > 0 {
                let region = &mut regions[region_id - 1];
                region.bounds.merge(&bounds);
                region.volume += tetrahedron_volume;
                volume_tables[region_id].push(i, tetrahedron_volume);
            }

            grid.insert(i, &bounds.loosened(tolerance));
            tetrahedra.push(Tetrahedron::new(&vertices, region_id));
        }

        Ok(Self {
            center,
            size,
            volume,
            regions,
            tetrahedra,
            volume_tables,
            grid,
            tolerance,
        })
    }

    /// Center of the bounding box, in the coordinates of the mesh file
    pub fn center(&self) -> Vector3<f64> {
        self.center
    }

    /// Size of the bounding box in x, y, z
    pub fn size(&self) -> Vector3<f64> {
        self.size
    }

    /// Sum of the cell volumes
    pub fn volume(&self) -> f64 {
        self.volume
    }

    /// The regions of the mesh, in the order of their values.
    /// Region N of a tetrahedron or a family is `regions()[N - 1]`.
    pub fn regions(&self) -> &[MeshRegion] {
        &self.regions
    }

    /// True if the point is in a cell of the region.
    /// `region_id` 0 means a cell of any region.
    pub fn contains_point(&self, point: &Point3<f64>, region_id: usize) -> bool {
        let cell = self.grid.index_of(self.grid.cell_of(point));

        self.grid.tetrahedra_of_cell[cell].iter().any(|&i| {
            let tetrahedron = &self.tetrahedra[i];
            (region_id == 0 || tetrahedron.region_id == region_id)
                && tetrahedron.contains_point(point, self.tolerance)
        })
    }

    /// Uniform random point in the cells of the region.
    /// `region_id` 0 means the cells of all regions.
    pub fn sample_point(&self, region_id: usize, rng: &mut impl Rng) -> Point3<f64> {
        let tetrahedron = &self.tetrahedra[self.volume_tables[region_id].select(rng.random())];
        let point = tetrahedron.point_at(rng.random(), rng.random(), rng.random());
        debug_assert!(self.contains_point(&point, region_id));

        point
    }

    /// True if the convex polygon has a point in a cell.
    pub fn intersects_polygon(&self, vertices: &[Point3<f64>]) -> bool {
        let mut near_tetrahedra: Vec<usize> = self
            .grid
            .cells_in(&Aabb::from_points(vertices))
            .flat_map(|cell| &self.grid.tetrahedra_of_cell[cell])
            .copied()
            .collect();
        near_tetrahedra.sort();
        near_tetrahedra.dedup();

        near_tetrahedra
            .iter()
            .any(|&i| self.tetrahedra[i].intersects_polygon(vertices))
    }
}

#[cfg(test)]
pub mod test_meshes {
    use super::*;

    /// The 6 tetrahedra of a box with the given region
    pub fn box_tetrahedra(mins: [f64; 3], maxs: [f64; 3], region: u64) -> Vec<MeshTetrahedron> {
        let corner = |x: usize, y: usize, z: usize| {
            let select = |axis: usize, side: usize| [mins[axis], maxs[axis]][side];
            Point3::new(select(0, x), select(1, y), select(2, z))
        };
        // The 6 tetrahedra have the diagonal from corner 000 to corner 111
        let ring = [
            corner(1, 0, 0),
            corner(1, 1, 0),
            corner(0, 1, 0),
            corner(0, 1, 1),
            corner(0, 0, 1),
            corner(1, 0, 1),
        ];

        (0..6)
            .map(|i| MeshTetrahedron {
                vertices: [corner(0, 0, 0), ring[i], ring[(i + 1) % 6], corner(1, 1, 1)],
                region,
            })
            .collect()
    }

    /// An L shape in the x-y plane, 10 high in z, with its corner at (10, 20, 30):
    /// two boxes of 5 x 5 with region 1 along x, and one box with region 2 in +y.
    /// The bounding box is 10 x 10 x 10 with its center at (15, 25, 35).
    pub fn l_shaped_mesh() -> DomainMesh {
        let mut tetrahedra = box_tetrahedra([10., 20., 30.], [15., 25., 40.], 1);
        tetrahedra.extend(box_tetrahedra([15., 20., 30.], [20., 25., 40.], 1));
        tetrahedra.extend(box_tetrahedra([10., 25., 30.], [15., 30., 40.], 2));

        DomainMesh::new(&tetrahedra).unwrap()
    }
}

#[cfg(test)]
mod tests {
    use rand_mt::Mt64;

    use super::test_meshes::*;
    use super::*;

    fn assert_close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() < 1e-9,
            "{} is not {}",
            actual,
            expected
        );
    }

    /// A square with sides of 2 in a plane of constant z
    fn square(center: [f64; 3]) -> Vec<Point3<f64>> {
        [(-1., -1.), (1., -1.), (1., 1.), (-1., 1.)]
            .iter()
            .map(|(x, y)| Point3::new(center[0] + x, center[1] + y, center[2]))
            .collect()
    }

    #[test]
    fn test_center_is_the_center_of_the_bounding_box() {
        assert_eq!(l_shaped_mesh().center(), Vector3::new(15., 25., 35.));
    }

    #[test]
    fn test_size_is_the_size_of_the_bounding_box() {
        assert_eq!(l_shaped_mesh().size(), Vector3::new(10., 10., 10.));
    }

    #[test]
    fn test_volume_is_the_sum_of_the_cell_volumes() {
        assert_close(l_shaped_mesh().volume(), 750.);
    }

    #[test]
    fn test_regions_are_in_the_order_of_their_values() {
        let mesh = l_shaped_mesh();
        let values: Vec<u64> = mesh.regions().iter().map(|region| region.value).collect();

        assert_eq!(values, [1, 2]);
    }

    #[test]
    fn test_region_volume_is_the_sum_of_its_cell_volumes() {
        let mesh = l_shaped_mesh();

        assert_close(mesh.regions()[0].volume, 500.);
        assert_close(mesh.regions()[1].volume, 250.);
    }

    #[test]
    fn test_region_bounds_have_the_mesh_center_at_the_origin() {
        let mesh = l_shaped_mesh();

        assert_eq!(
            mesh.regions()[1].bounds,
            Aabb::new(Point3::new(-5., 0., -5.), Point3::new(0., 5., 5.))
        );
    }

    #[test]
    fn test_mesh_without_regions_has_no_region() {
        let mesh = DomainMesh::new(&box_tetrahedra([0.; 3], [1.; 3], 0)).unwrap();

        assert!(mesh.regions().is_empty());
    }

    #[test]
    fn test_mesh_without_volume_is_rejected() {
        let flat_box = box_tetrahedra([0.; 3], [1., 1., 0.], 0);

        assert!(DomainMesh::new(&flat_box).is_err());
    }

    #[test]
    fn test_point_in_a_cell_is_in_the_mesh() {
        assert!(l_shaped_mesh().contains_point(&Point3::new(2.5, -2.5, 0.), 0));
    }

    #[test]
    fn test_point_in_the_empty_corner_of_the_bounding_box_is_not_in_the_mesh() {
        assert!(!l_shaped_mesh().contains_point(&Point3::new(2.5, 2.5, 0.), 0));
    }

    #[test]
    fn test_point_outside_the_bounding_box_is_not_in_the_mesh() {
        assert!(!l_shaped_mesh().contains_point(&Point3::new(-6., -2.5, 0.), 0));
    }

    #[test]
    fn test_point_on_the_face_between_two_cells_is_in_the_mesh() {
        assert!(l_shaped_mesh().contains_point(&Point3::new(0., -2.5, 1.), 0));
    }

    #[test]
    fn test_point_is_in_the_region_of_its_cell() {
        assert!(l_shaped_mesh().contains_point(&Point3::new(-2.5, 2.5, 0.), 2));
    }

    #[test]
    fn test_point_is_not_in_the_region_of_another_cell() {
        assert!(!l_shaped_mesh().contains_point(&Point3::new(-2.5, 2.5, 0.), 1));
    }

    fn sample_points(region_id: usize, count: usize) -> Vec<Point3<f64>> {
        let mesh = l_shaped_mesh();
        let mut rng = Mt64::new(1);

        (0..count)
            .map(|_| mesh.sample_point(region_id, &mut rng))
            .collect()
    }

    #[test]
    fn test_sampled_points_are_in_the_mesh() {
        let mesh = l_shaped_mesh();

        assert!(sample_points(0, 1000)
            .iter()
            .all(|point| mesh.contains_point(point, 0)));
    }

    #[test]
    fn test_sampled_points_of_a_region_are_in_that_region() {
        let mesh = l_shaped_mesh();

        assert!(sample_points(2, 1000)
            .iter()
            .all(|point| mesh.contains_point(point, 2)));
    }

    #[test]
    fn test_share_of_sampled_points_in_a_region_is_its_volume_share() {
        // Region 2 has 1/3 of the volume of the L shape
        let points = sample_points(0, 30000);
        let points_in_region_2 = points.iter().filter(|point| point.y > 0.).count();

        assert!(
            (9700..=10300).contains(&points_in_region_2),
            "{points_in_region_2} of 30000 points are in region 2"
        );
    }

    #[test]
    fn test_sampled_points_are_uniform_in_a_tetrahedron() {
        // The center of mass of a tetrahedron is the mean of its vertices
        let vertices = [
            Point3::origin(),
            Point3::new(4., 0., 0.),
            Point3::new(0., 4., 0.),
            Point3::new(0., 0., 4.),
        ];
        let mesh = DomainMesh::new(&[MeshTetrahedron {
            vertices,
            region: 0,
        }])
        .unwrap();
        let mut rng = Mt64::new(1);
        let count = 30000;
        let sum = (0..count).fold(Vector3::zeros(), |sum, _| {
            sum + mesh.sample_point(0, &mut rng).coords
        });

        // The points have the center of the bounding box (2, 2, 2) at the origin
        let center_of_mass = sum / count as f64 + Vector3::new(2., 2., 2.);
        assert!(
            (center_of_mass - Vector3::new(1., 1., 1.)).norm() < 0.03,
            "the center of mass is {center_of_mass}"
        );
    }

    #[test]
    fn test_polygon_with_its_center_outside_the_mesh_intersects_a_cell() {
        // The center is in the empty corner, and the square goes 0.5 into the cells
        assert!(l_shaped_mesh().intersects_polygon(&square([0.5, 0.5, 0.])));
    }

    #[test]
    fn test_polygon_in_the_empty_corner_intersects_no_cell() {
        assert!(!l_shaped_mesh().intersects_polygon(&square([2.5, 2.5, 0.])));
    }

    #[test]
    fn test_polygon_larger_than_the_mesh_intersects_a_cell() {
        let large_square: Vec<Point3<f64>> = square([0., 0., 0.])
            .iter()
            .map(|vertex| Point3::new(vertex.x * 100., vertex.y * 100., 0.))
            .collect();

        assert!(l_shaped_mesh().intersects_polygon(&large_square));
    }
}
