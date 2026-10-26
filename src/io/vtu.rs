use std::path::Path;
use std::str::FromStr;

use parry3d_f64::na::Point3;
use roxmltree::{Document, Node};

use crate::domain_mesh::{DomainMesh, MeshSurfaceTriangle, MeshTetrahedron};
use crate::error::DfngenError;

/// Prefix of the name of a cell data array which gives one region:
/// `region__<name>`
const REGION_ARRAY_PREFIX: &str = "region__";

/// Prefix of the name of a cell data array which gives one boundary surface:
/// `boundary__<name>`
const BOUNDARY_ARRAY_PREFIX: &str = "boundary__";

/// Name of the boundary of the mesh surface in the output
/// (`boundary__mesh_surface`): not a valid name of a boundary surface
const MESH_SURFACE_BOUNDARY_NAME: &str = "mesh_surface";

/// One cell data array `<prefix><name>`: `region__<name>` or `boundary__<name>`
pub(super) struct MarkerArray {
    /// `<name>` of the array name
    pub name: String,
    /// For each cell of the piece: true if its value is above 0
    pub is_cell_marked: Vec<bool>,
}

/// The cells of a piece, split into tetrahedra and triangles
#[derive(Default)]
struct MeshCells {
    tetrahedra: Vec<MeshTetrahedron>,
    surface_triangles: Vec<MeshSurfaceTriangle>,
}

/// The parts of one cell
enum CellParts {
    /// A volume cell: the 4 vertices of each tetrahedron
    Volume(&'static [[usize; 4]]),
    /// A surface cell: the 3 vertices of each triangle
    Surface(&'static [[usize; 3]]),
}

/// Reads the domain from a VTK XML UnstructuredGrid file (`.vtu`, ASCII)
///
/// The volume cells can be tetrahedra, voxels, hexahedra, wedges and
/// pyramids. The cell data arrays `region__<name>` (optional) give the
/// regions of each volume cell. The surface cells (triangles and quads,
/// optional) with the cell data arrays `boundary__<name>` give the boundary
/// surfaces.
pub fn read_domain_mesh(path: &Path) -> Result<DomainMesh, DfngenError> {
    let invalid_mesh =
        |message: String| DfngenError::InvalidInput(format!("{}: {}", path.display(), message));

    let text = std::fs::read_to_string(path).map_err(|e| {
        invalid_mesh(format!(
            "cannot read the file: {e}. The file must be an ASCII .vtu file."
        ))
    })?;

    parse_domain_mesh(&text).map_err(invalid_mesh)
}

pub(super) fn parse_domain_mesh(text: &str) -> Result<DomainMesh, String> {
    let document = Document::parse(text).map_err(|e| e.to_string())?;
    let root = document.root_element();

    if !root.has_tag_name("VTKFile") || root.attribute("type") != Some("UnstructuredGrid") {
        return Err("the file is not a VTK UnstructuredGrid (.vtu) file".to_string());
    }

    let mut cells = MeshCells::default();
    let mut names = MarkerNames::default();

    for piece in root.descendants().filter(|node| node.has_tag_name("Piece")) {
        read_piece(piece, &mut names, &mut cells)?;
    }

    check_boundary_names(&names)?;

    Ok(DomainMesh::new(&cells.tetrahedra, &names.regions)?
        .with_boundary_surfaces(&cells.surface_triangles, &names.boundaries))
}

/// Names of the regions and of the boundary surfaces of the pieces, in the
/// order of the file
#[derive(Default)]
struct MarkerNames {
    regions: Vec<String>,
    boundaries: Vec<String>,
}

/// A boundary surface must not have the name of a region or of the mesh
/// surface: each of them is an array `boundary__<name>` of the output.
fn check_boundary_names(names: &MarkerNames) -> Result<(), String> {
    for name in &names.boundaries {
        if name == MESH_SURFACE_BOUNDARY_NAME {
            return Err(format!(
                "the array {BOUNDARY_ARRAY_PREFIX}{name} has the name of the mesh surface boundary of the output. Use another name."
            ));
        }
        if names.regions.contains(name) {
            return Err(format!(
                "the boundary surface {name} (array {BOUNDARY_ARRAY_PREFIX}{name}) has the name of the region {name} (array {REGION_ARRAY_PREFIX}{name}): both give the output array {BOUNDARY_ARRAY_PREFIX}{name}. Rename one of them."
            ));
        }
    }

    Ok(())
}

/// Position of each array in the names, with the new names added
fn positions_in_names(
    arrays: Vec<MarkerArray>,
    names: &mut Vec<String>,
) -> Vec<(usize, MarkerArray)> {
    arrays
        .into_iter()
        .map(|array| {
            let position = match names.iter().position(|name| *name == array.name) {
                Some(position) => position,
                None => {
                    names.push(array.name.clone());
                    names.len() - 1
                }
            };
            (position, array)
        })
        .collect()
}

/// The positions in the names of the arrays which mark the cell
fn marked_positions(arrays: &[(usize, MarkerArray)], cell: usize) -> Vec<usize> {
    arrays
        .iter()
        .filter(|(_, array)| array.is_cell_marked[cell])
        .map(|&(position, _)| position)
        .collect()
}

/// Reads the cells of one `Piece` element and splits them into tetrahedra
/// (volume cells) and triangles (surface cells).
///
/// # Arguments
///
/// * `piece` - `Piece` element
/// * `names` - Names of the regions and boundaries of the pieces before this one. The new names of this piece are added.
/// * `cells` - The cells of the pieces before this one. The cells of this piece are added.
fn read_piece(piece: Node, names: &mut MarkerNames, cells: &mut MeshCells) -> Result<(), String> {
    let points_element = find_child(piece, "Points")?;
    let cells_element = find_child(piece, "Cells")?;

    let coordinates: Vec<f64> = read_values(find_child(points_element, "DataArray")?)?;
    let connectivity: Vec<usize> = read_values(find_data_array(cells_element, "connectivity")?)?;
    let offsets: Vec<usize> = read_values(find_data_array(cells_element, "offsets")?)?;
    let types: Vec<u8> = read_values(find_data_array(cells_element, "types")?)?;

    let (point_coordinates, remaining_coordinates) = coordinates.as_chunks::<3>();

    if !remaining_coordinates.is_empty() {
        return Err("the points must have 3 coordinates each".to_string());
    }
    if offsets.len() != types.len() {
        return Err("offsets and types must have one value for each cell".to_string());
    }

    let points: Vec<Point3<f64>> = point_coordinates
        .iter()
        .map(|&coordinates| Point3::from(coordinates))
        .collect();

    // Each array with its position in the names
    let region_arrays =
        positions_in_names(read_region_arrays(piece, types.len())?, &mut names.regions);
    let boundary_arrays = positions_in_names(
        read_marker_arrays(piece, types.len(), BOUNDARY_ARRAY_PREFIX, "boundary")?,
        &mut names.boundaries,
    );

    let mut cell_start = 0;

    for (cell, (&cell_end, &cell_type)) in offsets.iter().zip(&types).enumerate() {
        let point_ids = connectivity
            .get(cell_start..cell_end)
            .ok_or(format!("cell {cell}: the offset {cell_end} is not valid"))?;

        if let Some(point_id) = point_ids.iter().find(|&&id| id >= points.len()) {
            return Err(format!("cell {cell}: there is no point {point_id}"));
        }

        let parts = split_cell(cell_type, point_ids.len())
            .map_err(|message| format!("cell {cell}: {message}"))?;
        let cell_regions = marked_positions(&region_arrays, cell);
        let cell_boundaries = marked_positions(&boundary_arrays, cell);

        match parts {
            CellParts::Volume(local_tetrahedra) => {
                if let Some(&position) = cell_boundaries.first() {
                    return Err(format!(
                        "cell {cell}: the array {BOUNDARY_ARRAY_PREFIX}{} has a value above 0 for a volume cell. Only surface cells (triangles, quads) can be on a boundary surface.",
                        names.boundaries[position]
                    ));
                }
                cells
                    .tetrahedra
                    .extend(local_tetrahedra.iter().map(|local_ids| MeshTetrahedron {
                        vertices: local_ids.map(|local_id| points[point_ids[local_id]]),
                        regions: cell_regions.clone(),
                    }));
            }
            CellParts::Surface(local_triangles) => {
                if let Some(&position) = cell_regions.first() {
                    return Err(format!(
                        "cell {cell}: the array {REGION_ARRAY_PREFIX}{} has a value above 0 for a surface cell. Only volume cells can be in a region.",
                        names.regions[position]
                    ));
                }
                cells
                    .surface_triangles
                    .extend(local_triangles.iter().map(|local_ids| MeshSurfaceTriangle {
                        vertices: local_ids.map(|local_id| points[point_ids[local_id]]),
                        boundaries: cell_boundaries.clone(),
                    }));
            }
        }

        cell_start = cell_end;
    }

    Ok(())
}

/// Reads the `NumberOfComponents` attribute of a `DataArray` element: the
/// number of values for each cell. Without the attribute it is 1.
pub(super) fn read_component_count(data_array: Node) -> Result<usize, String> {
    let Some(text) = data_array.attribute("NumberOfComponents") else {
        return Ok(1);
    };

    match text.trim().parse() {
        Ok(component_count) if component_count > 0 => Ok(component_count),
        _ => Err(format!(
            "the DataArray {} has the invalid NumberOfComponents \"{}\"",
            data_array.attribute("Name").unwrap_or_default(),
            text
        )),
    }
}

/// Reads the cell data arrays `region__<name>` of a piece, in the order of
/// the file. Each array is one region, with one value for each cell: a cell
/// with a value above 0 is in the region.
///
/// # Arguments
///
/// * `piece` - `Piece` element
/// * `cell_count` - Number of cells of the piece
pub(super) fn read_region_arrays(
    piece: Node,
    cell_count: usize,
) -> Result<Vec<MarkerArray>, String> {
    read_marker_arrays(piece, cell_count, REGION_ARRAY_PREFIX, "region")
}

/// Reads the cell data arrays `<prefix><name>` of a piece, in the order of
/// the file, with one value for each cell: a cell with a value above 0 is
/// marked.
///
/// # Arguments
///
/// * `piece` - `Piece` element
/// * `cell_count` - Number of cells of the piece
/// * `prefix` - Prefix of the array names
/// * `kind` - What an array gives, for the error messages: "region" or "boundary"
fn read_marker_arrays(
    piece: Node,
    cell_count: usize,
    prefix: &str,
    kind: &str,
) -> Result<Vec<MarkerArray>, String> {
    let Some(cell_data) = piece.children().find(|node| node.has_tag_name("CellData")) else {
        return Ok(Vec::new());
    };

    let mut arrays = Vec::new();

    for data_array in cell_data
        .children()
        .filter(|node| node.has_tag_name("DataArray"))
    {
        let array_name = data_array.attribute("Name").unwrap_or_default();
        let Some(name) = array_name.strip_prefix(prefix) else {
            continue;
        };

        if name.is_empty() {
            return Err(format!(
                "the array {array_name} has no {kind} name after the prefix"
            ));
        }
        if read_component_count(data_array)? != 1 {
            return Err(format!(
                "the array {array_name} must have one value for each cell (NumberOfComponents=\"1\")"
            ));
        }

        let values: Vec<f64> = read_values(data_array)?;

        if values.len() != cell_count {
            return Err(format!(
                "the array {} must have one value for each cell, but it has {} values for {} cells",
                array_name,
                values.len(),
                cell_count
            ));
        }

        arrays.push(MarkerArray {
            name: name.to_string(),
            is_cell_marked: values.iter().map(|&value| value > 0.).collect(),
        });
    }

    Ok(arrays)
}

/// Splits a volume cell into tetrahedra, or a surface cell into triangles.
///
/// # Arguments
///
/// * `cell_type` - VTK cell type number
/// * `point_count` - Number of points which the file gives for the cell
///
/// # Returns
///
/// The vertices of each tetrahedron or triangle, as positions in the point list of the cell
fn split_cell(cell_type: u8, point_count: usize) -> Result<CellParts, String> {
    let (name, necessary_point_count, parts) = match cell_type {
        5 => ("triangle", 3, CellParts::Surface(&[[0, 1, 2]])),
        9 => ("quad", 4, CellParts::Surface(&[[0, 1, 2], [0, 2, 3]])),
        10 => ("tetrahedron", 4, CellParts::Volume(&[[0, 1, 2, 3]])),
        // The 6 tetrahedra of a voxel or hexahedron have one diagonal in common
        11 => (
            "voxel",
            8,
            CellParts::Volume(&[
                [0, 1, 3, 7],
                [0, 3, 2, 7],
                [0, 2, 6, 7],
                [0, 6, 4, 7],
                [0, 4, 5, 7],
                [0, 5, 1, 7],
            ]),
        ),
        12 => (
            "hexahedron",
            8,
            CellParts::Volume(&[
                [0, 1, 2, 6],
                [0, 2, 3, 6],
                [0, 3, 7, 6],
                [0, 7, 4, 6],
                [0, 4, 5, 6],
                [0, 5, 1, 6],
            ]),
        ),
        13 => (
            "wedge",
            6,
            CellParts::Volume(&[[0, 1, 2, 3], [1, 2, 4, 3], [2, 5, 4, 3]]),
        ),
        14 => (
            "pyramid",
            5,
            CellParts::Volume(&[[0, 1, 2, 4], [0, 2, 3, 4]]),
        ),
        _ => {
            return Err(format!(
                "the cell type {cell_type} is not supported. Use tetrahedra (10), voxels (11), hexahedra (12), wedges (13) or pyramids (14), and for boundary surfaces triangles (5) or quads (9)."
            ));
        }
    };

    if point_count != necessary_point_count {
        return Err(format!(
            "a {name} must have {necessary_point_count} points, but it has {point_count}"
        ));
    }

    Ok(parts)
}

pub(super) fn find_child<'a, 'input>(
    parent: Node<'a, 'input>,
    tag_name: &str,
) -> Result<Node<'a, 'input>, String> {
    parent
        .children()
        .find(|node| node.has_tag_name(tag_name))
        .ok_or(format!(
            "{} has no {} element",
            parent.tag_name().name(),
            tag_name
        ))
}

pub(super) fn find_data_array<'a, 'input>(
    parent: Node<'a, 'input>,
    name: &str,
) -> Result<Node<'a, 'input>, String> {
    parent
        .children()
        .find(|node| node.has_tag_name("DataArray") && node.attribute("Name") == Some(name))
        .ok_or(format!(
            "{} has no DataArray with the name {}",
            parent.tag_name().name(),
            name
        ))
}

/// Reads the numbers of an ASCII `DataArray` element.
pub(super) fn read_values<T: FromStr>(data_array: Node) -> Result<Vec<T>, String> {
    let name = data_array.attribute("Name").unwrap_or("Points");
    let format = data_array.attribute("format").unwrap_or("ascii");

    if format != "ascii" {
        return Err(format!(
            "the DataArray {name} has the format \"{format}\". Save the file as ASCII."
        ));
    }

    data_array
        .text()
        .unwrap_or_default()
        .split_whitespace()
        .map(|value| {
            value
                .parse()
                .map_err(|_| format!("the DataArray {name} has the invalid value {value}"))
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    const UNIT_CUBE_POINTS: &str = "0 0 0  1 0 0  1 1 0  0 1 0  0 0 1  1 0 1  1 1 1  0 1 1";

    /// A `.vtu` file with one piece.
    /// `cell_data` is the content of the `CellData` element.
    fn vtu_file(
        points: &str,
        connectivity: &str,
        offsets: &str,
        types: &str,
        cell_data: &str,
    ) -> String {
        format!(
            r#"<?xml version="1.0"?>
            <VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian">
            <UnstructuredGrid>
            <Piece NumberOfPoints="0" NumberOfCells="0">
            <CellData>{cell_data}</CellData>
            <Points>
            <DataArray type="Float64" Name="Points" NumberOfComponents="3" format="ascii">{points}</DataArray>
            </Points>
            <Cells>
            <DataArray type="Int64" Name="connectivity" format="ascii">{connectivity}</DataArray>
            <DataArray type="Int64" Name="offsets" format="ascii">{offsets}</DataArray>
            <DataArray type="UInt8" Name="types" format="ascii">{types}</DataArray>
            </Cells>
            </Piece>
            </UnstructuredGrid>
            </VTKFile>"#
        )
    }

    /// A file with one cell which uses the first points of the unit cube.
    fn one_cell_file(cell_type: u8, connectivity: &str) -> String {
        let point_count = connectivity.split_whitespace().count();

        vtu_file(
            UNIT_CUBE_POINTS,
            connectivity,
            &point_count.to_string(),
            &cell_type.to_string(),
            "",
        )
    }

    /// Two unit cubes, one on the other, with the given `CellData` content.
    fn two_cubes_file(cell_data: &str) -> String {
        vtu_file(
            &format!("{UNIT_CUBE_POINTS}  0 0 2  1 0 2  1 1 2  0 1 2"),
            "0 1 2 3 4 5 6 7  4 5 6 7 8 9 10 11",
            "8 16",
            "12 12",
            cell_data,
        )
    }

    fn region_array(name: &str, values: &str) -> String {
        format!(
            r#"<DataArray type="Int32" Name="region__{name}" format="ascii">{values}</DataArray>"#
        )
    }

    fn assert_volume(file: &str, expected: f64) {
        let volume = parse_domain_mesh(file).unwrap().volume();

        assert!(
            (volume - expected).abs() < 1e-12,
            "{volume} is not {expected}"
        );
    }

    fn error_message(file: &str) -> String {
        parse_domain_mesh(file).expect_err("the file must be rejected")
    }

    #[test]
    fn test_tetrahedron_has_its_volume() {
        assert_volume(&one_cell_file(10, "0 1 3 4"), 1. / 6.);
    }

    #[test]
    fn test_voxel_has_its_volume() {
        // The points of a voxel are in the order of x, then y, then z
        assert_volume(&one_cell_file(11, "0 1 3 2 4 5 7 6"), 1.);
    }

    #[test]
    fn test_hexahedron_has_its_volume() {
        assert_volume(&one_cell_file(12, "0 1 2 3 4 5 6 7"), 1.);
    }

    #[test]
    fn test_wedge_has_its_volume() {
        assert_volume(&one_cell_file(13, "0 1 3 4 5 7"), 0.5);
    }

    #[test]
    fn test_pyramid_has_its_volume() {
        assert_volume(&one_cell_file(14, "0 1 2 3 4"), 1. / 3.);
    }

    #[test]
    fn test_hexahedron_contains_its_center() {
        let mesh = parse_domain_mesh(&one_cell_file(12, "0 1 2 3 4 5 6 7")).unwrap();

        assert!(mesh.contains_point(&Point3::origin(), 0));
    }

    #[test]
    fn test_mesh_without_region_array_has_no_region() {
        let mesh = parse_domain_mesh(&two_cubes_file("")).unwrap();

        assert!(mesh.regions().is_empty());
    }

    /// Name and volume (rounded to 9 decimals) of each region of the two cubes
    fn region_volumes(cell_data: &str) -> Vec<(String, f64)> {
        parse_domain_mesh(&two_cubes_file(cell_data))
            .unwrap()
            .regions()
            .iter()
            .map(|region| (region.name.clone(), (region.volume * 1e9).round() / 1e9))
            .collect()
    }

    fn named_volumes(volumes: &[(&str, f64)]) -> Vec<(String, f64)> {
        volumes
            .iter()
            .map(|&(name, volume)| (name.to_string(), volume))
            .collect()
    }

    #[test]
    fn test_each_region_array_is_one_region_in_the_order_of_the_file() {
        let cell_data = region_array("rock", "1 0") + &region_array("clay", "0 1");

        assert_eq!(
            region_volumes(&cell_data),
            named_volumes(&[("rock", 1.), ("clay", 1.)])
        );
    }

    #[test]
    fn test_cell_with_a_value_above_0_in_two_arrays_is_in_two_regions() {
        let cell_data = region_array("rock", "1 1") + &region_array("clay", "0 5");

        assert_eq!(
            region_volumes(&cell_data),
            named_volumes(&[("rock", 2.), ("clay", 1.)])
        );
    }

    #[test]
    fn test_region_array_with_real_type_is_read() {
        let cell_data = region_array("rock", "1.0 0.0").replace("Int32", "Float64");

        assert_eq!(region_volumes(&cell_data), named_volumes(&[("rock", 1.)]));
    }

    #[test]
    fn test_array_with_the_name_region_is_ignored() {
        let cell_data = region_array("rock", "1 2").replace("region__rock", "region");
        let mesh = parse_domain_mesh(&two_cubes_file(&cell_data)).unwrap();

        assert!(mesh.regions().is_empty());
    }

    #[test]
    fn test_region_array_with_too_few_values_is_rejected() {
        let message = error_message(&two_cubes_file(&region_array("rock", "1")));

        assert!(
            message.contains("region__rock must have one value for each cell"),
            "{message}"
        );
    }

    #[test]
    fn test_region_array_with_two_components_is_rejected() {
        let cell_data =
            region_array("rock", "1 2 3 4").replace("format=", r#"NumberOfComponents="2" format="#);
        let message = error_message(&two_cubes_file(&cell_data));

        assert!(
            message.contains("region__rock must have one value for each cell"),
            "{message}"
        );
    }

    #[test]
    fn test_region_array_without_a_region_name_is_rejected() {
        let cell_data = region_array("", "1 0");
        let message = error_message(&two_cubes_file(&cell_data));

        assert!(message.contains("has no region name"), "{message}");
    }

    #[test]
    fn test_data_array_with_0_components_is_rejected() {
        let cell_data =
            region_array("rock", "1 2").replace("format=", r#"NumberOfComponents="0" format="#);
        let message = error_message(&two_cubes_file(&cell_data));

        assert!(message.contains("invalid NumberOfComponents"), "{message}");
    }

    #[test]
    fn test_binary_data_array_is_rejected() {
        let file = one_cell_file(12, "0 1 2 3 4 5 6 7").replace(
            r#"Name="offsets" format="ascii""#,
            r#"Name="offsets" format="binary""#,
        );

        assert!(error_message(&file).contains("Save the file as ASCII"));
    }

    #[test]
    fn test_unknown_cell_type_is_rejected() {
        // 24 is a quadratic tetrahedron
        let message = error_message(&one_cell_file(24, "0 1 2 3"));

        assert!(
            message.contains("the cell type 24 is not supported"),
            "{message}"
        );
    }

    #[test]
    fn test_cell_with_too_few_points_is_rejected() {
        let message = error_message(&one_cell_file(12, "0 1 2 3"));

        assert!(
            message.contains("a hexahedron must have 8 points"),
            "{message}"
        );
    }

    #[test]
    fn test_connectivity_with_unknown_point_is_rejected() {
        let message = error_message(&one_cell_file(10, "0 1 3 8"));

        assert!(message.contains("there is no point 8"), "{message}");
    }

    #[test]
    fn test_offset_beyond_the_connectivity_is_rejected() {
        let file = vtu_file(UNIT_CUBE_POINTS, "0 1 3 4", "5", "10", "");

        assert!(error_message(&file).contains("the offset 5 is not valid"));
    }

    /// A unit cube (cell 0) with its face at z = 0 as a quad (cell 1) and its
    /// face at z = 1 as two triangles (cells 2 and 3), with the given `CellData`
    /// content
    fn cube_with_surface_cells_file(cell_data: &str) -> String {
        vtu_file(
            UNIT_CUBE_POINTS,
            "0 1 2 3 4 5 6 7  0 1 2 3  4 5 6  4 6 7",
            "8 12 15 18",
            "12 9 5 5",
            cell_data,
        )
    }

    fn boundary_array(name: &str, values: &str) -> String {
        format!(
            r#"<DataArray type="Int32" Name="boundary__{name}" format="ascii">{values}</DataArray>"#
        )
    }

    #[test]
    fn test_surface_cells_do_not_change_the_volume() {
        assert_volume(&cube_with_surface_cells_file(""), 1.);
    }

    #[test]
    fn test_each_boundary_array_is_one_boundary_surface_in_the_order_of_the_file() {
        let cell_data = boundary_array("top", "0 0 1 1") + &boundary_array("bottom", "0 1 0 0");
        let mesh = parse_domain_mesh(&cube_with_surface_cells_file(&cell_data)).unwrap();
        let names: Vec<&str> = mesh
            .boundary_surfaces()
            .iter()
            .map(|surface| surface.name.as_str())
            .collect();

        assert_eq!(names, ["top", "bottom"]);
    }

    #[test]
    fn test_point_on_a_marked_quad_has_its_boundary_surface() {
        let cell_data = boundary_array("top", "0 0 1 1") + &boundary_array("bottom", "0 1 0 0");
        let mesh = parse_domain_mesh(&cube_with_surface_cells_file(&cell_data)).unwrap();

        // The center of the cube is at the origin
        assert_eq!(
            mesh.boundary_surfaces_at(&Point3::new(0.2, 0.3, -0.5), 1e-9),
            [2]
        );
    }

    #[test]
    fn test_boundary_value_above_0_on_a_volume_cell_is_rejected() {
        let file = cube_with_surface_cells_file(&boundary_array("top", "1 0 1 1"));
        let message = error_message(&file);

        assert!(
            message.contains("boundary__top has a value above 0 for a volume cell"),
            "{message}"
        );
    }

    #[test]
    fn test_region_value_above_0_on_a_surface_cell_is_rejected() {
        let file = cube_with_surface_cells_file(&region_array("rock", "1 1 0 0"));
        let message = error_message(&file);

        assert!(
            message.contains("region__rock has a value above 0 for a surface cell"),
            "{message}"
        );
    }

    #[test]
    fn test_boundary_surface_with_the_name_of_a_region_is_rejected() {
        let cell_data = region_array("west", "1 0 0 0") + &boundary_array("west", "0 1 0 0");
        let message = error_message(&cube_with_surface_cells_file(&cell_data));

        assert!(message.contains("Rename one of them"), "{message}");
    }

    #[test]
    fn test_boundary_surface_with_the_name_mesh_surface_is_rejected() {
        let file = cube_with_surface_cells_file(&boundary_array("mesh_surface", "0 1 0 0"));
        let message = error_message(&file);

        assert!(
            message.contains("name of the mesh surface boundary"),
            "{message}"
        );
    }

    #[test]
    fn test_boundary_array_without_a_boundary_name_is_rejected() {
        let file = cube_with_surface_cells_file(&boundary_array("", "0 1 0 0"));
        let message = error_message(&file);

        assert!(message.contains("has no boundary name"), "{message}");
    }

    #[test]
    fn test_quad_with_too_few_points_is_rejected() {
        let message = error_message(&one_cell_file(9, "0 1 2"));

        assert!(message.contains("a quad must have 4 points"), "{message}");
    }

    #[test]
    fn test_polygon_cell_is_rejected() {
        let message = error_message(&one_cell_file(7, "0 1 2 3 4"));

        assert!(
            message.contains("the cell type 7 is not supported"),
            "{message}"
        );
    }

    #[test]
    fn test_mesh_of_surface_cells_only_is_rejected() {
        let message = error_message(&one_cell_file(9, "0 1 2 3"));

        assert!(message.contains("no cell with a volume"), "{message}");
    }

    #[test]
    fn test_file_of_another_vtk_type_is_rejected() {
        let file =
            one_cell_file(12, "0 1 2 3 4 5 6 7").replace("UnstructuredGrid", "StructuredGrid");

        assert!(error_message(&file).contains("not a VTK UnstructuredGrid"));
    }
}
