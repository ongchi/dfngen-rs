use std::path::Path;
use std::str::FromStr;

use parry3d_f64::na::Point3;
use roxmltree::{Document, Node};

use crate::domain_mesh::{DomainMesh, MeshTetrahedron};
use crate::error::DfngenError;

/// Prefix of the name of a cell data array which gives one region:
/// `region__<name>`
const REGION_ARRAY_PREFIX: &str = "region__";

/// One cell data array `region__<name>`
pub(super) struct RegionArray {
    /// `<name>` of the array name
    pub name: String,
    /// For each cell of the piece: true if its value is above 0
    pub is_cell_in_region: Vec<bool>,
}

/// Reads the domain from a VTK XML UnstructuredGrid file (`.vtu`, ASCII)
///
/// The cells can be tetrahedra, voxels, hexahedra, wedges and pyramids. The
/// cell data arrays `region__<name>` (optional) give the regions of each cell.
pub fn read_domain_mesh(path: &Path) -> Result<DomainMesh, DfngenError> {
    let invalid_mesh =
        |message: String| DfngenError::InvalidInput(format!("{}: {}", path.display(), message));

    let text = std::fs::read_to_string(path).map_err(|e| {
        invalid_mesh(format!(
            "cannot read the file: {}. The file must be an ASCII .vtu file.",
            e
        ))
    })?;

    parse_domain_mesh(&text).map_err(invalid_mesh)
}

fn parse_domain_mesh(text: &str) -> Result<DomainMesh, String> {
    let document = Document::parse(text).map_err(|e| e.to_string())?;
    let root = document.root_element();

    if !root.has_tag_name("VTKFile") || root.attribute("type") != Some("UnstructuredGrid") {
        return Err("the file is not a VTK UnstructuredGrid (.vtu) file".to_string());
    }

    let mut tetrahedra = Vec::new();
    let mut region_names = Vec::new();

    for piece in root.descendants().filter(|node| node.has_tag_name("Piece")) {
        tetrahedra.extend(read_piece(piece, &mut region_names)?);
    }

    DomainMesh::new(&tetrahedra, &region_names)
}

/// Reads the cells of one `Piece` element and splits them into tetrahedra.
///
/// # Arguments
///
/// * `piece` - `Piece` element
/// * `region_names` - Names of the regions of the pieces before this one. The new names of this piece are added.
fn read_piece(piece: Node, region_names: &mut Vec<String>) -> Result<Vec<MeshTetrahedron>, String> {
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

    // Each region array with its position in the region names
    let mut region_arrays = Vec::new();

    for region_array in read_region_arrays(piece, types.len())? {
        let position = match region_names
            .iter()
            .position(|name| *name == region_array.name)
        {
            Some(position) => position,
            None => {
                region_names.push(region_array.name.clone());
                region_names.len() - 1
            }
        };
        region_arrays.push((position, region_array));
    }

    let mut tetrahedra = Vec::new();
    let mut cell_start = 0;

    for (cell, (&cell_end, &cell_type)) in offsets.iter().zip(&types).enumerate() {
        let point_ids = connectivity.get(cell_start..cell_end).ok_or(format!(
            "cell {}: the offset {} is not valid",
            cell, cell_end
        ))?;

        if let Some(point_id) = point_ids.iter().find(|&&id| id >= points.len()) {
            return Err(format!("cell {}: there is no point {}", cell, point_id));
        }

        let local_tetrahedra = split_cell_into_tetrahedra(cell_type, point_ids.len())
            .map_err(|message| format!("cell {}: {}", cell, message))?;
        let cell_regions: Vec<usize> = region_arrays
            .iter()
            .filter(|(_, region_array)| region_array.is_cell_in_region[cell])
            .map(|&(position, _)| position)
            .collect();

        tetrahedra.extend(local_tetrahedra.iter().map(|local_ids| MeshTetrahedron {
            vertices: local_ids.map(|local_id| points[point_ids[local_id]]),
            regions: cell_regions.clone(),
        }));

        cell_start = cell_end;
    }

    Ok(tetrahedra)
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
) -> Result<Vec<RegionArray>, String> {
    let Some(cell_data) = piece.children().find(|node| node.has_tag_name("CellData")) else {
        return Ok(Vec::new());
    };

    let mut region_arrays = Vec::new();

    for data_array in cell_data
        .children()
        .filter(|node| node.has_tag_name("DataArray"))
    {
        let array_name = data_array.attribute("Name").unwrap_or_default();
        let Some(name) = array_name.strip_prefix(REGION_ARRAY_PREFIX) else {
            continue;
        };

        if name.is_empty() {
            return Err(format!(
                "the array {} has no region name after the prefix",
                array_name
            ));
        }
        if read_component_count(data_array)? != 1 {
            return Err(format!(
                "the array {} must have one value for each cell (NumberOfComponents=\"1\")",
                array_name
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

        region_arrays.push(RegionArray {
            name: name.to_string(),
            is_cell_in_region: values.iter().map(|&value| value > 0.).collect(),
        });
    }

    Ok(region_arrays)
}

/// Splits a cell into tetrahedra.
///
/// # Arguments
///
/// * `cell_type` - VTK cell type number
/// * `point_count` - Number of points which the file gives for the cell
///
/// # Returns
///
/// The 4 vertices of each tetrahedron, as positions in the point list of the cell
fn split_cell_into_tetrahedra(
    cell_type: u8,
    point_count: usize,
) -> Result<&'static [[usize; 4]], String> {
    let (name, necessary_point_count, tetrahedra): (_, _, &[[usize; 4]]) = match cell_type {
        10 => ("tetrahedron", 4, &[[0, 1, 2, 3]]),
        // The 6 tetrahedra of a voxel or hexahedron have one diagonal in common
        11 => (
            "voxel",
            8,
            &[
                [0, 1, 3, 7],
                [0, 3, 2, 7],
                [0, 2, 6, 7],
                [0, 6, 4, 7],
                [0, 4, 5, 7],
                [0, 5, 1, 7],
            ],
        ),
        12 => (
            "hexahedron",
            8,
            &[
                [0, 1, 2, 6],
                [0, 2, 3, 6],
                [0, 3, 7, 6],
                [0, 7, 4, 6],
                [0, 4, 5, 6],
                [0, 5, 1, 6],
            ],
        ),
        13 => ("wedge", 6, &[[0, 1, 2, 3], [1, 2, 4, 3], [2, 5, 4, 3]]),
        14 => ("pyramid", 5, &[[0, 1, 2, 4], [0, 2, 3, 4]]),
        _ => {
            return Err(format!(
                "the cell type {} is not supported. Use tetrahedra (10), voxels (11), hexahedra (12), wedges (13) or pyramids (14).",
                cell_type
            ));
        }
    };

    if point_count != necessary_point_count {
        return Err(format!(
            "a {} must have {} points, but it has {}",
            name, necessary_point_count, point_count
        ));
    }

    Ok(tetrahedra)
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
            "the DataArray {} has the format \"{}\". Save the file as ASCII.",
            name, format
        ));
    }

    data_array
        .text()
        .unwrap_or_default()
        .split_whitespace()
        .map(|value| {
            value
                .parse()
                .map_err(|_| format!("the DataArray {} has the invalid value {}", name, value))
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
            &format!("{}  0 0 2  1 0 2  1 1 2  0 1 2", UNIT_CUBE_POINTS),
            "0 1 2 3 4 5 6 7  4 5 6 7 8 9 10 11",
            "8 16",
            "12 12",
            cell_data,
        )
    }

    fn region_array(name: &str, values: &str) -> String {
        format!(
            r#"<DataArray type="Int32" Name="region__{}" format="ascii">{}</DataArray>"#,
            name, values
        )
    }

    fn assert_volume(file: &str, expected: f64) {
        let volume = parse_domain_mesh(file).unwrap().volume();

        assert!(
            (volume - expected).abs() < 1e-12,
            "{} is not {}",
            volume,
            expected
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

    #[test]
    fn test_file_of_another_vtk_type_is_rejected() {
        let file =
            one_cell_file(12, "0 1 2 3 4 5 6 7").replace("UnstructuredGrid", "StructuredGrid");

        assert!(error_message(&file).contains("not a VTK UnstructuredGrid"));
    }
}
