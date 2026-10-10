use std::path::Path;
use std::str::FromStr;

use parry3d_f64::na::Point3;
use roxmltree::{Document, Node};

use crate::domain_mesh::{DomainMesh, MeshTetrahedron};
use crate::error::DfngenError;

/// Name of the cell data array which gives the region of each cell
pub(super) const REGION_ARRAY: &str = "region";

/// Reads the domain from a VTK XML UnstructuredGrid file (`.vtu`, ASCII)
///
/// The cells can be tetrahedra, voxels, hexahedra, wedges and pyramids. The
/// cell data array `region` (optional) gives the region of each cell.
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

    for piece in root.descendants().filter(|node| node.has_tag_name("Piece")) {
        tetrahedra.extend(read_piece(piece)?);
    }

    DomainMesh::new(&tetrahedra)
}

/// Reads the cells of one `Piece` element and splits them into tetrahedra.
fn read_piece(piece: Node) -> Result<Vec<MeshTetrahedron>, String> {
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

    let regions = match piece
        .children()
        .find(|node| node.has_tag_name("CellData"))
        .and_then(|cell_data| find_data_array(cell_data, REGION_ARRAY).ok())
    {
        Some(region_array) => read_regions(region_array, types.len())?,
        None => vec![0; types.len()],
    };

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

        tetrahedra.extend(local_tetrahedra.iter().map(|local_ids| MeshTetrahedron {
            vertices: local_ids.map(|local_id| points[point_ids[local_id]]),
            region: regions[cell],
        }));

        cell_start = cell_end;
    }

    Ok(tetrahedra)
}

/// Reads the `region` array: one integer of 0 or more for each cell.
/// The array can have a real type, but each value must be a whole number.
pub(super) fn read_regions(region_array: Node, cell_count: usize) -> Result<Vec<u64>, String> {
    let values: Vec<f64> = read_values(region_array)?;

    if values.len() != cell_count {
        return Err(format!(
            "the {} array must have one value for each cell",
            REGION_ARRAY
        ));
    }

    values
        .iter()
        .map(|&value| {
            if value >= 0. && value.fract() == 0. {
                Ok(value as u64)
            } else {
                Err(format!(
                    "the {} value {} is not an integer of 0 or more",
                    REGION_ARRAY, value
                ))
            }
        })
        .collect()
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

    fn region_array(vtk_type: &str, values: &str) -> String {
        format!(
            r#"<DataArray type="{}" Name="region" format="ascii">{}</DataArray>"#,
            vtk_type, values
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

    #[test]
    fn test_region_array_gives_the_regions() {
        let mesh = parse_domain_mesh(&two_cubes_file(&region_array("Int32", "7 3"))).unwrap();
        let values: Vec<u64> = mesh.regions().iter().map(|region| region.value).collect();

        assert_eq!(values, [3, 7]);
    }

    #[test]
    fn test_region_array_with_real_type_is_read() {
        let mesh = parse_domain_mesh(&two_cubes_file(&region_array("Float64", "1.0 2.0"))).unwrap();

        assert_eq!(mesh.regions().len(), 2);
    }

    #[test]
    fn test_region_value_which_is_not_an_integer_is_rejected() {
        let message = error_message(&two_cubes_file(&region_array("Float64", "1.5 2")));

        assert!(message.contains("not an integer of 0 or more"), "{message}");
    }

    #[test]
    fn test_negative_region_value_is_rejected() {
        let message = error_message(&two_cubes_file(&region_array("Int32", "-1 2")));

        assert!(message.contains("not an integer of 0 or more"), "{message}");
    }

    #[test]
    fn test_region_array_with_too_few_values_is_rejected() {
        let message = error_message(&two_cubes_file(&region_array("Int32", "1")));

        assert!(message.contains("one value for each cell"), "{message}");
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
