use std::path::Path;

use roxmltree::{Document, Node};
use tracing::warn;

use crate::error::DfngenError;
use crate::io::vtu::{
    find_child, find_data_array, read_component_count, read_region_arrays, read_values,
};

/// Cell kinds of a PolyData file which are not fractures. The cell data of
/// a file is in the order of these kinds, with the polygons before the strips.
const CELL_ELEMENTS_BEFORE_POLYGONS: [&str; 2] = ["Verts", "Lines"];
const CELL_ELEMENTS_AFTER_POLYGONS: [&str; 1] = ["Strips"];

/// Name of the cell data array which gives the fracture of each polygon
const FRACTURE_ID_ARRAY: &str = "fracture_id";

/// A polygon cell of the user-defined fracture file
#[derive(Debug)]
pub struct FracturePolygon {
    /// Vertices in the coordinates of the file
    pub vertices: Vec<[f64; 3]>,

    /// Names of the regions of the polygon: `<name>` of each cell data array
    /// `region__<name>` in which the polygon has a value above 0.
    pub regions: Vec<String>,

    /// Value of the cell data array `fracture_id`, `None` if the file has no
    /// such array. Polygons with the same value are pieces of one fracture.
    pub fracture_id: Option<u64>,
}

/// Reads the user-defined fractures from a VTK XML PolyData file (`.vtp`, ASCII)
///
/// Each polygon cell is one fracture, with the points of the cell as its
/// vertices. The fractures are in the order of the cells. The cell data array
/// arrays `region__<name>` (optional) give the regions of each fracture. The cell data array
/// `fracture_id` (optional) gives the fracture of each polygon, for a fracture
/// of more than one polygon. The point data and the other cell data of the
/// file are ignored.
pub fn read_user_fractures(path: &Path) -> Result<Vec<FracturePolygon>, DfngenError> {
    let invalid_file =
        |message: String| DfngenError::InvalidInput(format!("{}: {}", path.display(), message));

    let text = std::fs::read_to_string(path).map_err(|e| {
        invalid_file(format!(
            "cannot read the file: {}. The file must be an ASCII .vtp file.",
            e
        ))
    })?;

    parse_user_fractures(&text).map_err(invalid_file)
}

fn parse_user_fractures(text: &str) -> Result<Vec<FracturePolygon>, String> {
    let document = Document::parse(text).map_err(|e| e.to_string())?;
    let root = document.root_element();

    if !root.has_tag_name("VTKFile") || root.attribute("type") != Some("PolyData") {
        return Err("the file is not a VTK PolyData (.vtp) file".to_string());
    }

    let mut fractures = Vec::new();

    for piece in root.descendants().filter(|node| node.has_tag_name("Piece")) {
        fractures.extend(read_piece(piece, fractures.len())?);
    }

    Ok(fractures)
}

/// Reads the polygons of one `Piece` element.
///
/// # Arguments
///
/// * `piece` - `Piece` element
/// * `fracture_count` - Number of fractures in the pieces before this one, for the error messages
fn read_piece(piece: Node, fracture_count: usize) -> Result<Vec<FracturePolygon>, String> {
    let points_element = find_child(piece, "Points")?;
    let polys_element = find_child(piece, "Polys")?;

    let coordinates: Vec<f64> = read_values(find_child(points_element, "DataArray")?)?;
    let connectivity: Vec<usize> = read_values(find_data_array(polys_element, "connectivity")?)?;
    let offsets: Vec<usize> = read_values(find_data_array(polys_element, "offsets")?)?;

    let (points, remaining_coordinates) = coordinates.as_chunks::<3>();

    if !remaining_coordinates.is_empty() {
        return Err("the points must have 3 coordinates each".to_string());
    }

    let cell_count_before_polygons = count_ignored_cells(piece, &CELL_ELEMENTS_BEFORE_POLYGONS)?;
    let cell_count_after_polygons = count_ignored_cells(piece, &CELL_ELEMENTS_AFTER_POLYGONS)?;

    let cell_count = cell_count_before_polygons + offsets.len() + cell_count_after_polygons;
    let cell_data = piece.children().find(|node| node.has_tag_name("CellData"));

    let region_arrays = read_region_arrays(piece, cell_count)?;

    let user_fracture_ids = cell_data
        .and_then(|cell_data| find_data_array(cell_data, FRACTURE_ID_ARRAY).ok())
        .map(|id_array| read_user_fracture_ids(id_array, cell_count))
        .transpose()?;

    let mut fractures = Vec::with_capacity(offsets.len());
    let mut cell_start = 0;

    for (cell, &cell_end) in offsets.iter().enumerate() {
        let fracture_id = fracture_count + cell + 1;

        let point_ids = connectivity.get(cell_start..cell_end).ok_or(format!(
            "fracture {}: the offset {} is not valid",
            fracture_id, cell_end
        ))?;

        if point_ids.len() < 3 {
            return Err(format!(
                "fracture {}: a fracture must have 3 vertices or more",
                fracture_id
            ));
        }
        if let Some(point_id) = point_ids.iter().find(|&&id| id >= points.len()) {
            return Err(format!(
                "fracture {}: there is no point {}",
                fracture_id, point_id
            ));
        }

        fractures.push(FracturePolygon {
            vertices: point_ids.iter().map(|&id| points[id]).collect(),
            regions: region_arrays
                .iter()
                .filter(|region_array| {
                    region_array.is_cell_in_region[cell_count_before_polygons + cell]
                })
                .map(|region_array| region_array.name.clone())
                .collect(),
            fracture_id: user_fracture_ids
                .as_ref()
                .map(|ids| ids[cell_count_before_polygons + cell]),
        });

        cell_start = cell_end;
    }

    Ok(fractures)
}

/// Reads the `fracture_id` array: one integer of 0 or more for each cell.
fn read_user_fracture_ids(id_array: Node, cell_count: usize) -> Result<Vec<u64>, String> {
    if read_component_count(id_array)? != 1 {
        return Err(format!(
            "the {} array must have one value for each fracture (NumberOfComponents=\"1\")",
            FRACTURE_ID_ARRAY
        ));
    }

    let ids: Vec<u64> = read_values(id_array)?;

    if ids.len() != cell_count {
        return Err(format!(
            "the {} array must have one value for each cell, but it has {} values for {} cells",
            FRACTURE_ID_ARRAY,
            ids.len(),
            cell_count
        ));
    }

    Ok(ids)
}

/// Number of cells of the given kinds in the piece. These cells are not
/// fractures, and each kind which has cells gives a warning.
fn count_ignored_cells(piece: Node, tag_names: &[&str]) -> Result<usize, String> {
    let mut cell_count = 0;

    for &tag_name in tag_names {
        let Some(offsets_array) = piece
            .children()
            .find(|node| node.has_tag_name(tag_name))
            .and_then(|element| find_data_array(element, "offsets").ok())
        else {
            continue;
        };
        let offsets: Vec<usize> = read_values(offsets_array)?;

        if !offsets.is_empty() {
            warn!(
                "The {} cells of the user-defined fracture file are ignored. Only polygons are fractures.",
                tag_name
            );
        }
        cell_count += offsets.len();
    }

    Ok(cell_count)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A square and a triangle which have two points in common.
    const SQUARE_POINTS: &str = "0 0 0  1 0 0  1 1 0  0 1 0  2 0 5";

    fn vtp_file(
        file_type: &str,
        connectivity: &str,
        offsets: &str,
        other_elements: &str,
    ) -> String {
        format!(
            r#"<?xml version="1.0"?>
            <VTKFile type="{file_type}" version="1.0" byte_order="LittleEndian">
            <PolyData>
            <Piece NumberOfPoints="5" NumberOfPolys="2">
            <Points>
            <DataArray type="Float64" Name="points" NumberOfComponents="3" format="ascii">
            {SQUARE_POINTS}
            </DataArray>
            </Points>
            <Polys>
            <DataArray type="Int64" Name="connectivity" format="ascii">{connectivity}</DataArray>
            <DataArray type="Int64" Name="offsets" format="ascii">{offsets}</DataArray>
            </Polys>
            {other_elements}
            </Piece>
            </PolyData>
            </VTKFile>"#
        )
    }

    fn square_and_triangle_file(other_elements: &str) -> String {
        vtp_file("PolyData", "0 1 2 3  1 4 2", "4 7", other_elements)
    }

    fn error_message(file: &str) -> String {
        parse_user_fractures(file).expect_err("the file must be rejected")
    }

    #[test]
    fn test_each_polygon_is_one_fracture_in_cell_order() {
        let fractures = parse_user_fractures(&square_and_triangle_file("")).unwrap();

        let vertex_counts: Vec<usize> = fractures
            .iter()
            .map(|fracture| fracture.vertices.len())
            .collect();
        assert_eq!(vertex_counts, [4, 3]);
    }

    #[test]
    fn test_vertices_are_the_points_of_the_cell_in_its_order() {
        let fractures = parse_user_fractures(&square_and_triangle_file("")).unwrap();

        assert_eq!(
            fractures[1].vertices,
            [[1., 0., 0.], [2., 0., 5.], [1., 1., 0.]]
        );
    }

    #[test]
    fn test_cell_data_is_ignored() {
        let cell_data = r#"<CellData>
            <DataArray type="Int64" Name="family" format="ascii">0 0</DataArray>
            </CellData>"#;

        let fractures = parse_user_fractures(&square_and_triangle_file(cell_data)).unwrap();

        assert_eq!(fractures.len(), 2);
    }

    #[test]
    fn test_file_without_polygons_gives_no_fractures() {
        let fractures = parse_user_fractures(&vtp_file("PolyData", "", "", "")).unwrap();

        assert!(fractures.is_empty());
    }

    #[test]
    fn test_polygon_with_two_points_is_rejected() {
        let message = error_message(&vtp_file("PolyData", "0 1 2 3  1 4", "4 6", ""));

        assert!(
            message.contains("fracture 2: a fracture must have 3 vertices or more"),
            "{message}"
        );
    }

    #[test]
    fn test_point_which_is_not_in_the_file_is_rejected() {
        let message = error_message(&vtp_file("PolyData", "0 1 5", "3", ""));

        assert!(
            message.contains("fracture 1: there is no point 5"),
            "{message}"
        );
    }

    #[test]
    fn test_offset_beyond_the_connectivity_is_rejected() {
        let message = error_message(&vtp_file("PolyData", "0 1 2", "4", ""));

        assert!(message.contains("the offset 4 is not valid"), "{message}");
    }

    #[test]
    fn test_unstructured_grid_file_is_rejected() {
        let message = error_message(&vtp_file("UnstructuredGrid", "0 1 2", "3", ""));

        assert!(message.contains("not a VTK PolyData"), "{message}");
    }

    /// `CellData` element with a `region__<name>` array for each name and its values
    fn region_cell_data(arrays: &[(&str, &str)]) -> String {
        let data_arrays: String = arrays
            .iter()
            .map(|(name, values)| {
                format!(
                    r#"<DataArray type="Int64" Name="region__{name}" format="ascii">{values}</DataArray>"#
                )
            })
            .collect();

        format!("<CellData>{data_arrays}</CellData>")
    }

    fn regions(file: &str) -> Vec<Vec<String>> {
        parse_user_fractures(file)
            .unwrap()
            .into_iter()
            .map(|fracture| fracture.regions)
            .collect()
    }

    #[test]
    fn test_region_array_with_a_value_above_0_gives_a_region_of_the_fracture() {
        let file = square_and_triangle_file(&region_cell_data(&[("rock", "1 0"), ("clay", "0 1")]));

        assert_eq!(regions(&file), [["rock"], ["clay"]]);
    }

    #[test]
    fn test_fracture_with_a_value_above_0_in_two_region_arrays_has_two_regions() {
        let file = square_and_triangle_file(&region_cell_data(&[("rock", "1 1"), ("clay", "0 1")]));

        assert_eq!(regions(&file), [vec!["rock"], vec!["rock", "clay"]]);
    }

    #[test]
    fn test_file_without_region_arrays_gives_no_region() {
        assert!(regions(&square_and_triangle_file(""))
            .iter()
            .all(Vec::is_empty));
    }

    #[test]
    fn test_region_array_with_one_value_for_two_cells_is_rejected() {
        let message = error_message(&square_and_triangle_file(&region_cell_data(&[(
            "rock", "1",
        )])));

        assert!(
            message.contains("region__rock must have one value for each cell"),
            "{message}"
        );
    }

    fn fracture_id_cell_data(values: &str) -> String {
        format!(
            r#"<CellData>
            <DataArray type="Int64" Name="fracture_id" format="ascii">{values}</DataArray>
            </CellData>"#
        )
    }

    fn fracture_ids(file: &str) -> Vec<Option<u64>> {
        parse_user_fractures(file)
            .unwrap()
            .iter()
            .map(|fracture| fracture.fracture_id)
            .collect()
    }

    #[test]
    fn test_fracture_id_array_gives_the_fracture_of_each_polygon() {
        let file = square_and_triangle_file(&fracture_id_cell_data("0 7"));

        assert_eq!(fracture_ids(&file), [Some(0), Some(7)]);
    }

    #[test]
    fn test_file_without_fracture_id_array_gives_no_fracture_id() {
        assert_eq!(fracture_ids(&square_and_triangle_file("")), [None, None]);
    }

    #[test]
    fn test_fracture_id_array_with_one_value_for_two_cells_is_rejected() {
        let message = error_message(&square_and_triangle_file(&fracture_id_cell_data("2")));

        assert!(
            message.contains("the fracture_id array must have one value for each cell"),
            "{message}"
        );
    }

    #[test]
    fn test_fracture_id_array_with_two_components_is_rejected() {
        let cell_data = fracture_id_cell_data("2 1 0 0").replace(
            r#"Name="fracture_id""#,
            r#"Name="fracture_id" NumberOfComponents="2""#,
        );
        let message = error_message(&square_and_triangle_file(&cell_data));

        assert!(
            message.contains("the fracture_id array must have one value for each fracture"),
            "{message}"
        );
    }

    #[test]
    fn test_negative_fracture_id_is_rejected() {
        let message = error_message(&square_and_triangle_file(&fracture_id_cell_data("-1 0")));

        assert!(
            message.contains("the DataArray fracture_id has the invalid value -1"),
            "{message}"
        );
    }

    #[test]
    fn test_region_of_a_polygon_is_after_the_values_of_the_line_cells() {
        let line_and_regions = format!(
            r#"<Lines>
            <DataArray type="Int64" Name="connectivity" format="ascii">0 1</DataArray>
            <DataArray type="Int64" Name="offsets" format="ascii">2</DataArray>
            </Lines>
            {}"#,
            region_cell_data(&[("rock", "1 1 0")])
        );
        let expected: [&[&str]; 2] = [&["rock"], &[]];

        assert_eq!(
            regions(&square_and_triangle_file(&line_and_regions)),
            expected
        );
    }
}
