use std::fmt::Display;
use std::io::{self, Write};

use crate::fracture::poly::Poly;

/// Writes all accepted fractures as polygons in VTK XML PolyData format (`.vtp`, ASCII)
///
/// Must be called before the polygons are rotated to the x-y plane for the
/// intersection files, because it writes the 3D vertices.
///
/// One polygon cell is written for each accepted fracture, with the cell data:
/// * `fracture_id` - Number of the fracture in the other output files (starts at 1),
///   0 for a fracture which was removed
/// * `accepted_id` - Number of the fracture in the order of acceptance (starts at 1)
/// * `removed` - 1 if the fracture was removed as isolated, or does not connect
///   the boundary faces, else 0
/// * `family` - Family number as in radii.dat: -2 user polygon, -1 user rectangle,
///   0 user ellipse, > 0 stochastic family
/// * `cluster` - Number of the connected group of fractures
/// * `truncated` - 1 if the fracture was cut by the domain, else 0
/// * `intersection_count` - Number of intersections with other accepted fractures
/// * `x_radius`, `y_radius`, `aspect_ratio` - Size before truncation
/// * `area` - Area after truncation
/// * `normal` - Unit normal vector
///
/// # Arguments
///
/// * `final_fractures` - Sorted indices of fractures left after isolated fracture removal
/// * `accepted_poly` - All accepted fractures
/// * `output` - File to write to
pub fn write_fractures_vtk(
    final_fractures: &[usize],
    accepted_poly: &[Poly],
    output: &mut impl Write,
) -> io::Result<()> {
    let point_count: usize = accepted_poly.iter().map(|poly| poly.vertices.len()).sum();

    writeln!(output, r#"<?xml version="1.0"?>"#)?;
    writeln!(
        output,
        r#"<VTKFile type="PolyData" version="1.0" byte_order="LittleEndian">"#
    )?;
    writeln!(output, "<PolyData>")?;
    writeln!(
        output,
        r#"<Piece NumberOfPoints="{}" NumberOfVerts="0" NumberOfLines="0" NumberOfStrips="0" NumberOfPolys="{}">"#,
        point_count,
        accepted_poly.len()
    )?;

    writeln!(output, "<Points>")?;
    write_data_array(
        output,
        "points",
        "Float64",
        3,
        accepted_poly
            .iter()
            .flat_map(|poly| &poly.vertices)
            .map(|vertex| format!("{} {} {}", vertex.x, vertex.y, vertex.z)),
    )?;
    writeln!(output, "</Points>")?;

    // The points are in polygon order, so the connectivity is 0, 1, 2, ...
    // and the offset of a polygon is the number of points up to its end
    writeln!(output, "<Polys>")?;
    write_data_array(output, "connectivity", "Int64", 1, 0..point_count)?;
    write_data_array(
        output,
        "offsets",
        "Int64",
        1,
        accepted_poly.iter().scan(0, |offset, poly| {
            *offset += poly.vertices.len();
            Some(*offset)
        }),
    )?;
    writeln!(output, "</Polys>")?;

    // Number of each fracture in the other output files, 0 if it was removed
    let fracture_ids: Vec<usize> = (0..accepted_poly.len())
        .map(|i| match final_fractures.binary_search(&i) {
            Ok(final_index) => final_index + 1,
            Err(_) => 0,
        })
        .collect();

    writeln!(output, r#"<CellData Normals="normal">"#)?;
    write_data_array(output, "fracture_id", "Int32", 1, fracture_ids.iter())?;
    write_data_array(output, "accepted_id", "Int32", 1, 1..=accepted_poly.len())?;
    write_data_array(
        output,
        "removed",
        "Int32",
        1,
        fracture_ids.iter().map(|&id| (id == 0) as u8),
    )?;
    write_data_array(
        output,
        "family",
        "Int32",
        1,
        accepted_poly.iter().map(|poly| poly.family_num + 1),
    )?;
    write_data_array(
        output,
        "cluster",
        "Int32",
        1,
        accepted_poly.iter().map(|poly| poly.group_num),
    )?;
    write_data_array(
        output,
        "truncated",
        "Int32",
        1,
        accepted_poly.iter().map(|poly| poly.truncated as u8),
    )?;
    write_data_array(
        output,
        "intersection_count",
        "Int32",
        1,
        accepted_poly
            .iter()
            .map(|poly| poly.intersection_index.len()),
    )?;

    let mut write_real_property = |name: &str, value_of: fn(&Poly) -> f64| {
        write_data_array(
            output,
            name,
            "Float64",
            1,
            accepted_poly.iter().map(value_of),
        )
    };
    write_real_property("x_radius", |poly| poly.xradius)?;
    write_real_property("y_radius", |poly| poly.yradius)?;
    write_real_property("aspect_ratio", |poly| poly.aspect_ratio)?;
    write_real_property("area", |poly| poly.area)?;

    write_data_array(
        output,
        "normal",
        "Float64",
        3,
        accepted_poly
            .iter()
            .map(|poly| format!("{} {} {}", poly.normal.x, poly.normal.y, poly.normal.z)),
    )?;
    writeln!(output, "</CellData>")?;

    writeln!(output, "</Piece>")?;
    writeln!(output, "</PolyData>")?;
    writeln!(output, "</VTKFile>")
}

/// Writes one ASCII `DataArray` element, one tuple on each line.
fn write_data_array<T: Display>(
    output: &mut impl Write,
    name: &str,
    vtk_type: &str,
    number_of_components: usize,
    tuples: impl Iterator<Item = T>,
) -> io::Result<()> {
    writeln!(
        output,
        r#"<DataArray type="{}" Name="{}" NumberOfComponents="{}" format="ascii">"#,
        vtk_type, name, number_of_components
    )?;

    for tuple in tuples {
        writeln!(output, "{}", tuple)?;
    }

    writeln!(output, "</DataArray>")
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Writes a rectangle and an ellipse with 6 nodes. Only the ellipse is a final fracture.
    fn write_two_fractures() -> String {
        let mut rectangle = Poly::new_rect(1., 2.);
        rectangle.family_num = -2;
        let mut ellipse = Poly::new_ell(6, 3., 1.);
        ellipse.family_num = 0;
        ellipse.area = 20.;

        let mut output = Vec::new();
        write_fractures_vtk(&[1], &[rectangle, ellipse], &mut output).unwrap();
        String::from_utf8(output).unwrap()
    }

    /// The lines between the start tag of the named data array and its end tag.
    fn data_array<'a>(file: &'a str, name: &str) -> Vec<&'a str> {
        file.lines()
            .skip_while(|line| !line.contains(&format!(r#"Name="{}""#, name)))
            .skip(1)
            .take_while(|line| *line != "</DataArray>")
            .collect()
    }

    #[test]
    fn test_each_accepted_fracture_is_one_polygon_cell() {
        let file = write_two_fractures();

        assert!(file.contains(r#"NumberOfPoints="10""#));
        assert!(file.contains(r#"NumberOfPolys="2""#));
        assert_eq!(data_array(&file, "points").len(), 10);
        assert_eq!(data_array(&file, "offsets"), ["4", "10"]);
        assert_eq!(data_array(&file, "connectivity").last(), Some(&"9"));
    }

    #[test]
    fn test_points_are_the_polygon_vertices() {
        let file = write_two_fractures();

        assert_eq!(data_array(&file, "points")[0], "1 2 0");
    }

    #[test]
    fn test_removed_fracture_has_no_fracture_id() {
        let file = write_two_fractures();

        assert_eq!(data_array(&file, "fracture_id"), ["0", "1"]);
        assert_eq!(data_array(&file, "removed"), ["1", "0"]);
    }

    #[test]
    fn test_family_numbers_match_the_radii_files() {
        let file = write_two_fractures();

        assert_eq!(data_array(&file, "family"), ["-1", "1"]);
    }

    #[test]
    fn test_properties_are_written_for_each_fracture() {
        let file = write_two_fractures();

        assert_eq!(data_array(&file, "x_radius"), ["1", "3"]);
        assert_eq!(data_array(&file, "area"), ["0", "20"]);
        assert_eq!(data_array(&file, "normal"), ["0 0 1", "0 0 1"]);
    }
}
