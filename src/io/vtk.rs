use std::fmt::Display;
use std::io::{self, Write};

use parry3d_f64::na::Vector3;

use crate::fracture::poly::Poly;
use crate::structures::RejectedFracture;

/// Writes all accepted fractures as polygons in VTK XML PolyData format (`.vtp`, ASCII)
///
/// One polygon cell is written for each accepted fracture, with the cell data:
/// * `fracture_id` - Number of the fracture in the other output files (starts at 1),
///   0 for a fracture which was removed
/// * `accepted_id` - Number of the fracture in the order of acceptance (starts at 1)
/// * `removed` - 1 if the fracture was removed as isolated, or does not connect
///   the boundary faces, else 0
/// * `family` - Family number: 0 user-defined fracture, > 0 stochastic family
/// * `cluster` - Number of the connected group of fractures
/// * `truncated` - 1 if the fracture was cut by the domain, else 0
/// * `intersection_count` - Number of intersections with other accepted fractures
/// * `boundary_face_count` - Number of domain boundary faces which the fracture meets
/// * `x_radius`, `y_radius`, `aspect_ratio` - Size before truncation
/// * `area` - Area after truncation
/// * `normal` - Unit normal vector
/// * `translation` - Translation of the fracture center from the origin
///
/// # Arguments
///
/// * `final_fractures` - Sorted indices of fractures left after isolated fracture removal
/// * `accepted_poly` - All accepted fractures
/// * `boundary_face_counts` - `boundary_face_count` of each accepted fracture
/// * `output` - File to write to
pub fn write_fractures_vtk(
    final_fractures: &[usize],
    accepted_poly: &[Poly],
    boundary_face_counts: &[usize],
    output: &mut impl Write,
) -> io::Result<()> {
    let polys: Vec<&Poly> = accepted_poly.iter().collect();
    write_geometry(output, &polys)?;

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
    write_data_array(
        output,
        "boundary_face_count",
        "Int32",
        1,
        boundary_face_counts.iter(),
    )?;
    write_real_array(output, "x_radius", &polys, |poly| poly.xradius)?;
    write_real_array(output, "y_radius", &polys, |poly| poly.yradius)?;
    write_real_array(output, "aspect_ratio", &polys, |poly| poly.aspect_ratio)?;
    write_real_array(output, "area", &polys, |poly| poly.area)?;
    write_vector_array(output, "normal", &polys, |poly| poly.normal)?;
    write_vector_array(output, "translation", &polys, |poly| poly.translation)?;
    writeln!(output, "</CellData>")?;

    write_file_end(output)
}

/// Writes the rejected fractures as polygons in VTK XML PolyData format (`.vtp`, ASCII)
///
/// One polygon cell is written for each rejected fracture, at the position of its
/// last attempt, with the cell data:
/// * `family` - Family number: 0 user-defined fracture, > 0 stochastic family
/// * `reject_code` - FRAM reject code of the last attempt
/// * `truncated` - 1 if the fracture was cut by the domain, else 0
/// * `x_radius`, `y_radius`, `aspect_ratio` - Size before truncation
/// * `normal` - Unit normal vector
/// * `translation` - Translation of the fracture center from the origin
///
/// # Arguments
///
/// * `rejected_poly` - All rejected fractures
/// * `output` - File to write to
pub fn write_rejected_fractures_vtk(
    rejected_poly: &[RejectedFracture],
    output: &mut impl Write,
) -> io::Result<()> {
    let polys: Vec<&Poly> = rejected_poly
        .iter()
        .map(|rejected| &rejected.poly)
        .collect();
    write_geometry(output, &polys)?;

    writeln!(output, r#"<CellData Normals="normal">"#)?;
    write_data_array(
        output,
        "family",
        "Int32",
        1,
        polys.iter().map(|poly| poly.family_num + 1),
    )?;
    write_data_array(
        output,
        "reject_code",
        "Int32",
        1,
        rejected_poly.iter().map(|rejected| rejected.reject_code),
    )?;
    write_data_array(
        output,
        "truncated",
        "Int32",
        1,
        polys.iter().map(|poly| poly.truncated as u8),
    )?;
    write_real_array(output, "x_radius", &polys, |poly| poly.xradius)?;
    write_real_array(output, "y_radius", &polys, |poly| poly.yradius)?;
    write_real_array(output, "aspect_ratio", &polys, |poly| poly.aspect_ratio)?;
    write_vector_array(output, "normal", &polys, |poly| poly.normal)?;
    write_vector_array(output, "translation", &polys, |poly| poly.translation)?;
    writeln!(output, "</CellData>")?;

    write_file_end(output)
}

/// Writes the start of the file, and one polygon cell for each fracture.
fn write_geometry(output: &mut impl Write, polys: &[&Poly]) -> io::Result<()> {
    let point_count: usize = polys.iter().map(|poly| poly.vertices.len()).sum();

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
        polys.len()
    )?;

    writeln!(output, "<Points>")?;
    write_data_array(
        output,
        "points",
        "Float64",
        3,
        polys
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
        polys.iter().scan(0, |offset, poly| {
            *offset += poly.vertices.len();
            Some(*offset)
        }),
    )?;
    writeln!(output, "</Polys>")
}

/// Closes the elements which `write_geometry` opened.
fn write_file_end(output: &mut impl Write) -> io::Result<()> {
    writeln!(output, "</Piece>")?;
    writeln!(output, "</PolyData>")?;
    writeln!(output, "</VTKFile>")
}

/// Writes one real number for each fracture.
fn write_real_array(
    output: &mut impl Write,
    name: &str,
    polys: &[&Poly],
    value_of: fn(&Poly) -> f64,
) -> io::Result<()> {
    write_data_array(
        output,
        name,
        "Float64",
        1,
        polys.iter().map(|poly| value_of(poly)),
    )
}

/// Writes one vector for each fracture.
fn write_vector_array(
    output: &mut impl Write,
    name: &str,
    polys: &[&Poly],
    value_of: fn(&Poly) -> Vector3<f64>,
) -> io::Result<()> {
    write_data_array(
        output,
        name,
        "Float64",
        3,
        polys.iter().map(|poly| {
            let vector = value_of(poly);
            format!("{} {} {}", vector.x, vector.y, vector.z)
        }),
    )
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
    use crate::fracture::user_defined_fractures::USER_FRACTURE_FAMILY;

    /// Writes a rectangle and an ellipse with 6 nodes. Only the ellipse is a final fracture.
    fn write_two_fractures() -> String {
        let mut rectangle = Poly::new(4, 1., 2.);
        rectangle.family_num = USER_FRACTURE_FAMILY;
        let mut ellipse = Poly::new(6, 3., 1.);
        ellipse.family_num = 0;
        ellipse.area = 20.;
        ellipse.translation = Vector3::new(1., 2., 3.);

        let mut output = Vec::new();
        write_fractures_vtk(&[1], &[rectangle, ellipse], &[0, 2], &mut output).unwrap();
        String::from_utf8(output).unwrap()
    }

    /// Writes one rejected rectangle of the first stochastic family.
    fn write_rejected_rectangle() -> String {
        let mut rectangle = Poly::new(4, 1., 2.);
        rectangle.family_num = 0;
        let rejected = RejectedFracture {
            poly: rectangle,
            reject_code: -6,
        };

        let mut output = Vec::new();
        write_rejected_fractures_vtk(&[rejected], &mut output).unwrap();
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

        let first_vertex = Poly::new(4, 1., 2.).vertices[0];
        assert_eq!(
            data_array(&file, "points")[0],
            format!("{} {} 0", first_vertex.x, first_vertex.y)
        );
    }

    #[test]
    fn test_removed_fracture_has_no_fracture_id() {
        let file = write_two_fractures();

        assert_eq!(data_array(&file, "fracture_id"), ["0", "1"]);
        assert_eq!(data_array(&file, "removed"), ["1", "0"]);
    }

    #[test]
    fn test_family_number_is_one_more_than_the_internal_number() {
        let file = write_two_fractures();

        assert_eq!(data_array(&file, "family"), ["0", "1"]);
    }

    #[test]
    fn test_properties_are_written_for_each_fracture() {
        let file = write_two_fractures();

        assert_eq!(data_array(&file, "x_radius"), ["1", "3"]);
        assert_eq!(data_array(&file, "area"), ["0", "20"]);
        assert_eq!(data_array(&file, "boundary_face_count"), ["0", "2"]);
        assert_eq!(data_array(&file, "normal"), ["0 0 1", "0 0 1"]);
    }

    #[test]
    fn test_translation_is_written_for_each_fracture() {
        let file = write_two_fractures();

        assert_eq!(data_array(&file, "translation"), ["0 0 0", "1 2 3"]);
    }

    #[test]
    fn test_each_rejected_fracture_is_one_polygon_cell() {
        let file = write_rejected_rectangle();

        assert!(file.contains(r#"NumberOfPolys="1""#));
        assert_eq!(data_array(&file, "offsets"), ["4"]);
    }

    #[test]
    fn test_rejected_fracture_has_its_reject_code_and_family() {
        let file = write_rejected_rectangle();

        assert_eq!(data_array(&file, "reject_code"), ["-6"]);
        assert_eq!(data_array(&file, "family"), ["1"]);
    }

    #[test]
    fn test_rejected_fracture_has_no_fracture_id() {
        let file = write_rejected_rectangle();

        assert!(!file.contains(r#"Name="fracture_id""#));
    }
}
