use std::fmt::Display;
use std::io::{self, Write};

use parry3d_f64::na::Vector3;

use crate::fracture::poly::Poly;
use crate::io::output::to_output_coordinates;
use crate::structures::RejectedFracture;

/// Writes the final fractures as polygons in VTK XML PolyData format (`.vtp`, ASCII)
///
/// One polygon cell is written for each final fracture, in the order of the
/// fracture numbers, with the cell data:
/// * `fracture_id` - Number of the fracture in the other output files (starts at 1)
/// * the cell data of `write_accepted_cell_data`
///
/// # Arguments
///
/// * `final_fractures` - Sorted indices of fractures left after isolated fracture removal
/// * `accepted_poly` - All accepted fractures
/// * `mesh_center` - Center of the domain mesh
/// * `output` - File to write to
pub fn write_fractures_vtk(
    final_fractures: &[usize],
    accepted_poly: &[Poly],
    mesh_center: &Vector3<f64>,
    output: &mut impl Write,
) -> io::Result<()> {
    let polys: Vec<&Poly> = final_fractures.iter().map(|&i| &accepted_poly[i]).collect();
    write_geometry(output, &polys, mesh_center)?;

    writeln!(output, r#"<CellData Normals="normal">"#)?;
    write_data_array(output, "fracture_id", "Int32", 1, 1..=final_fractures.len())?;
    write_accepted_cell_data(output, final_fractures, accepted_poly, mesh_center)?;
    writeln!(output, "</CellData>")?;

    write_file_end(output)
}

/// Writes the removed fractures as polygons in VTK XML PolyData format (`.vtp`, ASCII)
///
/// A removed fracture is an accepted fracture which is isolated, or does not
/// connect the boundary faces. One polygon cell is written for each of them, with
/// the cell data of `write_accepted_cell_data`.
///
/// # Arguments
///
/// * `removed_fractures` - Sorted indices of the accepted fractures which were removed
/// * `accepted_poly` - All accepted fractures
/// * `mesh_center` - Center of the domain mesh
/// * `output` - File to write to
pub fn write_removed_fractures_vtk(
    removed_fractures: &[usize],
    accepted_poly: &[Poly],
    mesh_center: &Vector3<f64>,
    output: &mut impl Write,
) -> io::Result<()> {
    let polys: Vec<&Poly> = removed_fractures
        .iter()
        .map(|&i| &accepted_poly[i])
        .collect();
    write_geometry(output, &polys, mesh_center)?;

    writeln!(output, r#"<CellData Normals="normal">"#)?;
    write_accepted_cell_data(output, removed_fractures, accepted_poly, mesh_center)?;
    writeln!(output, "</CellData>")?;

    write_file_end(output)
}

/// Writes the cell data which the final and the removed fractures have in common
///
/// * `accepted_id` - Number of the fracture in the order of acceptance (starts at 1)
/// * `family` - Family number: 0 user-defined fracture, > 0 stochastic family
/// * `user_fracture_id` - Value of the `fracture_id` array of the user-defined fracture file, -1 if none
/// * `cluster` - Number of the connected group of fractures
/// * `intersection_count` - Number of intersections with other accepted fractures
/// * `x_radius`, `y_radius`, `aspect_ratio` - Size
/// * `area` - Area
/// * `normal` - Unit normal vector
/// * `translation` - Translation of the fracture center from the origin
///
/// # Arguments
///
/// * `output` - File to write to
/// * `fractures` - Indices of the accepted fractures to write
/// * `accepted_poly` - All accepted fractures
/// * `mesh_center` - Center of the domain mesh
fn write_accepted_cell_data(
    output: &mut impl Write,
    fractures: &[usize],
    accepted_poly: &[Poly],
    mesh_center: &Vector3<f64>,
) -> io::Result<()> {
    let polys: Vec<&Poly> = fractures.iter().map(|&i| &accepted_poly[i]).collect();

    write_data_array(
        output,
        "accepted_id",
        "Int32",
        1,
        fractures.iter().map(|&i| i + 1),
    )?;
    write_data_array(
        output,
        "family",
        "Int32",
        1,
        polys.iter().map(|poly| poly.family_num + 1),
    )?;
    write_user_fracture_id_array(output, &polys)?;
    write_data_array(
        output,
        "cluster",
        "Int32",
        1,
        polys.iter().map(|poly| poly.group_num),
    )?;
    write_data_array(
        output,
        "intersection_count",
        "Int32",
        1,
        polys.iter().map(|poly| poly.intersection_index.len()),
    )?;
    write_real_array(output, "x_radius", &polys, |poly| poly.xradius)?;
    write_real_array(output, "y_radius", &polys, |poly| poly.yradius)?;
    write_real_array(output, "aspect_ratio", &polys, |poly| poly.aspect_ratio)?;
    write_real_array(output, "area", &polys, |poly| poly.area)?;
    write_vector_array(output, "normal", &polys, |poly| poly.normal)?;
    write_translation_array(output, &polys, mesh_center)
}

/// Writes the rejected fractures as polygons in VTK XML PolyData format (`.vtp`, ASCII)
///
/// One polygon cell is written for each rejected fracture, at the position of its
/// last attempt, with the cell data:
/// * `family` - Family number: 0 user-defined fracture, > 0 stochastic family
/// * `user_fracture_id` - Value of the `fracture_id` array of the user-defined fracture file, -1 if none
/// * `reject_code` - FRAM reject code of the last attempt
/// * `x_radius`, `y_radius`, `aspect_ratio` - Size
/// * `normal` - Unit normal vector
/// * `translation` - Translation of the fracture center from the origin
///
/// # Arguments
///
/// * `rejected_poly` - All rejected fractures
/// * `mesh_center` - Center of the domain mesh
/// * `output` - File to write to
pub fn write_rejected_fractures_vtk(
    rejected_poly: &[RejectedFracture],
    mesh_center: &Vector3<f64>,
    output: &mut impl Write,
) -> io::Result<()> {
    let polys: Vec<&Poly> = rejected_poly
        .iter()
        .map(|rejected| &rejected.poly)
        .collect();
    write_geometry(output, &polys, mesh_center)?;

    writeln!(output, r#"<CellData Normals="normal">"#)?;
    write_data_array(
        output,
        "family",
        "Int32",
        1,
        polys.iter().map(|poly| poly.family_num + 1),
    )?;
    write_user_fracture_id_array(output, &polys)?;
    write_data_array(
        output,
        "reject_code",
        "Int32",
        1,
        rejected_poly.iter().map(|rejected| rejected.reject_code),
    )?;
    write_real_array(output, "x_radius", &polys, |poly| poly.xradius)?;
    write_real_array(output, "y_radius", &polys, |poly| poly.yradius)?;
    write_real_array(output, "aspect_ratio", &polys, |poly| poly.aspect_ratio)?;
    write_vector_array(output, "normal", &polys, |poly| poly.normal)?;
    write_translation_array(output, &polys, mesh_center)?;
    writeln!(output, "</CellData>")?;

    write_file_end(output)
}

/// Writes the start of the file, and one polygon cell for each fracture.
fn write_geometry(
    output: &mut impl Write,
    polys: &[&Poly],
    mesh_center: &Vector3<f64>,
) -> io::Result<()> {
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
            .map(|vertex| to_output_coordinates(vertex, mesh_center))
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

/// Writes the translation of each fracture, in the coordinates of the output.
fn write_translation_array(
    output: &mut impl Write,
    polys: &[&Poly],
    mesh_center: &Vector3<f64>,
) -> io::Result<()> {
    write_vector_array(output, "translation", polys, |poly| {
        to_output_coordinates(&poly.translation.into(), mesh_center).coords
    })
}

/// Writes one vector for each fracture.
fn write_vector_array(
    output: &mut impl Write,
    name: &str,
    polys: &[&Poly],
    value_of: impl Fn(&Poly) -> Vector3<f64>,
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

/// Writes the cell data array `user_fracture_id`: the value of the `fracture_id`
/// array of the user-defined fracture file, -1 for a fracture without it.
fn write_user_fracture_id_array(output: &mut impl Write, polys: &[&Poly]) -> io::Result<()> {
    write_data_array(
        output,
        "user_fracture_id",
        "Int64",
        1,
        polys
            .iter()
            .map(|poly| poly.user_fracture_id.map_or(-1, |id| id as i64)),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fracture::user_defined_fractures::USER_FRACTURE_FAMILY;

    /// A rectangle, an ellipse with 6 nodes and an ellipse with 8 nodes.
    /// The rectangle was removed; the ellipses are final fractures.
    fn three_accepted_fractures() -> Vec<Poly> {
        let mut rectangle = Poly::new(4, 1., 2.);
        rectangle.family_num = USER_FRACTURE_FAMILY;
        rectangle.user_fracture_id = Some(7);
        let mut small_ellipse = Poly::new(6, 3., 1.);
        small_ellipse.family_num = 0;
        small_ellipse.area = 20.;
        small_ellipse.translation = Vector3::new(1., 2., 3.);
        let mut large_ellipse = Poly::new(8, 5., 1.);
        large_ellipse.family_num = 1;

        vec![rectangle, small_ellipse, large_ellipse]
    }

    fn write_final_fractures() -> String {
        write_final_fractures_of_mesh(&Vector3::zeros())
    }

    fn write_final_fractures_of_mesh(mesh_center: &Vector3<f64>) -> String {
        let accepted_poly = three_accepted_fractures();

        let mut output = Vec::new();
        write_fractures_vtk(&[1, 2], &accepted_poly, mesh_center, &mut output).unwrap();
        String::from_utf8(output).unwrap()
    }

    fn write_removed_fractures() -> String {
        let accepted_poly = three_accepted_fractures();

        let mut output = Vec::new();
        write_removed_fractures_vtk(&[0], &accepted_poly, &Vector3::zeros(), &mut output).unwrap();
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
        write_rejected_fractures_vtk(&[rejected], &Vector3::zeros(), &mut output).unwrap();
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
    fn test_each_final_fracture_is_one_polygon_cell() {
        let file = write_final_fractures();

        assert!(file.contains(r#"NumberOfPoints="14""#));
        assert!(file.contains(r#"NumberOfPolys="2""#));
        assert_eq!(data_array(&file, "points").len(), 14);
        assert_eq!(data_array(&file, "offsets"), ["6", "14"]);
        assert_eq!(data_array(&file, "connectivity").last(), Some(&"13"));
    }

    #[test]
    fn test_points_are_the_polygon_vertices() {
        let file = write_final_fractures();

        let first_vertex = Poly::new(6, 3., 1.).vertices[0];
        assert_eq!(
            data_array(&file, "points")[0],
            format!("{} {} 0", first_vertex.x, first_vertex.y)
        );
    }

    #[test]
    fn test_fracture_id_is_the_cell_number() {
        let file = write_final_fractures();

        assert_eq!(data_array(&file, "fracture_id"), ["1", "2"]);
    }

    #[test]
    fn test_final_fracture_keeps_its_accepted_id() {
        let file = write_final_fractures();

        assert_eq!(data_array(&file, "accepted_id"), ["2", "3"]);
    }

    #[test]
    fn test_fractures_file_has_no_removed_array() {
        let file = write_final_fractures();

        assert!(!file.contains(r#"Name="removed""#));
    }

    #[test]
    fn test_family_number_is_one_more_than_the_internal_number() {
        let file = write_final_fractures();

        assert_eq!(data_array(&file, "family"), ["1", "2"]);
    }

    #[test]
    fn test_properties_are_written_for_each_fracture() {
        let file = write_final_fractures();

        assert_eq!(data_array(&file, "x_radius"), ["3", "5"]);
        assert_eq!(data_array(&file, "area"), ["20", "0"]);
        assert_eq!(data_array(&file, "normal"), ["0 0 1", "0 0 1"]);
    }

    #[test]
    fn test_points_with_mesh_have_the_mesh_coordinates() {
        let file = write_final_fractures_of_mesh(&Vector3::new(10., 20., 30.));
        // The first final fracture has the translation (1, 2, 3) only in its cell data
        let first_vertex = Poly::new(6, 3., 1.).vertices[0];

        assert_eq!(
            data_array(&file, "points")[0],
            format!("{} {} 30", first_vertex.x + 10., first_vertex.y + 20.)
        );
    }

    #[test]
    fn test_translation_with_mesh_has_the_mesh_coordinates() {
        let file = write_final_fractures_of_mesh(&Vector3::new(10., 20., 30.));

        assert_eq!(data_array(&file, "translation"), ["11 22 33", "10 20 30"]);
    }

    #[test]
    fn test_translation_is_written_for_each_fracture() {
        let file = write_final_fractures();

        assert_eq!(data_array(&file, "translation"), ["1 2 3", "0 0 0"]);
    }

    #[test]
    fn test_each_removed_fracture_is_one_polygon_cell() {
        let file = write_removed_fractures();

        assert!(file.contains(r#"NumberOfPolys="1""#));
        assert_eq!(data_array(&file, "offsets"), ["4"]);
    }

    #[test]
    fn test_removed_fracture_keeps_its_accepted_id() {
        let file = write_removed_fractures();

        assert_eq!(data_array(&file, "accepted_id"), ["1"]);
    }

    #[test]
    fn test_removed_user_fracture_has_its_user_fracture_id() {
        let file = write_removed_fractures();

        assert_eq!(data_array(&file, "user_fracture_id"), ["7"]);
    }

    #[test]
    fn test_fracture_without_user_fracture_id_has_the_value_minus_one() {
        let file = write_final_fractures();

        assert_eq!(data_array(&file, "user_fracture_id"), ["-1", "-1"]);
    }

    #[test]
    fn test_removed_fractures_file_has_no_fracture_id() {
        let file = write_removed_fractures();

        assert!(!file.contains(r#" Name="fracture_id""#));
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

        assert!(!file.contains(r#" Name="fracture_id""#));
    }
}
