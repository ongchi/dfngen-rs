//! One test case of `tests/data` for each feature of the input. Each test
//! checks what its feature must do with the network.

use std::cell::RefCell;
use std::fs;
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::rc::Rc;

use dfngen::fracture::fracture_family::FractureFamilyCollection;
use dfngen::fracture::poly::{FractureOrigin, Poly};
use dfngen::io::input::{Input, StopCondition};
use dfngen::mesh::FractureMesh;
use dfngen::region_boundaries::{BoundarySegment, OtherSide, SegmentSide};
use dfngen::stats::RejectionReasons;
use dfngen::{
    build_input_from_text, generate, read_input, read_input_file_names, write_output, Network,
};
use parry3d_f64::na::{Point3, Vector3};

fn case_folder(case: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(case)
}

/// Network of a case
fn run(case: &str) -> Network {
    let input_file = case_folder(case).join(format!("{case}.toml"));
    let (input, families) = read_input(input_file.to_str().unwrap()).unwrap();

    generate(input, families, &mut |_| {}).unwrap()
}

/// Network of a case, with one text of its input file replaced
fn run_with_replaced_text(case: &str, from: &str, to: &str) -> Network {
    run_with_replacements(case, &[(from, to)])
}

/// Network of a case, with texts of its input file replaced (from, to)
fn run_with_replacements(case: &str, replacements: &[(&str, &str)]) -> Network {
    let (input, families) = input_with_replacements(case, replacements);

    generate(input, families, &mut |_| {}).unwrap()
}

/// Input of a case, with texts of its input file replaced (from, to)
fn input_with_replacements(
    case: &str,
    replacements: &[(&str, &str)],
) -> (Input, FractureFamilyCollection) {
    let folder = case_folder(case);
    let mut toml_text = fs::read_to_string(folder.join(format!("{case}.toml"))).unwrap();

    for (from, to) in replacements {
        assert!(
            toml_text.contains(from),
            "the input of {case} has no {from:?}"
        );
        toml_text = toml_text.replace(from, to);
    }

    let names = read_input_file_names(&toml_text).unwrap();
    let read_text = |name: &str| fs::read_to_string(folder.join(name)).unwrap();
    let mesh_text = read_text(&names.domain);
    let fractures_text = names.user_fractures.as_deref().map(read_text);

    build_input_from_text(&toml_text, &mesh_text, fractures_text.as_deref()).unwrap()
}

/// Names of the files which `write_output` writes for the network
fn output_file_names(network: &Network) -> Vec<String> {
    let mut names = Vec::new();
    write_output(network, &mut |name| {
        names.push(name.to_string());
        Ok(Box::new(io::sink()) as Box<dyn Write>)
    })
    .unwrap();

    names
}

/// Accepted fractures of a stochastic family, by its index in the input file
fn family_fractures(network: &Network, family_index: usize) -> Vec<&Poly> {
    let fractures: Vec<&Poly> = network
        .dfn
        .accepted_polys
        .iter()
        .filter(|poly| poly.origin == FractureOrigin::Family(family_index))
        .collect();
    assert!(
        !fractures.is_empty(),
        "family {family_index} has no accepted fracture"
    );

    fractures
}

fn center(poly: &Poly) -> Point3<f64> {
    Point3::from(poly.translation)
}

/// True if the center of each fracture of the family is in a cell of the
/// region with the name, or of the whole mesh for `None`
fn are_centers_in_region(network: &Network, family_index: usize, region: Option<&str>) -> bool {
    let mesh = &network.input.domain_mesh;
    let region_id = region.map_or(0, |name| {
        1 + mesh
            .regions()
            .iter()
            .position(|mesh_region| mesh_region.name == name)
            .unwrap()
    });

    family_fractures(network, family_index)
        .iter()
        .all(|poly| mesh.contains_point(&center(poly), region_id))
}

/// Number of fractures which FRAM rejected (the rejects outside the domain are not included)
fn fram_reject_count(reasons: &RejectionReasons) -> usize {
    reasons.short_intersection
        + reasons.close_to_node
        + reasons.close_to_edge
        + reasons.close_point_to_edge
        + reasons.triple
        + reasons.inter_close_to_inter
}

fn angle_degrees(a: &Vector3<f64>, b: &Vector3<f64>) -> f64 {
    a.normalize()
        .dot(&b.normalize())
        .clamp(-1., 1.)
        .acos()
        .to_degrees()
}

/// Mean of unit vectors which point to the same side as `side`
fn mean_direction(
    vectors: impl Iterator<Item = Vector3<f64>>,
    side: &Vector3<f64>,
) -> Vector3<f64> {
    vectors
        .map(|vector| {
            let unit = vector.normalize();
            if unit.dot(side) < 0. {
                -unit
            } else {
                unit
            }
        })
        .sum::<Vector3<f64>>()
        .normalize()
}

mod base {
    use super::*;

    #[test]
    fn test_n_poly_stop_condition_accepts_n_poly_fractures() {
        let network = run("base");

        assert_eq!(
            network.input.stop_condition,
            StopCondition::NPoly(network.dfn.accepted_polys.len())
        );
    }

    #[test]
    fn test_each_final_fracture_has_an_intersection() {
        let network = run("base");

        assert!(network
            .final_fractures
            .iter()
            .all(|&index| !network.dfn.accepted_polys[index]
                .intersection_indices
                .is_empty()));
    }

    #[test]
    fn test_isolated_fractures_are_removed() {
        let network = run("base");

        assert!(network.final_fractures.len() < network.dfn.accepted_polys.len());
    }
}

mod stop_p32 {
    use super::*;

    #[test]
    fn test_each_family_reaches_its_p32_target() {
        let network = run("stop_p32");

        for family in &network.families.families {
            assert!(family.current_p32 >= family.p32_target);
        }
    }

    #[test]
    fn test_p32_of_a_family_is_twice_its_area_per_volume() {
        let network = run("stop_p32");

        for (index, family) in network.families.families.iter().enumerate() {
            let area: f64 = family_fractures(&network, index)
                .iter()
                .map(|poly| poly.area)
                .sum();
            let p32 = 2. * area / network.input.domain_volume;
            assert!((p32 - family.current_p32).abs() < 1e-9 * p32);
        }
    }

    #[test]
    fn test_each_family_stops_at_the_fracture_which_reaches_its_target() {
        let network = run("stop_p32");

        for (index, family) in network.families.families.iter().enumerate() {
            let last_fracture = family_fractures(&network, index).pop().unwrap();
            let p32_without_last =
                family.current_p32 - 2. * last_fracture.area / network.input.domain_volume;
            assert!(p32_without_last < family.p32_target);
        }
    }
}

mod radius_distributions {
    use super::*;

    #[test]
    fn test_each_radius_is_in_the_limits_of_its_distribution() {
        let network = run("radius_distributions");

        for family in 0..5 {
            for poly in family_fractures(&network, family) {
                assert!(
                    (1.0..=3.0).contains(&poly.x_radius),
                    "family {family}: {}",
                    poly.x_radius
                );
            }
        }
    }

    #[test]
    fn test_constant_distribution_gives_its_value() {
        let network = run("radius_distributions");

        assert!(family_fractures(&network, 4)
            .iter()
            .all(|poly| poly.x_radius == 2.0));
    }
}

mod orientation {
    use super::*;

    #[test]
    fn test_the_three_orientation_forms_give_the_same_pole() {
        let network = run("orientation");
        let pole = Vector3::new(0.5, 0., 0.75_f64.sqrt());

        for family in 0..3 {
            let normals = family_fractures(&network, family)
                .into_iter()
                .map(|poly| poly.normal);
            let mean_normal = mean_direction(normals, &pole);
            assert!(angle_degrees(&mean_normal, &pole) < 5., "family {family}");
        }
    }

    /// Largest angle between the direction to the first vertex of a fracture and its mean
    fn first_vertex_spread(network: &Network, family: usize) -> f64 {
        let directions: Vec<Vector3<f64>> = family_fractures(network, family)
            .iter()
            .map(|poly| poly.vertices[0] - center(poly))
            .collect();
        let mean = mean_direction(directions.iter().copied(), &directions[0]);

        directions
            .iter()
            .map(|direction| angle_degrees(direction, &mean))
            .fold(0., f64::max)
    }

    #[test]
    fn test_beta_fixes_the_rotation_in_the_plane() {
        let network = run("orientation");

        assert!(first_vertex_spread(&network, 2) < 30.);
    }

    #[test]
    fn test_rotation_in_the_plane_is_random_without_beta() {
        let network = run("orientation");

        assert!(first_vertex_spread(&network, 0) > 90.);
    }
}

mod fracture_shape {
    use super::*;

    #[test]
    fn test_fracture_has_the_number_of_points_of_its_family() {
        let network = run("fracture_shape");

        for (family, number_of_points) in [(0, 3), (1, 4), (2, 12)] {
            assert!(family_fractures(&network, family)
                .iter()
                .all(|poly| poly.vertices.len() == number_of_points));
        }
    }

    #[test]
    fn test_rectangle_sides_are_radius_and_aspect_ratio_times_square_root_of_2() {
        let network = run("fracture_shape");

        for poly in family_fractures(&network, 1) {
            let v = &poly.vertices;
            let mut sides = [(v[1] - v[0]).norm(), (v[2] - v[1]).norm()];
            sides.sort_by(f64::total_cmp);
            // Radius 2, aspect ratio 0.5
            assert!((sides[0] - 2_f64.sqrt()).abs() < 1e-9);
            assert!((sides[1] - 2. * 2_f64.sqrt()).abs() < 1e-9);
        }
    }

    #[test]
    fn test_y_radius_is_aspect_ratio_times_x_radius() {
        let network = run("fracture_shape");

        for (family, aspect_ratio) in [(0, 1.0), (1, 0.5), (2, 2.0)] {
            for poly in family_fractures(&network, family) {
                assert!((poly.y_radius - aspect_ratio * poly.x_radius).abs() < 1e-12);
            }
        }
    }
}

mod mesh_domain {
    use super::*;

    #[test]
    fn test_fracture_center_is_in_a_cell_of_the_region_of_its_family() {
        let network = run("mesh_domain");

        assert!(are_centers_in_region(&network, 0, None));
        assert!(are_centers_in_region(&network, 1, Some("north")));
    }

    #[test]
    fn test_user_fracture_which_intersects_no_cell_is_rejected() {
        let network = run("mesh_domain");

        assert_eq!(network.dfn.stats.rejected_user_fractures, [2]);
    }

    #[test]
    fn test_n_poly_does_not_include_the_user_fractures() {
        let network = run("mesh_domain");
        let StopCondition::NPoly(n_poly) = network.input.stop_condition else {
            panic!("the case must set generation.n_poly");
        };
        let stochastic_count = network
            .dfn
            .accepted_polys
            .iter()
            .filter(|poly| poly.origin != FractureOrigin::UserDefined)
            .count();

        assert_eq!(stochastic_count, n_poly);
    }
}

mod mesh_cell_types {
    use super::*;

    #[test]
    fn test_domain_volume_is_the_sum_of_the_cell_volumes() {
        let network = run("mesh_cell_types");

        assert!((network.input.domain_volume - 24.).abs() < 1e-9);
    }

    #[test]
    fn test_each_cell_type_gets_fracture_centers() {
        let network = run("mesh_cell_types");
        let centers: Vec<Point3<f64>> = network.dfn.accepted_polys.iter().map(center).collect();

        assert!(centers
            .iter()
            .all(|point| network.input.domain_mesh.contains_point(point, 0)));
        // The mesh is from x = 10 to 20 with its center at x = 15. Cell N is
        // from x = -5 + 2 N to -3 + 2 N in the coordinates of the generation.
        for cell in 0..5 {
            let start = -5. + 2. * f64::from(cell);
            assert!(
                centers
                    .iter()
                    .any(|point| (start..start + 2.).contains(&point.x)),
                "cell {cell} has no fracture center"
            );
        }
    }
}

mod family_regions {
    use super::*;

    #[test]
    fn test_fracture_center_is_in_a_cell_of_the_region_of_its_family() {
        let network = run("family_regions");

        assert!(are_centers_in_region(&network, 0, None));
        assert!(are_centers_in_region(&network, 1, Some("center")));
        assert!(are_centers_in_region(&network, 2, Some("center_and_east")));
    }

    #[test]
    fn test_region_of_two_cells_gets_centers_in_both_cells() {
        let network = run("family_regions");
        let fractures = family_fractures(&network, 2);

        // The center cell is from -3 to 3, the cell at +x of it from 3 to 5
        assert!(fractures.iter().any(|poly| poly.translation.x < 3.));
        assert!(fractures.iter().any(|poly| poly.translation.x > 3.));
    }
}

mod user_fractures {
    use super::*;

    #[test]
    fn test_user_fractures_keep_their_vertices() {
        let network = run("user_fractures");
        let vertex_counts: Vec<usize> = network
            .dfn
            .accepted_polys
            .iter()
            .map(|poly| poly.vertices.len())
            .collect();

        assert_eq!(vertex_counts, [3, 4, 5, 8, 12]);
    }

    #[test]
    fn test_user_fracture_outside_the_domain_is_rejected() {
        let network = run("user_fractures");

        assert_eq!(network.dfn.stats.rejected_user_fractures, [6]);
    }
}

mod user_fracture_regions {
    use super::*;

    #[test]
    fn test_user_fracture_outside_its_regions_is_rejected() {
        let network = run("user_fracture_regions");

        assert_eq!(network.dfn.stats.rejected_user_fractures, [2]);
    }
}

mod user_fracture_pieces {
    use super::*;

    #[test]
    fn test_piece_outside_the_domain_is_accepted() {
        let network = run("user_fracture_pieces");

        assert_eq!(network.dfn.accepted_polys.len(), 3);
    }

    #[test]
    fn test_pieces_of_one_fracture_are_one_cluster() {
        let network = run("user_fracture_pieces");
        let polys = &network.dfn.accepted_polys;

        assert_eq!(polys[0].group, polys[1].group);
    }

    #[test]
    fn test_pieces_of_one_fracture_have_no_intersection() {
        let network = run("user_fracture_pieces");

        assert!(network.dfn.intersections.is_empty());
    }

    #[test]
    fn test_fracture_of_two_pieces_is_not_removed_as_isolated() {
        let network = run("user_fracture_pieces");
        let mut final_fractures = network.final_fractures.clone();
        final_fractures.sort_unstable();

        assert_eq!(final_fractures, [0, 1]);
    }
}

mod user_fractures_without_fram {
    use super::*;

    #[test]
    fn test_user_fractures_are_accepted_without_fram() {
        let network = run("user_fractures_without_fram");

        assert!(network.dfn.stats.rejected_user_fractures.is_empty());
        assert_eq!(network.dfn.intersections.len(), 2);
    }

    #[test]
    fn test_fram_rejects_the_same_user_fracture_by_default() {
        let network = run_with_replaced_text(
            "user_fractures_without_fram",
            "check_user_fractures = false",
            "check_user_fractures = true",
        );

        assert_eq!(network.dfn.stats.rejected_user_fractures, [2]);
    }
}

mod fram_rejections {
    use super::*;

    #[test]
    fn test_fram_rejects_fractures_of_a_dense_network() {
        let network = run("fram_rejections");

        assert!(fram_reject_count(&network.dfn.stats.rejection_reasons) > 0);
    }

    #[test]
    fn test_each_rejected_fracture_is_kept() {
        let network = run("fram_rejections");

        assert!(network.dfn.stats.rejected_count > 0);
        assert_eq!(
            network.dfn.rejected_polys.len(),
            network.dfn.stats.rejected_count
        );
    }

    #[test]
    fn test_rejected_fractures_file_is_written() {
        let network = run("fram_rejections");

        assert!(output_file_names(&network).contains(&"fractures_rejected.vtp".to_string()));
    }

    #[test]
    fn test_no_intersection_is_shorter_than_min_feature_size() {
        let network = run("fram_rejections");
        let min_feature_size = network.input.min_feature_size;

        for intersection in &network.dfn.intersections {
            assert!((intersection.p1 - intersection.p2).norm() >= min_feature_size);
        }
    }
}

mod fram_disabled {
    use super::*;

    #[test]
    fn test_fram_rejects_no_fracture_when_disabled() {
        let network = run("fram_disabled");

        assert_eq!(fram_reject_count(&network.dfn.stats.rejection_reasons), 0);
    }
}

mod fram_relaxed {
    use super::*;

    #[test]
    fn test_relaxed_fram_rejects_only_short_intersections() {
        let network = run("fram_relaxed");
        let reasons = &network.dfn.stats.rejection_reasons;

        assert_eq!(fram_reject_count(reasons), reasons.short_intersection);
    }
}

mod no_triple_intersections {
    use super::*;

    #[test]
    fn test_network_has_no_triple_point() {
        let network = run("no_triple_intersections");

        assert!(network.dfn.triple_points.is_empty());
        assert!(!output_file_names(&network).contains(&"triple_points.dat".to_string()));
    }
}

mod clusters {
    use super::*;

    #[test]
    fn test_keep_isolated_fractures_keeps_all_accepted_fractures() {
        let network = run("clusters_keep_isolated");
        let all_fractures: Vec<usize> = (0..network.dfn.accepted_polys.len()).collect();

        assert_eq!(network.final_fractures, all_fractures);
        assert!(!output_file_names(&network).contains(&"fractures_removed.vtp".to_string()));
    }

    #[test]
    fn test_keep_only_largest_keeps_the_cluster_with_the_most_fractures() {
        let network = run("clusters_keep_largest");
        let polys = &network.dfn.accepted_polys;
        let fractures_in_group = |group| polys.iter().filter(|poly| poly.group == group).count();
        let largest_size = polys
            .iter()
            .map(|poly| fractures_in_group(poly.group))
            .max()
            .unwrap();

        let final_group = polys[network.final_fractures[0]].group;
        assert!(network
            .final_fractures
            .iter()
            .all(|&index| polys[index].group == final_group));
        assert_eq!(network.final_fractures.len(), largest_size);
    }

    /// Groups of `keep_only_connected_to` of the cases: west (1) and east (2)
    const WEST_AND_EAST: &str = r#"keep_only_connected_to = [["west"], ["east"]]"#;

    /// Each boundary surface in its own group
    fn each_surface(network: &Network) -> Vec<Vec<usize>> {
        let surface_count = network.input.domain_mesh.boundary_surfaces().len();
        (1..=surface_count)
            .map(|surface_id| vec![surface_id])
            .collect()
    }

    /// True if the segments touch a surface of each group
    fn touch_each_group(segments: &[&BoundarySegment], groups: &[Vec<usize>]) -> bool {
        groups.iter().all(|group| {
            segments.iter().any(|segment| {
                group
                    .iter()
                    .any(|surface_id| segment.side.surface_ids.contains(surface_id))
            })
        })
    }

    /// The clusters of the final fractures which do not touch a surface of each group
    fn clusters_not_connected_to(
        network: &Network,
        surface_groups: &[Vec<usize>],
    ) -> Vec<Option<usize>> {
        let polys = final_polys(network);
        let mut groups: Vec<Option<usize>> = polys.iter().map(|poly| poly.group).collect();
        groups.sort_unstable();
        groups.dedup();
        groups.retain(|&group| {
            let cluster: Vec<&Poly> = polys
                .iter()
                .copied()
                .filter(|poly| poly.group == group)
                .collect();
            !touch_each_group(&segments_of(&cluster), surface_groups)
        });
        groups
    }

    #[test]
    fn test_keep_only_connected_to_keeps_the_clusters_on_west_and_east() {
        let network = run("clusters_connect_surfaces");

        assert_eq!(
            clusters_not_connected_to(&network, &each_surface(&network)),
            []
        );
    }

    /// The case `clusters_connect_surfaces` with west and east in one group
    fn run_with_one_group() -> Network {
        run_with_replaced_text(
            "clusters_connect_surfaces",
            WEST_AND_EAST,
            r#"keep_only_connected_to = [["west", "east"]]"#,
        )
    }

    #[test]
    fn test_one_group_keeps_the_clusters_on_one_of_its_surfaces() {
        let network = run_with_one_group();

        assert_eq!(clusters_not_connected_to(&network, &[vec![1, 2]]), []);
    }

    #[test]
    fn test_one_group_keeps_clusters_which_do_not_touch_each_surface() {
        let network = run_with_one_group();

        assert_ne!(
            clusters_not_connected_to(&network, &each_surface(&network)),
            []
        );
    }

    /// The case `clusters_connect_surfaces_mesh`, with the option on or off
    fn run_with_mesh(is_on: bool) -> Network {
        let groups = if is_on {
            WEST_AND_EAST
        } else {
            "keep_only_connected_to = []"
        };
        run_with_replaced_text("clusters_connect_surfaces_mesh", WEST_AND_EAST, groups)
    }

    /// The number of connected parts of the mesh which do not have a point
    /// on each boundary surface
    fn mesh_parts_not_on_all_surfaces(network: &Network) -> usize {
        let mesh = network.mesh.as_ref().unwrap();
        let mut parent: Vec<usize> = (0..mesh.points.len()).collect();
        fn root(parent: &[usize], mut point: usize) -> usize {
            while parent[point] != point {
                point = parent[point];
            }
            point
        }
        for triangle in &mesh.triangles {
            for pair in [[triangle[0], triangle[1]], [triangle[1], triangle[2]]] {
                let (first, second) = (root(&parent, pair[0]), root(&parent, pair[1]));
                parent[first] = second;
            }
        }

        let surface_count = network.input.domain_mesh.boundary_surfaces().len();
        let mut parts: Vec<usize> = mesh
            .triangles
            .iter()
            .map(|triangle| root(&parent, triangle[0]))
            .collect();
        parts.sort_unstable();
        parts.dedup();
        parts
            .into_iter()
            .filter(|&part| {
                (1..=surface_count).any(|surface_id| {
                    !mesh.point_sides.iter().enumerate().any(|(point, sides)| {
                        root(&parent, point) == part
                            && sides
                                .iter()
                                .any(|side| side.surface_ids.contains(&surface_id))
                    })
                })
            })
            .count()
    }

    #[test]
    fn test_each_part_of_the_mesh_touches_all_surfaces() {
        let network = run_with_mesh(true);

        assert_eq!(mesh_parts_not_on_all_surfaces(&network), 0);
    }

    #[test]
    fn test_without_the_option_the_mesh_has_parts_not_on_all_surfaces() {
        let network = run_with_mesh(false);

        assert!(mesh_parts_not_on_all_surfaces(&network) > 0);
    }

    #[test]
    fn test_each_final_fracture_has_triangles_in_the_mesh() {
        let network = run_with_mesh(true);
        let mut fractures = network.mesh.as_ref().unwrap().triangle_fractures.clone();
        fractures.sort_unstable();
        fractures.dedup();

        assert_eq!(fractures.len(), network.final_fractures.len());
    }

    #[test]
    fn test_a_fracture_connected_only_outside_the_domain_is_not_final() {
        let network = run_with_mesh(true);

        // Fracture 2 meets fracture 1 only outside the box
        assert_eq!(network.final_fractures, [0]);
    }

    #[test]
    fn test_a_connecting_cluster_does_not_set_the_flag() {
        let network = run("clusters_connect_surfaces");

        assert!(!network.is_without_connecting_cluster);
    }

    #[test]
    fn test_with_no_connecting_cluster_all_clusters_stay_with_the_flag() {
        // With 100 fractures, no cluster connects the faces west and east
        let network = run_with_replaced_text(
            "clusters_connect_surfaces",
            "\nn_poly = 300",
            "\nn_poly = 100",
        );
        let unchecked = run_with_replacements(
            "clusters_connect_surfaces",
            &[
                ("\nn_poly = 300", "\nn_poly = 100"),
                (WEST_AND_EAST, "keep_only_connected_to = []"),
            ],
        );

        assert!(network.is_without_connecting_cluster);
        assert_eq!(network.final_fractures, unchecked.final_fractures);
    }

    #[test]
    fn test_without_keep_only_connected_to_other_clusters_stay() {
        let network = run_with_replaced_text(
            "clusters_connect_surfaces",
            WEST_AND_EAST,
            "keep_only_connected_to = []",
        );

        assert_ne!(
            clusters_not_connected_to(&network, &each_surface(&network)),
            []
        );
    }
}

mod generation_options {
    use super::*;

    #[test]
    fn test_small_fractures_are_removed() {
        let network = run("remove_small_fractures");
        let StopCondition::NPoly(n_poly) = network.input.stop_condition else {
            panic!("the case must set generation.n_poly");
        };

        assert!(network
            .dfn
            .accepted_polys
            .iter()
            .all(|poly| poly.x_radius >= 1.5));
        assert!(network.dfn.accepted_polys.len() < n_poly);
    }

    #[test]
    fn test_each_family_gets_a_fracture_with_its_maximum_radius() {
        let network = run("force_large_fractures");

        for family in 0..2 {
            assert!(family_fractures(&network, family)
                .iter()
                .any(|poly| poly.x_radius == 4.0));
        }
    }
}

fn final_polys(network: &Network) -> Vec<&Poly> {
    network
        .final_fractures
        .iter()
        .map(|&i| &network.dfn.accepted_polys[i])
        .collect()
}

fn segments_of<'a>(polys: &[&'a Poly]) -> Vec<&'a BoundarySegment> {
    polys
        .iter()
        .flat_map(|poly| &poly.boundary_lines)
        .flat_map(|line| &line.segments)
        .collect()
}

fn other_sides<'a>(polys: &[&'a Poly]) -> Vec<&'a OtherSide> {
    segments_of(polys)
        .into_iter()
        .map(|segment| &segment.side.other_side)
        .collect()
}

/// True if the point is on the surface of the box from -half to half
fn is_on_box_surface(point: &Point3<f64>, half: [f64; 3]) -> bool {
    const TOLERANCE: f64 = 1e-9;

    (0..3).all(|axis| point[axis].abs() <= half[axis] + TOLERANCE)
        && (0..3).any(|axis| (point[axis].abs() - half[axis]).abs() < TOLERANCE)
}

/// Text of an output file which `write_output` writes for the network
fn output_file_text(network: &Network, file_name: &str) -> String {
    let text = Rc::new(RefCell::new(Vec::new()));
    write_output(network, &mut |name| {
        if name == file_name {
            Ok(Box::new(SharedBuffer(Rc::clone(&text))) as Box<dyn Write>)
        } else {
            Ok(Box::new(io::sink()) as Box<dyn Write>)
        }
    })
    .unwrap();

    String::from_utf8(text.take()).unwrap()
}

/// A writer into a buffer which the test keeps (the writer of `write_output`
/// must own its target)
struct SharedBuffer(Rc<RefCell<Vec<u8>>>);

impl Write for SharedBuffer {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        self.0.borrow_mut().extend_from_slice(bytes);
        Ok(bytes.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

mod region_boundaries_mesh {
    use super::*;

    #[test]
    fn test_lines_are_on_the_surface_of_the_mesh() {
        let network = run("region_boundaries_mesh");
        let segments = segments_of(&final_polys(&network));

        // The box mesh goes from -5 to 5
        assert!(!segments.is_empty());
        assert!(segments.iter().all(|segment| {
            segment.side.other_side == OtherSide::MeshSurface
                && is_on_box_surface(&segment.start, [5.; 3])
                && is_on_box_surface(&segment.end, [5.; 3])
        }));
    }

    #[test]
    fn test_fram_rejects_fractures_for_their_boundary_lines() {
        let network = run("region_boundaries_mesh");

        assert!(network.dfn.stats.rejection_reasons.boundary_line > 0);
    }

    #[test]
    fn test_without_fram_the_option_adds_only_the_lines() {
        let fram_off = ("[fram]\nenable = true", "[fram]\nenable = false");
        let with_lines = run_with_replacements("region_boundaries_mesh", &[fram_off]);
        let without_lines = run_with_replacements(
            "region_boundaries_mesh",
            &[
                fram_off,
                (
                    "[region_boundaries]\nenable = true",
                    "[region_boundaries]\nenable = false",
                ),
            ],
        );
        let vertices = |network: &Network| -> Vec<Vec<Point3<f64>>> {
            final_polys(network)
                .iter()
                .map(|poly| poly.vertices.clone())
                .collect()
        };

        assert!(!segments_of(&final_polys(&with_lines)).is_empty());
        assert_eq!(vertices(&with_lines), vertices(&without_lines));
    }

    #[test]
    fn test_fracture_boundaries_file_has_a_line_for_each_segment() {
        let network = run("region_boundaries_mesh");
        let segment_count = segments_of(&final_polys(&network)).len();

        let text = output_file_text(&network, "fracture_boundaries.vtp");

        assert!(text.contains(&format!(r#"NumberOfLines="{segment_count}""#)));
    }
}

mod region_boundaries_family_region {
    use super::*;

    #[test]
    fn test_lines_of_a_family_with_a_region_are_on_the_mesh_surface() {
        let network = run("region_boundaries_family_region");
        let sides = other_sides(&family_fractures(&network, 1));

        assert!(!sides.is_empty());
        assert!(sides.iter().all(|side| **side == OtherSide::MeshSurface));
    }
}

mod region_boundaries_user_fractures {
    use super::*;

    /// The region `north` is region 2 of `mesh_domain.vtu`
    const NORTH: usize = 2;

    fn sides_of_fracture(network: &Network, index: usize) -> Vec<&OtherSide> {
        other_sides(&[&network.dfn.accepted_polys[index]])
    }

    #[test]
    fn test_fracture_in_two_regions_has_no_line_between_them() {
        let network = run("region_boundaries_user_fractures");

        // Only the line at the corner where the mesh has no cell
        assert_eq!(
            sides_of_fracture(&network, 0),
            [&OtherSide::MeshSurface, &OtherSide::MeshSurface]
        );
    }

    #[test]
    fn test_fracture_in_one_region_has_the_other_region_beside_its_line() {
        let network = run("region_boundaries_user_fractures");

        assert_eq!(
            sides_of_fracture(&network, 1),
            [&OtherSide::Cell(vec![NORTH])]
        );
    }

    #[test]
    fn test_fracture_through_both_arms_of_the_l_shape_has_two_lines() {
        let network = run("region_boundaries_user_fractures");

        assert_eq!(network.dfn.accepted_polys[2].boundary_lines.len(), 2);
    }
}

/// The part of the mesh of each point: parts which are not connected have
/// different numbers
fn mesh_parts(mesh: &FractureMesh) -> Vec<usize> {
    let mut parent: Vec<usize> = (0..mesh.points.len()).collect();
    fn root(parent: &mut [usize], mut point: usize) -> usize {
        while parent[point] != point {
            parent[point] = parent[parent[point]];
            point = parent[point];
        }
        point
    }

    for triangle in &mesh.triangles {
        for pair in [[triangle[0], triangle[1]], [triangle[1], triangle[2]]] {
            let (first, second) = (root(&mut parent, pair[0]), root(&mut parent, pair[1]));
            parent[first] = second;
        }
    }

    (0..mesh.points.len())
        .map(|point| root(&mut parent, point))
        .collect()
}

/// Number of parts of the mesh which are not connected
fn mesh_part_count(mesh: &FractureMesh) -> usize {
    let mut parts = mesh_parts(mesh);
    parts.sort_unstable();
    parts.dedup();
    parts.len()
}

/// True if each vertex of the final fracture is in the box of the cases
/// (−5 to 5)
fn is_in_box(network: &Network, position: usize) -> bool {
    network.dfn.accepted_polys[network.final_fractures[position]]
        .vertices
        .iter()
        .all(|vertex| vertex.coords.amax() <= 5.)
}

/// Positions of the final fractures which are in the mesh
fn fractures_in_mesh(mesh: &FractureMesh) -> Vec<usize> {
    let mut fractures = mesh.triangle_fractures.clone();
    fractures.sort_unstable();
    fractures.dedup();
    fractures
}

/// The smallest angle of each triangle (degrees)
fn smallest_angles(mesh: &FractureMesh) -> Vec<f64> {
    mesh.triangles
        .iter()
        .map(|triangle| {
            let corners = triangle.map(|point| mesh.points[point]);
            (0..3)
                .map(|i| {
                    let to_next = corners[(i + 1) % 3] - corners[i];
                    let to_previous = corners[(i + 2) % 3] - corners[i];
                    to_next.angle(&to_previous).to_degrees()
                })
                .fold(180., f64::min)
        })
        .collect()
}

/// The length of each edge of each triangle
fn edge_lengths(mesh: &FractureMesh) -> impl Iterator<Item = f64> + '_ {
    mesh.triangles.iter().flat_map(move |triangle| {
        (0..3).map(move |i| (mesh.points[triangle[(i + 1) % 3]] - mesh.points[triangle[i]]).norm())
    })
}

fn triangle_area(mesh: &FractureMesh, triangle: &[usize; 3]) -> f64 {
    let [a, b, c] = triangle.map(|point| mesh.points[point]);
    (b - a).cross(&(c - a)).norm() / 2.
}

mod mesh {
    use super::*;

    #[test]
    fn test_mesh_file_is_written() {
        let network = run("mesh");

        assert!(output_file_names(&network).contains(&"dfn_mesh.vtu".to_string()));
    }

    #[test]
    fn test_each_triangle_is_in_the_box() {
        let network = run("mesh");
        let mesh = network.mesh.as_ref().unwrap();

        for triangle in &mesh.triangles {
            let [a, b, c] = triangle.map(|point| mesh.points[point]);
            let centroid = (a.coords + b.coords + c.coords) / 3.;
            assert!(centroid.amax() < 5., "{centroid}");
        }
    }

    #[test]
    fn test_mesh_has_less_area_than_a_fracture_which_crosses_the_box_surface() {
        let network = run("mesh");
        let mesh = network.mesh.as_ref().unwrap();
        let crossing = fractures_in_mesh(mesh)
            .into_iter()
            .find(|&position| !is_in_box(&network, position))
            .expect("a fracture of the case crosses the box surface");

        let mesh_area: f64 = mesh
            .triangles
            .iter()
            .zip(&mesh.triangle_fractures)
            .filter(|(_, &fracture)| fracture == crossing)
            .map(|(triangle, _)| triangle_area(mesh, triangle))
            .sum();

        let area = network.dfn.accepted_polys[network.final_fractures[crossing]].area;
        assert!(
            mesh_area < 0.999 * area,
            "{mesh_area} is not less than {area}"
        );
    }

    #[test]
    fn test_each_part_of_the_mesh_has_a_boundary_point() {
        let network = run("mesh");
        let mesh = network.mesh.as_ref().unwrap();
        let parts = mesh_parts(mesh);
        let mut parts_with_boundary: Vec<usize> = parts
            .iter()
            .zip(&mesh.point_sides)
            .filter(|(_, sides)| !sides.is_empty())
            .map(|(&part, _)| part)
            .collect();
        parts_with_boundary.sort_unstable();
        parts_with_boundary.dedup();

        assert_eq!(parts_with_boundary.len(), mesh_part_count(mesh));
    }

    #[test]
    fn test_clusters_without_a_boundary_point_are_not_in_the_mesh() {
        let network = run("mesh");

        assert!(
            fractures_in_mesh(network.mesh.as_ref().unwrap()).len() < network.final_fractures.len()
        );
    }

    #[test]
    fn test_triangles_cover_each_fracture_in_the_box() {
        let network = run("mesh");
        let mesh = network.mesh.as_ref().unwrap();
        let mut mesh_areas = vec![0.; network.final_fractures.len()];
        for (triangle, &fracture) in mesh.triangles.iter().zip(&mesh.triangle_fractures) {
            mesh_areas[fracture] += triangle_area(mesh, triangle);
        }

        for position in fractures_in_mesh(mesh)
            .into_iter()
            .filter(|&position| is_in_box(&network, position))
        {
            let area = network.dfn.accepted_polys[network.final_fractures[position]].area;
            assert!(
                (mesh_areas[position] - area).abs() < 1e-9 * area,
                "fracture {}: mesh area {} is not {area}",
                position + 1,
                mesh_areas[position]
            );
        }
    }

    #[test]
    fn test_almost_no_triangle_has_an_angle_below_20_degrees() {
        let network = run("mesh");
        let angles = smallest_angles(network.mesh.as_ref().unwrap());

        let small_count = angles.iter().filter(|&&angle| angle < 20.).count();
        assert!(
            small_count * 1000 < angles.len(),
            "{small_count} of {} triangles",
            angles.len()
        );
    }

    #[test]
    fn test_aperture_of_family_2_is_correlated_with_the_radius() {
        let network = run("mesh");
        let mesh = network.mesh.as_ref().unwrap();

        for (position, &i) in network.final_fractures.iter().enumerate() {
            let poly = &network.dfn.accepted_polys[i];
            if poly.origin == FractureOrigin::Family(1) {
                let expected = 1e-4 * poly.x_radius.sqrt();
                assert!((mesh.apertures[position] - expected).abs() < 1e-15 * expected);
            }
        }
    }

    #[test]
    fn test_mesh_surface_points_are_on_the_box_surface() {
        let network = run("mesh");
        let mesh = network.mesh.as_ref().unwrap();
        let surface_points: Vec<&Point3<f64>> = mesh
            .points
            .iter()
            .zip(&mesh.point_sides)
            .filter(|(_, sides)| {
                sides
                    .iter()
                    .any(|side| side.other_side == OtherSide::MeshSurface)
            })
            .map(|(point, _)| point)
            .collect();

        assert!(!surface_points.is_empty());
        for point in surface_points {
            let largest_coordinate = point.coords.amax();
            assert!((largest_coordinate - 5.).abs() < 1e-9, "{point}");
        }
    }

    #[test]
    fn test_mesh_does_not_change_the_network() {
        let with_mesh = run("mesh");
        let without_mesh =
            run_with_replaced_text("mesh", "[mesh]\nenable = true", "[mesh]\nenable = false");
        let vertices = |network: &Network| -> Vec<Vec<Point3<f64>>> {
            network
                .dfn
                .accepted_polys
                .iter()
                .map(|poly| poly.vertices.clone())
                .collect()
        };

        assert_eq!(vertices(&with_mesh), vertices(&without_mesh));
    }
}

mod mesh_user_fractures {
    use super::*;

    #[test]
    fn test_pieces_and_the_fracture_on_their_common_edge_are_one_part() {
        let network = run("mesh_user_fractures");

        assert_eq!(mesh_part_count(network.mesh.as_ref().unwrap()), 1);
    }

    #[test]
    fn test_user_fractures_have_the_aperture_of_their_file() {
        let network = run("mesh_user_fractures");

        assert!(network
            .mesh
            .as_ref()
            .unwrap()
            .apertures
            .iter()
            .all(|&aperture| aperture == 2e-4));
    }
}

mod mesh_user_fractures_tiny {
    use super::*;

    #[test]
    fn test_tiny_features_give_no_edge_shorter_than_a_quarter_of_h() {
        let network = run("mesh_user_fractures_tiny");

        let shortest = edge_lengths(network.mesh.as_ref().unwrap()).fold(f64::INFINITY, f64::min);
        assert!(shortest >= 0.025, "{shortest}");
    }

    #[test]
    fn test_tiny_features_give_no_angle_below_20_degrees() {
        let network = run("mesh_user_fractures_tiny");

        let smallest = smallest_angles(network.mesh.as_ref().unwrap())
            .into_iter()
            .fold(180., f64::min);
        assert!(smallest > 19.9, "{smallest}");
    }

    #[test]
    fn test_fractures_with_tiny_features_are_one_part() {
        let network = run("mesh_user_fractures_tiny");

        assert_eq!(mesh_part_count(network.mesh.as_ref().unwrap()), 1);
    }
}

mod mesh_regions {
    use super::*;

    /// Region id (position in the regions of the domain mesh plus 1) of the name
    fn region_id(network: &Network, name: &str) -> usize {
        network
            .input
            .domain_mesh
            .regions()
            .iter()
            .position(|region| region.name == name)
            .expect("the mesh has the region")
            + 1
    }

    /// The regions of the cell at the centroid of each triangle of the
    /// fractures which the filter selects
    fn triangle_regions(network: &Network, select: impl Fn(&Poly) -> bool) -> Vec<Vec<usize>> {
        let mesh = network.mesh.as_ref().unwrap();
        mesh.triangles
            .iter()
            .zip(&mesh.triangle_fractures)
            .filter(|(_, &position)| {
                select(&network.dfn.accepted_polys[network.final_fractures[position]])
            })
            .map(|(triangle, _)| {
                let [a, b, c] = triangle.map(|point| mesh.points[point]);
                let centroid = Point3::from((a.coords + b.coords + c.coords) / 3.);
                network
                    .input
                    .domain_mesh
                    .regions_at(&centroid)
                    .expect("each triangle is in the domain mesh")
                    .to_vec()
            })
            .collect()
    }

    #[test]
    fn test_family_triangles_cross_the_family_region_boundary() {
        let network = run("mesh_regions");
        let north = region_id(&network, "north");

        let regions = triangle_regions(&network, |poly| poly.origin == FractureOrigin::Family(0));

        assert!(regions
            .iter()
            .any(|cell_regions| cell_regions.contains(&north)));
    }

    #[test]
    fn test_user_fracture_triangles_are_in_its_region() {
        let network = run("mesh_regions");
        let south = region_id(&network, "south");

        let regions = triangle_regions(&network, |poly| {
            poly.origin == FractureOrigin::UserDefined && poly.region_ids == [south]
        });

        assert!(!regions.is_empty());
        assert!(regions
            .iter()
            .all(|cell_regions| cell_regions.contains(&south)));
    }

    #[test]
    fn test_user_fracture_without_region_has_triangles_in_both_regions() {
        let network = run("mesh_regions");
        let (south, north) = (region_id(&network, "south"), region_id(&network, "north"));

        let regions = triangle_regions(&network, |poly| {
            poly.origin == FractureOrigin::UserDefined && poly.region_ids.is_empty()
        });

        assert!(regions
            .iter()
            .any(|cell_regions| cell_regions.contains(&south)));
        assert!(regions
            .iter()
            .any(|cell_regions| cell_regions.contains(&north)));
    }
}

mod boundary_surfaces {
    use super::*;

    const TOLERANCE: f64 = 1e-9;

    /// Boundary N of the boundary surface with the name
    fn surface_id(network: &Network, name: &str) -> usize {
        network
            .input
            .domain_mesh
            .boundary_surfaces()
            .iter()
            .position(|surface| surface.name == name)
            .expect("the mesh has the boundary surface")
            + 1
    }

    fn is_at_x(point: &Point3<f64>, x: f64) -> bool {
        (point.x - x).abs() < TOLERANCE
    }

    /// The segments of the final fractures with both ends at x
    fn segments_at_x(network: &Network, x: f64) -> Vec<&BoundarySegment> {
        segments_of(&final_polys(network))
            .into_iter()
            .filter(|segment| is_at_x(&segment.start, x) && is_at_x(&segment.end, x))
            .collect()
    }

    #[test]
    fn test_segments_on_the_west_face_are_on_the_boundary_west() {
        let network = run("boundary_surfaces");
        let west = surface_id(&network, "west");
        let segments = segments_at_x(&network, -5.);

        assert!(!segments.is_empty());
        assert!(segments
            .iter()
            .all(|segment| segment.side.surface_ids == [west]));
    }

    #[test]
    fn test_segments_on_the_east_face_are_on_the_boundary_east() {
        let network = run("boundary_surfaces");
        let east = surface_id(&network, "east");
        let segments = segments_at_x(&network, 5.);

        assert!(!segments.is_empty());
        assert!(segments
            .iter()
            .all(|segment| segment.side.surface_ids == [east]));
    }

    #[test]
    fn test_segments_on_the_other_faces_are_on_no_boundary_surface() {
        let network = run("boundary_surfaces");
        let segments: Vec<&BoundarySegment> = segments_of(&final_polys(&network))
            .into_iter()
            .filter(|segment| !is_at_x(&segment.start, -5.) && !is_at_x(&segment.start, 5.))
            .collect();

        assert!(!segments.is_empty());
        assert!(segments
            .iter()
            .all(|segment| segment.side.surface_ids.is_empty()));
    }

    #[test]
    fn test_box_mesh_without_surface_cells_gives_no_boundary_surface() {
        let network = run_with_replaced_text(
            "boundary_surfaces",
            r#"domain = "boundary_surfaces.vtu""#,
            r#"domain = "../box.vtu""#,
        );

        assert!(segments_of(&final_polys(&network))
            .iter()
            .all(|segment| segment.side.surface_ids.is_empty()));
    }

    #[test]
    fn test_mesh_points_on_a_line_of_the_west_face_are_on_the_boundary_west() {
        let network = run("boundary_surfaces");
        let west = surface_id(&network, "west");
        let mesh = network.mesh.as_ref().unwrap();
        let west_sides: Vec<&Vec<SegmentSide>> = mesh
            .points
            .iter()
            .zip(&mesh.point_sides)
            .filter(|(point, sides)| is_at_x(point, -5.) && !sides.is_empty())
            .map(|(_, sides)| sides)
            .collect();

        assert!(!west_sides.is_empty());
        assert!(west_sides
            .iter()
            .all(|sides| sides.iter().any(|side| side.surface_ids == [west])));
    }

    /// The values of the data array with the name in the text of a VTK file
    fn array_values(text: &str, name: &str) -> Vec<u8> {
        let start = text
            .find(&format!(r#"Name="{name}""#))
            .unwrap_or_else(|| panic!("the file has no array {name}"));
        let values_start = start + text[start..].find('>').unwrap() + 1;
        let values_end = values_start + text[values_start..].find('<').unwrap();

        text[values_start..values_end]
            .split_whitespace()
            .map(|value| value.parse().unwrap())
            .collect()
    }

    #[test]
    fn test_output_files_have_an_array_for_each_boundary_surface() {
        let network = run("boundary_surfaces");

        for file_name in ["fracture_boundaries.vtp", "dfn_mesh.vtu"] {
            let text = output_file_text(&network, file_name);
            for name in ["boundary__west", "boundary__east"] {
                assert!(
                    array_values(&text, name).contains(&1),
                    "{file_name}: {name} has no 1"
                );
            }
        }
    }
}

mod checkpoint {
    use dfngen::{resume, Checkpoint, DfngenError, Resume, ResumeMode};

    use super::*;

    /// Number of accepted stochastic fractures of the network
    fn stochastic_count(network: &Network) -> usize {
        network
            .dfn
            .accepted_polys
            .iter()
            .filter(|poly| poly.origin != FractureOrigin::UserDefined)
            .count()
    }

    /// The checkpoint of a run of a case after the accepted stochastic
    /// fracture `accepted_count`, through its bytes
    fn checkpoint_at(case: &str, accepted_count: usize) -> Checkpoint {
        let (input, families) = input_with_replacements(case, &[]);
        let mut bytes = None;
        generate(input, families, &mut |update| {
            if update.accepted_count == accepted_count {
                bytes = Some(update.checkpoint().to_bytes().unwrap());
            }
        })
        .unwrap();

        Checkpoint::from_bytes(&bytes.expect("no update with this count")).unwrap()
    }

    fn run_from(case: &str, checkpoint: &Checkpoint, mode: ResumeMode) -> Network {
        let (input, families) = input_with_replacements(case, &[]);
        let restart = Resume { checkpoint, mode };

        resume(input, families, restart, &mut |_| {}).unwrap()
    }

    /// A run with a restart from the checkpoint at half of the stochastic
    /// fractures gives the output files of a run with no stop
    fn assert_restart_gives_the_output_of_a_run_with_no_stop(case: &str) {
        let network = run(case);
        let checkpoint = checkpoint_at(case, stochastic_count(&network) / 2);
        let restarted = run_from(case, &checkpoint, ResumeMode::ContinueInsertion);

        let names = output_file_names(&network);
        assert_eq!(output_file_names(&restarted), names);
        for name in names {
            assert!(
                output_file_text(&restarted, &name) == output_file_text(&network, &name),
                "{case}: {name} differs"
            );
        }
    }

    #[test]
    fn test_restart_gives_the_network_of_a_run_with_no_stop() {
        assert_restart_gives_the_output_of_a_run_with_no_stop("base");
    }

    #[test]
    fn test_restart_with_the_p32_stop_condition_gives_the_network_of_a_run_with_no_stop() {
        assert_restart_gives_the_output_of_a_run_with_no_stop("stop_p32");
    }

    #[test]
    fn test_restart_keeps_the_rejected_fractures() {
        assert_restart_gives_the_output_of_a_run_with_no_stop("fram_rejections");
    }

    #[test]
    fn test_restart_gives_the_mesh_of_a_run_with_no_stop() {
        assert_restart_gives_the_output_of_a_run_with_no_stop("mesh");
    }

    #[test]
    fn test_restart_with_user_fractures_and_regions_gives_the_network_of_a_run_with_no_stop() {
        assert_restart_gives_the_output_of_a_run_with_no_stop("mesh_regions");
    }

    #[test]
    fn test_finish_keeps_only_the_fractures_of_the_checkpoint() {
        let network = run_from("base", &checkpoint_at("base", 40), ResumeMode::Finish);

        assert_eq!(stochastic_count(&network), 40);
    }

    #[test]
    fn test_restart_with_another_input_is_an_error() {
        let checkpoint = checkpoint_at("base", 40);
        let (input, families) = input_with_replacements("base", &[("seed = 12345", "seed = 1")]);
        let restart = Resume {
            checkpoint: &checkpoint,
            mode: ResumeMode::ContinueInsertion,
        };

        let result = resume(input, families, restart, &mut |_| {});

        assert!(matches!(result, Err(DfngenError::CheckpointInputChanged)));
    }
}
