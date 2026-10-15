//! One test case of `tests/data` for each feature of the input. Each test
//! checks what its feature must do with the network.

use std::fs;
use std::io::{self, Write};
use std::path::{Path, PathBuf};

use dfngen::fracture::poly::{FractureOrigin, Poly};
use dfngen::io::input::StopCondition;
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
    let folder = case_folder(case);
    let toml_text = fs::read_to_string(folder.join(format!("{case}.toml"))).unwrap();
    assert!(
        toml_text.contains(from),
        "the input of {case} has no {from:?}"
    );
    let toml_text = toml_text.replace(from, to);

    let names = read_input_file_names(&toml_text).unwrap();
    let read_text = |name: &str| fs::read_to_string(folder.join(name)).unwrap();
    let mesh_text = read_text(&names.domain);
    let fractures_text = names.user_fractures.as_deref().map(read_text);
    let (input, families) =
        build_input_from_text(&toml_text, &mesh_text, fractures_text.as_deref()).unwrap();

    generate(input, families, &mut |_| {}).unwrap()
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
