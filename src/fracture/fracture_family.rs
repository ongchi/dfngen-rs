use std::rc::Rc;

use parry3d_f64::na::Vector3;
use rand::Rng;
use rand_distr::Uniform;
use rand_mt::Mt64;

use crate::distribution::{Fisher, RadiusDistribution};
use crate::domain_mesh::DomainMesh;
use crate::error::{DfngenError, FamilyError};

use super::poly::{FractureOrigin, Poly};

/// FractureFamily holds the variables of one stochastic family of the input file.
pub struct FractureFamily {
    pub id: usize,

    /// Number of vertices of each fracture.
    pub number_of_points: usize,

    pub radius: RadiusDistribution,

    /// Index of the next radius to use in `radii_list`
    pub radii_index: usize,

    /// Initial list of fracture/polygon radii, sorted largest to smallest.
    pub radii_list: Vec<f64>,

    /// Region the family belongs to. 0 is entire domain, greater than 0 is
    /// the position of the region in `DomainMesh::regions` plus 1.
    pub region_id: usize,

    /// Aspect ratio for family.
    pub aspect_ratio: f64,

    /// Target p32 (fracture intensity) for the family when using p32 program-stopping option.
    pub p32_target: f64,

    /// Current P32 value for this family.
    pub current_p32: f64,

    /// Rotation, or twist, around the z normal before the 3d rotation, in radians.
    /// `None` gives each fracture a uniform random rotation.
    pub beta: Option<f64>,

    pub orientation: Fisher,

    /// The fracture centers are sampled in the cells of the family's region.
    pub domain_mesh: Rc<DomainMesh>,
}

impl FractureFamily {
    /// Volume of the family's region, or of the whole domain.
    pub fn region_volume(&self, domain_volume: f64, region_volumes: &[f64]) -> f64 {
        if self.region_id > 0 {
            region_volumes[self.region_id - 1]
        } else {
            domain_volume
        }
    }

    /// Adds the area of an accepted polygon to the family's current P32.
    /// The volume is that of the family's region, or of the whole domain.
    pub fn add_area_to_p32(&mut self, area: f64, domain_volume: f64, region_volumes: &[f64]) {
        self.current_p32 += area * 2. / self.region_volume(domain_volume, region_volumes);
    }

    /// Area of a fracture of the family with this radius: the polygon of
    /// `Poly::new` is a regular polygon, stretched in y by the aspect ratio.
    pub fn fracture_area(&self, radius: f64) -> f64 {
        let number_of_points = self.number_of_points as f64;
        let angle_step = 2. * std::f64::consts::PI / number_of_points;

        number_of_points / 2. * angle_step.sin() * self.aspect_ratio * radius * radius
    }

    /// Fills the radii list until the radii have the P32 target of the family.
    ///
    /// P32 counts both sides of a fracture. With `force_large_fractures`,
    /// the list starts with the maximum radius.
    fn generate_p32_radii(
        &mut self,
        force_large_fractures: bool,
        volume: f64,
        generator: &mut Mt64,
    ) -> Result<(), DfngenError> {
        let mut p32 = 0.;

        if force_large_fractures {
            self.radii_list.push(self.radius.max);
            p32 += self.fracture_area(self.radius.max) * 2. / volume;
        }

        while p32 < self.p32_target {
            let radius = self.sample_radius(generator)?;
            self.radii_list.push(radius);
            p32 += self.fracture_area(radius) * 2. / volume;
        }

        Ok(())
    }

    /// Uniform random position in the cells of the family's region
    pub fn sample_position(&self, generator: &mut Mt64) -> Vector3<f64> {
        self.domain_mesh
            .sample_point(self.region_id, generator)
            .coords
    }

    /// Samples a radius from the radius distribution of the family.
    fn sample_radius(&self, generator: &mut Mt64) -> Result<f64, DfngenError> {
        self.radius
            .sample(generator)
            .map_err(|source| DfngenError::FamilyRadius {
                family: self.id,
                source: Box::new(source),
            })
    }

    /// Translates the polygon to a new position in the region of the family.
    /// The size, the normal vector and the rotation of the polygon stay.
    /// Thus a fracture which was rejected gets another attempt, and the
    /// family keeps its target distributions.
    pub fn retranslate(&self, poly: &mut Poly, generator: &mut Mt64) {
        // The cluster group and the intersections are those of the old position
        poly.group = None;
        poly.intersection_indices.clear();

        // Move poly back to origin
        for vertex in &mut poly.vertices {
            *vertex -= poly.translation;
        }

        let position = self.sample_position(generator);
        poly.translate(position);
    }

    pub fn normal_vector(&self, rng: &mut Mt64) -> Vector3<f64> {
        rng.sample(&self.orientation)
    }

    pub fn create_poly(
        &mut self,
        eps: f64,
        family_index: usize,
        generator: &mut Mt64,
    ) -> Result<Poly, DfngenError> {
        let radius = if self.radii_index >= self.radii_list.len() {
            // If out of radii from list, insert random radius
            self.sample_radius(generator)?
        } else {
            // Insert radius from list
            let radius = self.radii_list[self.radii_index];
            self.radii_index += 1;
            radius
        };

        // New polygon to build
        let mut new_poly = Poly::new(self.number_of_points, radius, self.aspect_ratio);
        new_poly.origin = FractureOrigin::Family(family_index);

        let beta = self.beta.unwrap_or_else(|| {
            let uniform_dist = Uniform::new(0., 2. * std::f64::consts::PI)
                .expect("the limits are constants and are valid");
            generator.sample(uniform_dist)
        });

        // Apply 2d rotation matrix, twist around origin
        // assumes polygon on x-y plane
        // Angle must be in rad
        new_poly.rotation_2d(beta);

        // Rotate vertices to norm (new normal)
        let norm = self.normal_vector(generator).normalize();
        new_poly.rotation_3d(&norm, eps);

        // Save the new fracture's new normal vector
        new_poly.normal = norm;

        // The translation is also saved in the polygon
        let position = self.sample_position(generator);
        new_poly.translate(position);

        Ok(new_poly)
    }
}

#[derive(Default)]
pub struct FractureFamilyBuilder {
    id: Option<usize>,
    number_of_points: Option<usize>,
    radius: Option<RadiusDistribution>,
    orientation: Option<Fisher>,
    aspect_ratio: Option<f64>,
    beta: Option<f64>,
    p32_target: Option<f64>,
    region_id: Option<usize>,
    domain_mesh: Option<Rc<DomainMesh>>,
}

impl FractureFamilyBuilder {
    pub fn id(&mut self, id: usize) -> &mut Self {
        self.id = Some(id);
        self
    }

    pub fn number_of_points(&mut self, number_of_points: usize) -> &mut Self {
        self.number_of_points = Some(number_of_points);
        self
    }

    pub fn radius(&mut self, radius: RadiusDistribution) -> &mut Self {
        self.radius = Some(radius);
        self
    }

    pub fn aspect_ratio(&mut self, aspect_ratio: f64) -> &mut Self {
        self.aspect_ratio = Some(aspect_ratio);
        self
    }

    pub fn orientation(&mut self, orientation: Fisher) -> &mut Self {
        self.orientation = Some(orientation);
        self
    }

    pub fn beta(&mut self, beta: f64) -> &mut Self {
        self.beta = Some(beta);
        self
    }

    pub fn p32_target(&mut self, p32_target: f64) -> &mut Self {
        self.p32_target = Some(p32_target);
        self
    }

    pub fn region_id(&mut self, region: usize) -> &mut Self {
        self.region_id = Some(region);
        self
    }

    /// The fracture centers are sampled in the cells of the domain mesh.
    pub fn domain_mesh(&mut self, domain_mesh: Rc<DomainMesh>) -> &mut Self {
        self.domain_mesh = Some(domain_mesh);
        self
    }

    pub fn build(&mut self) -> Result<FractureFamily, FamilyError> {
        Ok(FractureFamily {
            id: self.id.take().ok_or(FamilyError::MissingValue("id"))?,
            number_of_points: self
                .number_of_points
                .ok_or(FamilyError::MissingValue("number of points"))?,
            radius: self
                .radius
                .take()
                .ok_or(FamilyError::MissingValue("radius"))?,
            radii_index: 0,
            radii_list: Vec::new(),
            region_id: self
                .region_id
                .ok_or(FamilyError::MissingValue("region_id"))?,
            aspect_ratio: self
                .aspect_ratio
                .ok_or(FamilyError::MissingValue("aspect ratio"))?,
            p32_target: self.p32_target.unwrap_or(0.),
            current_p32: 0.,
            beta: self.beta,
            orientation: self
                .orientation
                .take()
                .ok_or(FamilyError::MissingValue("orientation"))?,
            domain_mesh: self
                .domain_mesh
                .take()
                .ok_or(FamilyError::MissingValue("domain mesh"))?,
        })
    }
}

pub struct FractureFamilyCollection {
    pub families: Vec<FractureFamily>,

    // Each element is the probability of chosing a fracture from
    // the element's corresponding family to be inserted into the DFN.
    // The elements should add up to 1.0 (for %100).
    pub probabilities: Vec<f64>,
}

impl FractureFamilyCollection {
    /// Generates a list of radii by distribution function of each family.
    pub fn generate_radii(
        &mut self,
        force_large_fractures: bool,
        mut n_poly: usize,
        generator: &mut Mt64,
    ) -> Result<(), DfngenError> {
        if force_large_fractures {
            for family in &mut self.families {
                let radius = family.radius.max;
                family.radii_list.push(radius);
            }
            n_poly -= self.families.len();
        }

        for (family, probability) in self.families.iter_mut().zip(&self.probabilities) {
            for _ in 0..((probability * n_poly as f64).ceil() as usize) {
                let radius = family.sample_radius(generator)?;
                family.radii_list.push(radius);
            }
        }

        Ok(())
    }

    /// Generates the radii list of each family for the P32 stop condition:
    /// sampled radii until their area gives the P32 target of the family.
    ///
    /// # Arguments
    ///
    /// * `force_large_fractures` - If true, each list starts with the maximum radius
    /// * `domain_volume` - Volume of the domain mesh
    /// * `region_volumes` - Volume of each region of the domain mesh
    /// * `generator` - Random number generator
    pub fn generate_p32_radii(
        &mut self,
        force_large_fractures: bool,
        domain_volume: f64,
        region_volumes: &[f64],
        generator: &mut Mt64,
    ) -> Result<(), DfngenError> {
        for family in &mut self.families {
            let volume = family.region_volume(domain_volume, region_volumes);
            family.generate_p32_radii(force_large_fractures, volume, generator)?;
        }

        Ok(())
    }

    /// Sort each family's radii list from largest to smallest.
    /// This will allow the DFN gereration to start from largest to smallest
    /// fractures.
    pub fn sort_radii(&mut self) {
        for family in &mut self.families {
            family.radii_list.sort_by(|a, b| b.total_cmp(a));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::input::test_inputs::{parse, GENERAL_OPTIONS};
    use crate::io::input::Input;

    /// Input without n_poly and with one family of this radius and P32 target
    fn parse_p32_family(radius: &str, p32_target: f64) -> (Input, FractureFamilyCollection) {
        let family = format!(
            r#"
            [[family]]
            number_of_points = 8
            probability = 1.0
            aspect_ratio = 1.5
            p32_target = {p32_target}
            radius = {radius}
            orientation = {{ theta = 0.0, phi = 0.0, kappa = 8.0 }}
            "#
        );
        let options = GENERAL_OPTIONS.replace("n_poly = 60", "");

        parse(&[&options, &family]).unwrap()
    }

    fn list_p32(family: &FractureFamily, radii: &[f64], volume: f64) -> f64 {
        radii
            .iter()
            .map(|&radius| family.fracture_area(radius) * 2. / volume)
            .sum()
    }

    #[test]
    fn test_fracture_area_is_the_area_of_the_polygon() {
        let (_, mut families) =
            parse_p32_family(r#"{ distribution = "constant", value = 2.0 }"#, 0.1);
        let family = &mut families.families[0];

        for number_of_points in [3, 4, 8] {
            for aspect_ratio in [1.0, 2.0] {
                family.number_of_points = number_of_points;
                family.aspect_ratio = aspect_ratio;
                let mut poly = Poly::new(number_of_points, 2.0, aspect_ratio);
                poly.assign_area();

                let area = family.fracture_area(2.0);

                assert!(
                    (area - poly.area).abs() < 1e-12 * poly.area,
                    "{number_of_points} points, aspect ratio {aspect_ratio}: {area} != {}",
                    poly.area
                );
            }
        }
    }

    #[test]
    fn test_p32_radii_list_of_constant_radius_has_the_target_count() {
        let (input, mut families) =
            parse_p32_family(r#"{ distribution = "constant", value = 2.0 }"#, 0.5);
        let family_area = families.families[0].fracture_area(2.0);
        let expected = (0.5 * input.domain_volume / (2. * family_area)).ceil() as usize;

        families
            .generate_p32_radii(
                false,
                input.domain_volume,
                &input.region_volumes,
                &mut Mt64::new(1),
            )
            .unwrap();

        assert_eq!(families.families[0].radii_list.len(), expected);
    }

    #[test]
    fn test_p32_radii_list_stops_at_the_first_radius_which_reaches_the_target() {
        let (input, mut families) = parse_p32_family(
            r#"{ distribution = "truncated_power_law", alpha = 1.5, min = 1.0, max = 5.0 }"#,
            0.5,
        );

        families
            .generate_p32_radii(
                false,
                input.domain_volume,
                &input.region_volumes,
                &mut Mt64::new(1),
            )
            .unwrap();

        let family = &families.families[0];
        let radii = &family.radii_list;
        let all_but_last = &radii[..radii.len() - 1];
        assert!(
            list_p32(family, radii, input.domain_volume) >= 0.5
                && list_p32(family, all_but_last, input.domain_volume) < 0.5
        );
    }

    #[test]
    fn test_p32_radii_list_starts_with_the_maximum_radius_for_large_fractures() {
        let (input, mut families) = parse_p32_family(
            r#"{ distribution = "truncated_power_law", alpha = 1.5, min = 1.0, max = 5.0 }"#,
            0.5,
        );

        families
            .generate_p32_radii(
                true,
                input.domain_volume,
                &input.region_volumes,
                &mut Mt64::new(1),
            )
            .unwrap();

        assert_eq!(families.families[0].radii_list[0], 5.0);
    }
}
