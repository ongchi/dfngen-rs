use std::rc::Rc;

use parry3d_f64::na::Vector3;
use rand::Rng;
use rand_distr::Uniform;
use rand_mt::Mt64;

use crate::distribution::Fisher;
use crate::domain_mesh::DomainMesh;
use crate::error::DfngenError;
use crate::structures::RadiusDistribution;

use super::poly::{FractureOrigin, Poly};

/// FractureFamily holds the variables of one stochastic family of the input file.
pub struct FractureFamily {
    pub id: usize,

    /// Number of vertices of each fracture.
    pub number_of_points: usize,

    pub radius: RadiusDistribution,

    /// Current index to the radii list 'radiiList'.
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

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum RadiusOption {
    FromCacheOrRng,
    FromRng,
    MaxRadius,
}

impl FractureFamily {
    /// Adds the area of an accepted polygon to the family's current P32.
    /// The volume is that of the family's region, or of the whole domain.
    pub fn add_area_to_p32(&mut self, area: f64, domain_volume: f64, region_volumes: &[f64]) {
        let volume = if self.region_id > 0 {
            region_volumes[self.region_id - 1]
        } else {
            domain_volume
        };

        self.current_p32 += area * 2. / volume;
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

    pub fn normal_vector(&self, rng: &mut Mt64) -> Vector3<f64> {
        rng.sample(&self.orientation)
    }

    pub fn create_poly(
        &mut self,
        eps: f64,
        family_index: usize,
        radius_opt: RadiusOption,
        generator: &mut Mt64,
    ) -> Result<Poly, DfngenError> {
        let radius = match radius_opt {
            RadiusOption::FromCacheOrRng | RadiusOption::FromRng => {
                if self.radii_index >= self.radii_list.len() || radius_opt == RadiusOption::FromRng
                {
                    // If out of radii from list, insert random radius
                    self.sample_radius(generator)?
                } else {
                    // Insert radius from list
                    let radius = self.radii_list[self.radii_index];
                    self.radii_index += 1;
                    radius
                }
            }
            RadiusOption::MaxRadius => self.radius.max,
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

        // Save newPoly's new normal vector
        new_poly.normal = norm;

        // Translate - will also set translation vector in poly structure
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

    pub fn build(&mut self) -> Result<FractureFamily, String> {
        Ok(FractureFamily {
            id: self.id.take().ok_or("id is required".to_string())?,
            number_of_points: self
                .number_of_points
                .ok_or("number of points is required".to_string())?,
            radius: self.radius.take().ok_or("radius is required".to_string())?,
            radii_index: 0,
            radii_list: Vec::new(),
            region_id: self.region_id.ok_or("region_id is required".to_string())?,
            aspect_ratio: self
                .aspect_ratio
                .ok_or("aspect ratio is required".to_string())?,
            p32_target: self.p32_target.unwrap_or(0.),
            current_p32: 0.,
            beta: self.beta,
            orientation: self
                .orientation
                .take()
                .ok_or("orientation is required".to_string())?,
            domain_mesh: self
                .domain_mesh
                .take()
                .ok_or("domain mesh is required".to_string())?,
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

    /// Sort each family's radii list from largest to smallest.
    /// This will allow the DFN gereration to start from largest to smallest
    /// fractures.
    pub fn sort_radii(&mut self) {
        for family in &mut self.families {
            family.radii_list.sort_by(|a, b| b.total_cmp(a));
        }
    }
}
