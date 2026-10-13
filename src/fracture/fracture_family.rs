use std::{cell::RefCell, rc::Rc};

use parry3d_f64::na::Vector3;
use rand::Rng;
use rand_distr::Uniform;
use rand_mt::Mt64;
use tracing::error;

use crate::distribution::Fisher;
use crate::domain_mesh::DomainMesh;
use crate::error::DfngenError;
use crate::structures::RadiusDistribution;

use super::poly::Poly;

/// FractureFamily holds the variables of one stochastic family of the input file.
pub struct FractureFamily {
    pub id: usize,

    /// Number of vertices of each fracture.
    pub number_of_points: usize,

    pub radius: RadiusDistribution,

    /// Current index to the radii list 'radiiList'.
    pub radii_idx: usize,

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

    /// 'betaOption' is the rotation about the polygon's  normal vector
    /// True - User Specified Rotation
    /// False - Uniform Distribution
    pub beta_distribution: bool,

    /// 'beta' is the rotation, or twist, around z normal before 3d rotation in radians.
    pub beta: f64,

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
    pub fn sample_position(&self, generator: Rc<RefCell<Mt64>>) -> Vector3<f64> {
        self.domain_mesh
            .sample_point(self.region_id, &mut *generator.borrow_mut())
            .coords
    }

    pub fn normal_vector(&self, rng: Rc<RefCell<Mt64>>) -> Vector3<f64> {
        rng.borrow_mut().sample(&self.orientation)
    }

    pub fn create_poly(
        &mut self,
        eps: f64,
        fam_idx: usize,
        radius_opt: RadiusOption,
        generator: Rc<RefCell<Mt64>>,
    ) -> Poly {
        let radius = match radius_opt {
            RadiusOption::FromCacheOrRng | RadiusOption::FromRng => {
                if self.radii_idx >= self.radii_list.len() || radius_opt == RadiusOption::FromRng {
                    // If out of radii from list, insert random radius
                    match self.radius.sample(generator.clone()) {
                        Ok(radius) => radius,
                        Err(e) => {
                            error!("");
                            error!("{}", e);
                            panic!(
                                "distribution for family {} has been unable to generate a fracture with radius within set parameters after 1000 consecutive tries.",
                                self.id
                            );
                        }
                    }
                } else {
                    // Insert radius from list
                    let radius = self.radii_list[self.radii_idx];
                    self.radii_idx += 1;
                    radius
                }
            }
            RadiusOption::MaxRadius => self.radius.max,
        };

        // New polygon to build
        let mut new_poly = Poly::new(self.number_of_points, radius, self.aspect_ratio);
        new_poly.family_num = fam_idx as isize;

        // Initialize beta based on distrubution type: 0 = unifrom on [0,2PI], 1 = constant
        let beta = if !self.beta_distribution {
            let uniform_dist = Uniform::new(0., 2. * std::f64::consts::PI).unwrap();
            generator.borrow_mut().sample(uniform_dist)
        } else {
            self.beta
        };

        // Apply 2d rotation matrix, twist around origin
        // assumes polygon on x-y plane
        // Angle must be in rad
        new_poly.rotation_2d(beta);

        // Rotate vertices to norm (new normal)
        let norm = self.normal_vector(generator.clone()).normalize();
        new_poly.rotation_3d(&norm, eps);

        // Save newPoly's new normal vector
        new_poly.normal = norm;

        // Translate - will also set translation vector in poly structure
        let position = self.sample_position(generator.clone());
        new_poly.translate(position);

        new_poly
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
            radii_idx: 0,
            radii_list: Vec::new(),
            region_id: self.region_id.ok_or("region_id is required".to_string())?,
            aspect_ratio: self
                .aspect_ratio
                .ok_or("aspect ratio is required".to_string())?,
            p32_target: self.p32_target.unwrap_or(0.),
            current_p32: 0.,
            beta_distribution: self.beta.is_some(),
            beta: self.beta.unwrap_or(0.),
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
        generator: Rc<RefCell<Mt64>>,
    ) -> Result<(), DfngenError> {
        if force_large_fractures {
            for shape in self.families.iter_mut() {
                let radius = shape.radius.max;
                shape.radii_list.push(radius);
            }
            n_poly -= self.families.len();
        }

        for (frac_fam, probability) in self.families.iter_mut().zip(&self.probabilities) {
            for _ in 0..((probability * n_poly as f64).ceil() as usize) {
                match frac_fam.radius.sample(generator.clone()) {
                    Ok(radius) => frac_fam.radii_list.push(radius),
                    Err(e) => {
                        error!("Attempt to populate fracture radii: Family {}", frac_fam.id);
                        error!("{}", e);
                        return Err(e);
                    }
                }
            }
        }

        Ok(())
    }

    /// Sort each family's radii list from largest to smallest.
    /// This will allow the DFN gereration to start from largest to smallest
    /// fractures.
    pub fn sort_radii(&mut self) {
        for ff in self.families.iter_mut() {
            ff.radii_list.sort_by(|a, b| b.partial_cmp(a).unwrap())
        }
    }
}
