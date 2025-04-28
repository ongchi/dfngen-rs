use std::{cell::RefCell, rc::Rc};

use parry3d_f64::bounding_volume::Aabb;
use parry3d_f64::na::Vector3;
use rand::Rng;
use rand_distr::Uniform;
use rand_mt::Mt64;
use tracing::error;

use crate::distribution::generating_points::random_position;
use crate::distribution::Fisher;
use crate::error::DfngenError;
use crate::structures::{RadiusDistribution, Shape};

use super::poly::Poly;

/// FractureFamily is used to hold varibales for all types of stochastic shapes. During getInput(),
/// all stochastic families for both recaangles and ellipses are parsed from the user input
/// and are placed in a Shape structure array.
pub struct FractureFamily {
    pub id: usize,

    pub shape: Shape,

    pub radius: RadiusDistribution,

    /// Current index to the radii list 'radiiList'.
    pub radii_idx: usize,

    /// Initial list of fracture/polygon radii, sorted largest to smallest.
    pub radii_list: Vec<f64>,

    /// Layer the family belongs to. 0 is entire domain, greater than 0 is a layer.
    /// e.g. 2 would be the second layer listed in the input file under "layers:".
    pub layer_id: usize,

    /// Region the family belongs to. 0 is entire domain, greater than 0 is a region.
    /// e.g. 2 would be the second region listed in the input file under "regions:".
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

    /// 'beta' is the rotation, or twist, around z normal before 3d rotation in radians
    /// or degrees depending on 'angleOption'.
    pub beta: f64,

    pub orientation: Fisher,

    pub boundary: Aabb,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum RadiusOption {
    FromCacheOrRng,
    FromRng,
    MaxRadius,
}

impl FractureFamily {
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
                                "distribution for {} family {} has been unable to generate a fracture with radius within set parameters after 1000 consecutive tries.",
                                &self.shape,
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
        let mut new_poly = match self.shape {
            Shape::Ellipse(n) => Poly::new_ell(n as usize, radius, self.aspect_ratio),
            Shape::Rectangle => Poly::new_rect(radius, self.aspect_ratio),
        };
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
        let position = random_position(self.boundary, generator.clone());
        new_poly.translate(position);

        new_poly
    }
}

#[derive(Default)]
pub struct FractureFamilyBuilder<'a> {
    id: Option<usize>,
    number_of_nodes: Option<u8>,
    radius: Option<RadiusDistribution>,
    orientation: Option<Fisher>,
    aspect_ratio: Option<f64>,
    beta: Option<f64>,
    p32_target: Option<f64>,
    layer_id: Option<usize>,
    region_id: Option<usize>,
    layers: Option<&'a [f64]>,
    regions: Option<&'a [f64]>,
    domain_size: Option<&'a Vector3<f64>>,
    domain_size_inc: Option<&'a Vector3<f64>>,
}

impl<'a> FractureFamilyBuilder<'a> {
    pub fn new() -> Self {
        Self {
            ..Default::default()
        }
    }

    pub fn id(&mut self, id: usize) -> &mut Self {
        self.id = Some(id);
        self
    }

    pub fn number_of_nodes(&mut self, number_of_nodes: u8) -> &mut Self {
        self.number_of_nodes = Some(number_of_nodes);
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

    pub fn layer_id(&mut self, layer: usize) -> &mut Self {
        self.layer_id = Some(layer);
        self
    }

    pub fn region_id(&mut self, region: usize) -> &mut Self {
        self.region_id = Some(region);
        self
    }

    pub fn layers(&mut self, layers: &'a [f64]) -> &mut Self {
        self.layers = Some(layers);
        self
    }

    pub fn regions(&mut self, regions: &'a [f64]) -> &mut Self {
        self.regions = Some(regions);
        self
    }

    pub fn domain_size(&mut self, domain_size: &'a Vector3<f64>) -> &mut Self {
        self.domain_size = Some(domain_size);
        self
    }

    pub fn domain_size_increase(&mut self, domain_size_inc: &'a Vector3<f64>) -> &mut Self {
        self.domain_size_inc = Some(domain_size_inc);
        self
    }

    pub fn build(&mut self) -> Result<FractureFamily, String> {
        let shape = self
            .number_of_nodes
            .map(Shape::Ellipse)
            .unwrap_or(Shape::Rectangle);

        let layers = self.layers.ok_or("layers is required".to_string())?;
        let regions = self.regions.ok_or("regions is required".to_string())?;
        let layer_id = self.layer_id.ok_or("layer_id is required".to_string())?;
        let region_id = self.region_id.ok_or("region_id is required".to_string())?;

        let domain_size = self
            .domain_size
            .ok_or("domain size is required".to_string())?;
        let domain_size_increase = self
            .domain_size_inc
            .ok_or("domain size increase is required".to_string())?;

        let (mins, maxs) = if layer_id == 0 && region_id == 0 {
            // The family layer is the whole domain
            let mins = (-domain_size - domain_size_increase) / 2.;
            let maxs = (domain_size + domain_size_increase) / 2.;
            (mins, maxs)
        } else if layer_id > 0 && region_id == 0 {
            // Family belongs to a certain layer, frac_fam.layer is > zero
            // Layers start at 1, but the array of layers start at 0, hence
            // the subtraction by 1
            // Layer 0 is reservered to be the entire domain
            let layer_idx = (layer_id - 1) * 2;
            let _mins = (-domain_size - domain_size_increase) / 2.;
            let _maxs = (domain_size + domain_size_increase) / 2.;
            // Layers only apply to z coordinates
            let mins = Vector3::new(_mins.x, _mins.y, layers[layer_idx]);
            let maxs = Vector3::new(_maxs.x, _maxs.y, layers[layer_idx + 1]);
            (mins, maxs)
        } else if layer_id == 0 && region_id > 0 {
            let region_idx = (region_id - 1) * 6;
            let mins = Vector3::new(
                regions[region_idx],
                regions[region_idx + 2],
                regions[region_idx + 4],
            );
            let maxs = Vector3::new(
                regions[region_idx + 1],
                regions[region_idx + 3],
                regions[region_idx + 5],
            );
            (mins, maxs)
        } else {
            return Err("Layer and Region both defined for this Family.".to_string());
        };

        let boundary = Aabb::new(mins.into(), maxs.into());

        Ok(FractureFamily {
            id: self.id.take().ok_or("id is required".to_string())?,
            shape,
            radius: self.radius.take().ok_or("radius is required".to_string())?,
            radii_idx: 0,
            radii_list: Vec::new(),
            layer_id,
            region_id,
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
            boundary,
        })
    }
}

pub struct FractureFamilyOption {
    pub families: Vec<FractureFamily>,

    // Each element is the probability of chosing a fracture from
    // the element's corresponding family to be inserted into the DFN.
    // The elements should add up to 1.0 (for %100).
    pub probabilities: Vec<f64>,

    // Copy of the original values of probabilities
    pub original_probabilities: Vec<f64>,
}

impl FractureFamilyOption {
    /// Generates a list of radii by distribution function of each family.
    pub fn generate_radii(
        &mut self,
        force_large_fractures: bool,
        n_poly: usize,
        generator: Rc<RefCell<Mt64>>,
    ) -> Result<(), DfngenError> {
        let mut n_poly = n_poly;
        if force_large_fractures {
            for shape in self.families.iter_mut() {
                let radius = shape.radius.max;
                shape.radii_list.push(radius);
            }
            n_poly -= self.families.len();
        }

        let mut ell_id = 0;
        let mut rect_id = 0;

        for i in 0..self.families.len() {
            let frac_fam = &mut self.families[i];

            match frac_fam.shape {
                Shape::Ellipse(_) => ell_id += 1,
                Shape::Rectangle => rect_id += 1,
            };

            for _ in 0..((self.probabilities[i] * n_poly as f64).ceil() as usize) {
                match frac_fam.radius.sample(generator.clone()) {
                    Ok(radius) => frac_fam.radii_list.push(radius),
                    Err(e) => {
                        error!(
                            "Attempt to populate fracture radii: {} Family {}",
                            frac_fam.shape,
                            match frac_fam.shape {
                                Shape::Ellipse(_) => ell_id,
                                Shape::Rectangle => rect_id,
                            }
                        );
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
