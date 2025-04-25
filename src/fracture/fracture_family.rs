use std::{cell::RefCell, rc::Rc};

use parry3d_f64::bounding_volume::Aabb;
use parry3d_f64::na::Vector3;
use rand::Rng;
use rand_distr::Uniform;
use rand_mt::Mt64;

use crate::computational_geometry::{apply_rotation2_d, apply_rotation3_d, translate};
use crate::distribution::generating_points::random_translation;
use crate::distribution::{generating_points::generate_theta, Fisher};
use crate::distribution::{TruncExp, TruncLogNormal, TruncPowerLaw};
use crate::error::DfngenError;
use crate::structures::{Poly, RadiusDistribution, RadiusFunction, Shape};

use super::insert_shape::{
    get_family_number, initialize_ell_vertices, initialize_rect_vertices, poly_boundary,
};

/// FractureFamily is used to hold varibales for all types of stochastic shapes. During getInput(),
/// all stochastic families for both recaangles and ellipses are parsed from the user input
/// and are placed in a Shape structure array.
pub struct FractureFamily {
    pub shape: Shape,

    pub radius: RadiusDistribution,

    /// Array of thetas to build poly from, initialized while reading input and building shape structures
    pub theta_list: Vec<f64>,

    /// Current index to the radii list 'radiiList'.
    pub radii_idx: usize,

    /// Initial list of fracture/polygon radii, sorted largest to smallest.
    pub radii_list: Vec<f64>,

    /// Layer the family belongs to. 0 is entire domain, greater than 0 is a layer.
    /// e.g. 2 would be the second layer listed in the input file under "layers:".
    pub layer: usize,

    /// Region the family belongs to. 0 is entire domain, greater than 0 is a region.
    /// e.g. 2 would be the second region listed in the input file under "regions:".
    pub region: usize,

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
}

#[derive(Debug, Clone, Copy)]
pub enum RadiusOption {
    FromCacheOrRng,
    FromRng,
    MaxRadius,
}

impl FractureFamily {
    pub fn normal_vector(&self, rng: Rc<RefCell<Mt64>>) -> Vector3<f64> {
        rng.borrow_mut().sample(&self.orientation)
    }

    #[allow(clippy::too_many_arguments)]
    pub fn create_poly(
        &mut self,
        h: f64,
        eps: f64,
        n_fam_ell: usize,
        fam_idx: usize,
        radius_opt: RadiusOption,
        domain_size: &Vector3<f64>,
        domain_size_inc: &Vector3<f64>,
        layers: &[f64],
        regions: &[f64],
        generator: Rc<RefCell<Mt64>>,
    ) -> Poly {
        let radius = match radius_opt {
            RadiusOption::FromCacheOrRng => generate_radius(
                h,
                n_fam_ell,
                self,
                generator.clone(),
                fam_idx as isize,
                true,
            ),
            RadiusOption::FromRng => generate_radius(
                h,
                n_fam_ell,
                self,
                generator.clone(),
                fam_idx as isize,
                false,
            ),
            RadiusOption::MaxRadius => self.radius.max,
        };

        let boundary = poly_boundary(domain_size, domain_size_inc, layers, regions, self);

        generate_poly_with_radius(
            eps,
            radius,
            self,
            boundary,
            generator.clone(),
            fam_idx as isize,
        )
    }
}

#[derive(Default)]
pub struct FractureFamilyBuilder {
    number_of_nodes: Option<u8>,
    radius: Option<RadiusDistribution>,
    orientation: Option<Fisher>,
    aspect_ratio: Option<f64>,
    beta: Option<f64>,
    p32_target: Option<f64>,
    layer: Option<usize>,
    region: Option<usize>,
}

impl FractureFamilyBuilder {
    pub fn new() -> Self {
        Self {
            ..Default::default()
        }
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

    pub fn layer(&mut self, layer: usize) -> &mut Self {
        self.layer = Some(layer);
        self
    }

    pub fn region(&mut self, region: usize) -> &mut Self {
        self.region = Some(region);
        self
    }

    pub fn build(&mut self) -> Result<FractureFamily, String> {
        let frac_family = self
            .number_of_nodes
            .map(Shape::Ellipse)
            .unwrap_or(Shape::Rectangle);
        let aspect_ratio = self
            .aspect_ratio
            .ok_or("aspect ration is required".to_string())?;
        let theta_list = generate_theta(
            aspect_ratio,
            match frac_family {
                Shape::Ellipse(n) => n as usize,
                Shape::Rectangle => 4,
            },
        );

        Ok(FractureFamily {
            shape: frac_family,
            radius: self.radius.take().ok_or("radius is required".to_string())?,
            theta_list,
            radii_idx: 0,
            radii_list: Vec::new(),
            layer: self.layer.ok_or("layer is required".to_string())?,
            region: self.region.ok_or("region is required".to_string())?,
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

            match frac_fam.radius.sample(
                (self.probabilities[i] * n_poly as f64).ceil() as usize,
                generator.clone(),
            ) {
                Ok(radii) => frac_fam.radii_list.extend(radii),
                Err(_) => {
                    return Err(DfngenError::TooManySmallFractures {
                        shape: frac_fam.shape,
                        id: match frac_fam.shape {
                            Shape::Ellipse(_) => ell_id,
                            Shape::Rectangle => rect_id,
                        },
                    });
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

pub fn generate_radius(
    h: f64,
    n_fam_ell: usize,
    frac_fam: &mut FractureFamily,
    generator: Rc<RefCell<Mt64>>,
    family_index: isize,
    use_list: bool,
) -> f64 {
    let radius_distr = &frac_fam.radius;
    match radius_distr.function {
        RadiusFunction::Constant(c) => c,
        ref others => {
            if frac_fam.radii_idx >= frac_fam.radii_list.len() || !use_list {
                // If out of radii from list, insert random radius
                let distr = match others {
                    RadiusFunction::LogNormal { mu, sigma } => {
                        let min_val = f64::max(h, radius_distr.min);
                        let distr =
                            TruncLogNormal::new(min_val, radius_distr.max, *mu, *sigma).unwrap();
                        generator.clone().borrow_mut().sample(distr)
                    }
                    RadiusFunction::TruncatedPowerLaw { alpha } => {
                        let distr = TruncPowerLaw::new(radius_distr.min, radius_distr.max, *alpha);
                        Ok(generator.clone().borrow_mut().sample(distr))
                    }
                    RadiusFunction::Exponential { lambda } => {
                        let min_val = f64::max(h, radius_distr.min);
                        let distr = TruncExp::new(min_val, radius_distr.max, *lambda).unwrap();

                        generator.clone().borrow_mut().sample(distr)
                    }
                    RadiusFunction::Constant(_) => unreachable!(),
                };

                match distr {
                            Ok(radius) => radius,
                            Err(_) => panic!(
                                    "distribution for {} family {} has been unable to generate a fracture with radius within set parameters after 1000 consecutive tries.",
                                    &frac_fam.shape,
                                    get_family_number(n_fam_ell, family_index, frac_fam.shape),
                                )
                        }
            } else {
                // Insert radius from list
                let radius = frac_fam.radii_list[frac_fam.radii_idx];
                frac_fam.radii_idx += 1;
                radius
            }
        }
    }
}

/// Generate Polygon/Fracture With Given Radius
///
/// Similar to generatePoly() except the radius is passed to the function.
/// Generates a polygon
/// Shape (ell or rect) still comes from the fractures' familiy
///
/// NOTE: Function does not create bouding box. The bouding box has to be
///       created after fracture truncation
///
/// # Arguments
///
/// * `orientation_option` - Orientation option
/// * `eps` - Epsilon value for floating point comparisons
/// * `radius` - Radius for polygon
/// * `fracture_fam` - fracture family to generate fracture from
/// * `bbox` - Bounding box
/// * `generator` - Random generator
/// * `family_index` - Index of 'FractureFamily' in the frac_fam array in main()
///
/// # Returns
///
/// Polygon with radius passed in arg 1 and shape based on `frac_fam`
pub fn generate_poly_with_radius(
    eps: f64,
    radius: f64,
    frac_fam: &FractureFamily,
    bbox: Aabb,
    generator: Rc<RefCell<Mt64>>,
    family_index: isize,
) -> Poly {
    // New polygon to build
    let mut new_poly = Poly::default();
    // Initialize normal to {0,0,1}. ( All polys start on x-y plane )
    new_poly.normal = Vector3::new(0., 0., 1.);
    new_poly.number_of_nodes = frac_fam.shape.number_of_nodes() as isize;
    new_poly.vertices = Vec::with_capacity((3 * new_poly.number_of_nodes) as usize); //numPoints*{x,y,z}
    new_poly.family_num = family_index;

    match frac_fam.shape {
        Shape::Ellipse(n) => {
            initialize_ell_vertices(
                &mut new_poly,
                radius,
                frac_fam.aspect_ratio,
                &frac_fam.theta_list,
                n as usize,
            );
        }
        Shape::Rectangle => {
            initialize_rect_vertices(&mut new_poly, radius, frac_fam.aspect_ratio);
        }
    }

    // Initialize beta based on distrubution type: 0 = unifrom on [0,2PI], 1 = constant
    let beta = if !frac_fam.beta_distribution {
        let uniform_dist = Uniform::new(0., 2. * std::f64::consts::PI).unwrap();
        generator.borrow_mut().sample(uniform_dist)
    } else {
        frac_fam.beta
    };

    // Apply 2d rotation matrix, twist around origin
    // assumes polygon on x-y plane
    // Angle must be in rad
    apply_rotation2_d(&mut new_poly, beta);

    // Rotate vertices to norm (new normal)
    let norm = frac_fam.normal_vector(generator.clone()).normalize();
    apply_rotation3_d(&mut new_poly, &norm, eps);

    // Save newPoly's new normal vector
    new_poly.normal = norm;

    let translation = random_translation(
        generator.clone(),
        bbox.mins.x,
        bbox.maxs.x,
        bbox.mins.y,
        bbox.maxs.y,
        bbox.mins.z,
        bbox.maxs.z,
    );

    // Translate - will also set translation vector in poly structure
    translate(&mut new_poly, translation);

    new_poly
}
