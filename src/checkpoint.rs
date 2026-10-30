//! Checkpoint of the insertion of the stochastic fractures: a restart from
//! it (`generation::resume`) gives the same network as a run with no stop.
//!
//! A rejected attempt does not change the network, thus the intersections,
//! the triple points and the cluster groups are a function of the accepted
//! fractures and their order. The checkpoint has the accepted stochastic
//! fractures, and the restart checks them again against the network. The
//! radii lists and the user-defined fractures come again from the input.

use parry3d_f64::na::{Point3, Vector3};
use serde::{Deserialize, Serialize};

use crate::error::DfngenError;
use crate::fracture::family_probability::cumulative_sum;
use crate::fracture::fracture_family::{FractureFamily, FractureFamilyCollection};
use crate::fracture::poly::{FractureOrigin, Poly};
use crate::fram::RejectReason;
use crate::random::GeneratorState;
use crate::stats::Stats;

/// Version of the format: a checkpoint of another version is an error
const FORMAT_VERSION: u32 = 1;

/// State of the insertion of the stochastic fractures, after an accepted fracture
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Checkpoint {
    /// `Input::text_hash` of the input of the run
    pub(crate) input_hash: u64,
    /// The seed of the run (an input with `seed = 0` gets a seed from the binary)
    pub(crate) seed: u64,
    pub(crate) generator: GeneratorState,
    /// The accepted stochastic fractures, in their order
    pub(crate) accepted_fractures: Vec<SavedFracture>,
    /// The rejected stochastic fractures (only with `fram.write_rejected_fractures`)
    pub(crate) rejected_fractures: Vec<(SavedFracture, RejectReason)>,
    /// One element for each family, in the order of the families
    pub(crate) families: Vec<FamilyState>,
    pub(crate) insertion: InsertionState,
    pub(crate) stats: Stats,
}

/// The values of the loop of the insertion
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub(crate) struct InsertionState {
    /// Probabilities of the families which have not met their P32 target
    /// (P32 option). The family probabilities of the input stay as they are.
    pub probabilities: Vec<f64>,
    /// Cumulative sum of `probabilities`
    pub cdf: Vec<f64>,
    /// True for each family which has its P32 target (P32 option)
    pub p32_status: Vec<bool>,
    /// Number of accepted stochastic fractures: `n_poly` does not include the
    /// user-defined fractures
    pub stochastic_accepted_count: usize,
}

/// The part of the state of a family which the insertion changes
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub(crate) struct FamilyState {
    pub radii_index: usize,
    pub current_p32: f64,
}

/// A stochastic fracture, before the check against the accepted fractures:
/// the fields which the check does not set, and which do not come from its
/// vertices (bounding box, area, boundary lines)
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub(crate) struct SavedFracture {
    family_index: usize,
    x_radius: f64,
    y_radius: f64,
    aspect_ratio: f64,
    translation: [f64; 3],
    normal: [f64; 3],
    vertices: Vec<[f64; 3]>,
    region_ids: Vec<usize>,
}

impl SavedFracture {
    /// The saved fields of a stochastic fracture; `None` for a user-defined fracture
    pub fn of(poly: &Poly) -> Option<Self> {
        Some(Self {
            family_index: poly.origin.family_index()?,
            x_radius: poly.x_radius,
            y_radius: poly.y_radius,
            aspect_ratio: poly.aspect_ratio,
            translation: poly.translation.into(),
            normal: poly.normal.into(),
            vertices: poly
                .vertices
                .iter()
                .map(|vertex| vertex.coords.into())
                .collect(),
            region_ids: poly.region_ids.clone(),
        })
    }

    /// The fracture as it was before the check against the accepted
    /// fractures, with no bounding box, area or boundary lines
    pub fn to_poly(&self) -> Poly {
        Poly {
            origin: FractureOrigin::Family(self.family_index),
            x_radius: self.x_radius,
            y_radius: self.y_radius,
            aspect_ratio: self.aspect_ratio,
            translation: Vector3::from(self.translation),
            normal: Vector3::from(self.normal),
            vertices: self
                .vertices
                .iter()
                .map(|&vertex| Point3::from(vertex))
                .collect(),
            region_ids: self.region_ids.clone(),
            ..Poly::default()
        }
    }
}

impl InsertionState {
    /// The state at the start of the insertion
    pub fn new(family_collection: &FractureFamilyCollection) -> Self {
        // A copy: the family probabilities must stay in their original state
        let probabilities = family_collection.probabilities.clone();
        Self {
            cdf: cumulative_sum(&probabilities),
            probabilities,
            p32_status: vec![false; family_collection.families.len()],
            stochastic_accepted_count: 0,
        }
    }
}

impl FamilyState {
    pub fn of(family: &FractureFamily) -> Self {
        Self {
            radii_index: family.radii_index,
            current_p32: family.current_p32,
        }
    }

    pub fn apply_to(&self, family: &mut FractureFamily) {
        family.radii_index = self.radii_index;
        family.current_p32 = self.current_p32;
    }
}

impl Checkpoint {
    /// The checkpoint in its binary format
    pub fn to_bytes(&self) -> Result<Vec<u8>, DfngenError> {
        let mut bytes =
            postcard::to_allocvec(&FORMAT_VERSION).map_err(|error| format_error(&error))?;
        bytes.extend(postcard::to_allocvec(self).map_err(|error| format_error(&error))?);
        Ok(bytes)
    }

    /// Reads a checkpoint of `to_bytes`
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, DfngenError> {
        let (version, rest) =
            postcard::take_from_bytes::<u32>(bytes).map_err(|error| format_error(&error))?;
        if version != FORMAT_VERSION {
            return Err(DfngenError::CheckpointFormat(format!(
                "format version {version}; this version of dfngen reads version {FORMAT_VERSION}"
            )));
        }
        let (checkpoint, rest) =
            postcard::take_from_bytes(rest).map_err(|error| format_error(&error))?;
        if !rest.is_empty() {
            return Err(DfngenError::CheckpointFormat(format!(
                "{} bytes after the end of the checkpoint",
                rest.len()
            )));
        }
        Ok(checkpoint)
    }

    /// Number of accepted stochastic fractures
    pub fn accepted_count(&self) -> usize {
        self.insertion.stochastic_accepted_count
    }
}

fn format_error(error: &postcard::Error) -> DfngenError {
    DfngenError::CheckpointFormat(error.to_string())
}

/// FNV-1a hash of the texts of the input files. It is the same in each
/// version of Rust (the hash of `std` is not).
pub fn text_hash(texts: &[&str]) -> u64 {
    const OFFSET_BASIS: u64 = 0xcbf2_9ce4_8422_2325;
    const PRIME: u64 = 0x0000_0100_0000_01b3;

    // The length before each text: ("ab", "c") and ("a", "bc") differ
    texts
        .iter()
        .flat_map(|text| {
            (text.len() as u64)
                .to_le_bytes()
                .into_iter()
                .chain(text.bytes())
        })
        .fold(OFFSET_BASIS, |hash, byte| {
            (hash ^ u64::from(byte)).wrapping_mul(PRIME)
        })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fracture::poly::FractureOrigin;

    fn checkpoint() -> Checkpoint {
        let poly = Poly {
            origin: FractureOrigin::Family(1),
            x_radius: 1.5,
            y_radius: 0.75,
            aspect_ratio: 0.5,
            translation: Vector3::new(0.1, -0.2, 1. / 3.),
            normal: Vector3::new(0., 0.6, 0.8),
            vertices: vec![Point3::new(1., 2., 3.), Point3::new(-1., 0.1, 1e-300)],
            region_ids: vec![2],
            ..Poly::default()
        };
        let fracture = SavedFracture::of(&poly).unwrap();
        Checkpoint {
            input_hash: 42,
            seed: 12345,
            generator: GeneratorState::RecentOutputs((0..312).collect()),
            accepted_fractures: vec![fracture.clone()],
            rejected_fractures: vec![(fracture, RejectReason::CloseToEdge)],
            families: vec![FamilyState {
                radii_index: 3,
                current_p32: 0.1 + 0.2,
            }],
            insertion: InsertionState {
                probabilities: vec![0.4, 0.6],
                cdf: vec![0.4, 1.],
                p32_status: vec![false, true],
                stochastic_accepted_count: 1,
            },
            stats: Stats::default(),
        }
    }

    #[test]
    fn test_bytes_give_the_same_checkpoint() {
        let checkpoint = checkpoint();
        let bytes = checkpoint.to_bytes().unwrap();
        assert_eq!(Checkpoint::from_bytes(&bytes).unwrap(), checkpoint);
    }

    #[test]
    fn test_other_format_version_is_an_error() {
        let mut bytes = checkpoint().to_bytes().unwrap();
        bytes[0] = 2;
        assert!(matches!(
            Checkpoint::from_bytes(&bytes),
            Err(DfngenError::CheckpointFormat(_))
        ));
    }

    #[test]
    fn test_cut_bytes_are_an_error() {
        let bytes = checkpoint().to_bytes().unwrap();
        assert!(matches!(
            Checkpoint::from_bytes(&bytes[..bytes.len() / 2]),
            Err(DfngenError::CheckpointFormat(_))
        ));
    }

    #[test]
    fn test_saved_fracture_gives_the_fields_of_the_fracture() {
        let poly = SavedFracture::of(&Poly {
            origin: FractureOrigin::Family(0),
            x_radius: 2.,
            vertices: vec![Point3::new(1., 2., 3.)],
            ..Poly::default()
        })
        .unwrap()
        .to_poly();
        assert_eq!(
            (poly.origin, poly.x_radius, poly.vertices),
            (FractureOrigin::Family(0), 2., vec![Point3::new(1., 2., 3.)])
        );
    }

    #[test]
    fn test_user_defined_fracture_is_not_saved() {
        let poly = Poly {
            origin: FractureOrigin::UserDefined,
            ..Poly::default()
        };
        assert!(SavedFracture::of(&poly).is_none());
    }

    #[test]
    fn test_text_hash_depends_on_the_text_boundaries() {
        assert_ne!(text_hash(&["ab", "c"]), text_hash(&["a", "bc"]));
    }
}
