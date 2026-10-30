//! The random generator of the generation, with a state which a checkpoint
//! can save

use rand::RngCore;
use rand_mt::Mt64;
use serde::{Deserialize, Serialize};

/// Number of outputs of `Mt64` which give its state (`Mt64::recover`)
const STATE_OUTPUT_COUNT: usize = 312;

/// State of a `Generator`, from which `Generator::from_state` makes it again
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum GeneratorState {
    /// Fewer than 312 outputs: the seed and the number of outputs
    Start { seed: u64, output_count: usize },
    /// The last 312 outputs, oldest first
    RecentOutputs(Vec<u64>),
}

/// The outputs which a `Generator` keeps
#[derive(Clone)]
enum OutputHistory {
    /// Fewer than 312 outputs: the seed and all outputs
    Start { seed: u64, outputs: Vec<u64> },
    /// The last 312 outputs; `oldest_position` is the oldest one
    Full {
        outputs: Vec<u64>,
        oldest_position: usize,
    },
}

/// `Mt64` which keeps its last 312 outputs: they give its state, which
/// `Mt64` does not show. Each call uses whole 64-bit outputs, as `Mt64`
/// does, thus the random numbers are those of `Mt64`.
#[derive(Clone)]
pub struct Generator {
    mersenne_twister: Mt64,
    history: OutputHistory,
}

impl Generator {
    pub fn new(seed: u64) -> Self {
        Self {
            mersenne_twister: Mt64::new(seed),
            history: OutputHistory::Start {
                seed,
                outputs: Vec::with_capacity(STATE_OUTPUT_COUNT),
            },
        }
    }

    /// The state, from which `from_state` makes the same generator
    pub fn state(&self) -> GeneratorState {
        match &self.history {
            OutputHistory::Start { seed, outputs } => GeneratorState::Start {
                seed: *seed,
                output_count: outputs.len(),
            },
            OutputHistory::Full {
                outputs,
                oldest_position,
            } => {
                let (newest, oldest) = outputs.split_at(*oldest_position);
                GeneratorState::RecentOutputs([oldest, newest].concat())
            }
        }
    }

    /// The generator of a state from `state`. It gives the same next outputs
    /// as the generator which gave the state. `None` if the state does not
    /// have 312 outputs.
    pub fn from_state(state: &GeneratorState) -> Option<Self> {
        match state {
            GeneratorState::Start { seed, output_count } => {
                let mut generator = Self::new(*seed);
                for _ in 0..*output_count {
                    generator.next_u64();
                }
                Some(generator)
            }
            GeneratorState::RecentOutputs(outputs) => Some(Self {
                mersenne_twister: Mt64::recover(outputs.iter().copied()).ok()?,
                history: OutputHistory::Full {
                    outputs: outputs.clone(),
                    oldest_position: 0,
                },
            }),
        }
    }
}

impl RngCore for Generator {
    fn next_u64(&mut self) -> u64 {
        let output = self.mersenne_twister.next_u64();
        match &mut self.history {
            OutputHistory::Start { outputs, .. } if outputs.len() + 1 < STATE_OUTPUT_COUNT => {
                outputs.push(output);
            }
            OutputHistory::Start { outputs, .. } => {
                let mut outputs = std::mem::take(outputs);
                outputs.push(output);
                self.history = OutputHistory::Full {
                    outputs,
                    oldest_position: 0,
                };
            }
            OutputHistory::Full {
                outputs,
                oldest_position,
            } => {
                outputs[*oldest_position] = output;
                *oldest_position = (*oldest_position + 1) % STATE_OUTPUT_COUNT;
            }
        }
        output
    }

    /// The low 32 bits of a 64-bit output, as `Mt64::next_u32`
    fn next_u32(&mut self) -> u32 {
        self.next_u64() as u32
    }

    /// One 64-bit output for each 8 bytes (little endian), and one for the
    /// rest, as `Mt64::fill_bytes`
    fn fill_bytes(&mut self, destination: &mut [u8]) {
        for chunk in destination.chunks_mut(8) {
            let bytes = self.next_u64().to_le_bytes();
            chunk.copy_from_slice(&bytes[..chunk.len()]);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The next 1000 outputs of a generator after `output_count` outputs, and
    /// those of the generator of its state
    fn next_outputs_after_restore(output_count: usize) -> (Vec<u64>, Vec<u64>) {
        let mut generator = Generator::new(12345);
        for _ in 0..output_count {
            generator.next_u64();
        }
        let mut restored = Generator::from_state(&generator.state()).unwrap();
        let expected = (0..1000).map(|_| generator.next_u64()).collect();
        let actual = (0..1000).map(|_| restored.next_u64()).collect();
        (expected, actual)
    }

    #[test]
    fn test_outputs_are_those_of_mt64() {
        let mut mersenne_twister = Mt64::new(7);
        let mut generator = Generator::new(7);
        let mut bytes = [0; 13];
        generator.fill_bytes(&mut bytes);
        let first = mersenne_twister.next_u64().to_le_bytes();
        let second = mersenne_twister.next_u64().to_le_bytes();
        assert_eq!(bytes, [&first[..], &second[..5]].concat().as_slice());
        assert_eq!(generator.next_u32(), mersenne_twister.next_u32());
    }

    #[test]
    fn test_restore_before_312_outputs_gives_the_same_outputs() {
        let (expected, actual) = next_outputs_after_restore(100);
        assert_eq!(actual, expected);
    }

    #[test]
    fn test_restore_at_312_outputs_gives_the_same_outputs() {
        let (expected, actual) = next_outputs_after_restore(STATE_OUTPUT_COUNT);
        assert_eq!(actual, expected);
    }

    #[test]
    fn test_restore_after_many_outputs_gives_the_same_outputs() {
        // Not a multiple of 312: the oldest output is not at position 0
        let (expected, actual) = next_outputs_after_restore(5000);
        assert_eq!(actual, expected);
    }

    #[test]
    fn test_restore_of_a_restored_generator_gives_the_same_outputs() {
        let mut generator = Generator::new(3);
        for _ in 0..400 {
            generator.next_u64();
        }
        let mut restored = Generator::from_state(&generator.state()).unwrap();
        for _ in 0..50 {
            restored.next_u64();
            generator.next_u64();
        }
        let mut restored_again = Generator::from_state(&restored.state()).unwrap();
        assert_eq!(restored_again.next_u64(), generator.next_u64());
    }

    #[test]
    fn test_state_with_too_few_outputs_gives_none() {
        assert!(Generator::from_state(&GeneratorState::RecentOutputs(vec![1; 10])).is_none());
    }
}
