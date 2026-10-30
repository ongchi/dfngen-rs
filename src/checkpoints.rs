//! The checkpoint files of the binary in the output folder

use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::{Duration, Instant};

use dfngen::{Checkpoint, DfngenError, ProgressUpdate};

/// The checkpoint of the last write, for `--resume`
pub const CHECKPOINT_FILE: &str = "checkpoint.bin";
/// The accepted fractures at the last checkpoint, for a view of a stopped run
const PARTIAL_FRACTURES_FILE: &str = "fractures_partial.vtp";

/// Writes the checkpoints of a run into the output folder
pub struct CheckpointWriter {
    folder: PathBuf,
    interval: Duration,
    last_write: Instant,
}

impl CheckpointWriter {
    /// # Arguments
    ///
    /// * `folder` - The output folder
    /// * `interval` - Time between two checkpoints
    pub fn new(folder: &str, interval: Duration) -> Self {
        Self {
            folder: PathBuf::from(folder),
            interval,
            last_write: Instant::now(),
        }
    }

    /// True if the time from the last checkpoint is the interval or more
    pub fn is_due(&self) -> bool {
        self.last_write.elapsed() >= self.interval
    }

    /// Writes the checkpoint and the accepted fractures of the update.
    /// Each file is first written with another name, then renamed: an
    /// interruption does not leave a part of a file. Returns the path of the
    /// checkpoint.
    pub fn write(&mut self, update: &ProgressUpdate) -> Result<PathBuf, DfngenError> {
        fs::create_dir_all(&self.folder)?;

        let bytes = update.checkpoint().to_bytes()?;
        let checkpoint_path = self.folder.join(CHECKPOINT_FILE);
        write_then_rename(&checkpoint_path, |file| file.write_all(&bytes))?;
        write_then_rename(&self.folder.join(PARTIAL_FRACTURES_FILE), |file| {
            update.write_accepted_fractures(file)
        })?;

        self.last_write = Instant::now();
        Ok(checkpoint_path)
    }
}

/// Writes a file with the name plus `.part`, then renames it to the name
fn write_then_rename(
    path: &Path,
    write: impl FnOnce(&mut BufWriter<File>) -> std::io::Result<()>,
) -> Result<(), DfngenError> {
    let part_path = path.with_extension("part");
    let mut file = BufWriter::new(File::create(&part_path)?);
    write(&mut file)?;
    file.flush()?;
    drop(file);
    fs::rename(&part_path, path)?;
    Ok(())
}

/// Removes the checkpoint files of the output folder: the output of a
/// complete run does not need them
pub fn remove_checkpoint_files(folder: &str) -> Result<(), DfngenError> {
    for name in [CHECKPOINT_FILE, PARTIAL_FRACTURES_FILE] {
        let path = Path::new(folder).join(name);
        if path.exists() {
            fs::remove_file(path)?;
        }
    }
    Ok(())
}

/// Reads a checkpoint file
pub fn read_checkpoint(path: &str) -> Result<Checkpoint, DfngenError> {
    let bytes = fs::read(path).map_err(|source| DfngenError::ReadInput {
        file: path.to_string(),
        source,
    })?;
    Checkpoint::from_bytes(&bytes)
}
