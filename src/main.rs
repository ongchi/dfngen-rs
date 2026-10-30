use std::process::ExitCode;
use std::sync::atomic::{AtomicBool, Ordering};
use std::time::{Duration, SystemTime, UNIX_EPOCH};

use clap::Parser;
use dfngen::{
    generate, read_input, resume, write_output_to_folder, DfngenError, ProgressUpdate, Resume,
    ResumeMode,
};
use tracing::{error, info, warn};
use tracing_subscriber::{layer::SubscriberExt, util::SubscriberInitExt};

use crate::checkpoints::{read_checkpoint, remove_checkpoint_files, CheckpointWriter};
use crate::progress::{GenerationProgress, LogWriter};

mod checkpoints;
mod progress;

/// Exit code of a run which Ctrl-C stopped (128 + SIGINT, as a shell gives)
const STOPPED_EXIT_CODE: i32 = 130;

/// Set by the first Ctrl-C: the run writes a checkpoint at the next accepted
/// fracture, then stops
static IS_STOP_REQUESTED: AtomicBool = AtomicBool::new(false);
/// Set after the last accepted fracture: the last checkpoint has all of the
/// insertion, thus a Ctrl-C stops at once
static IS_INSERTION_COMPLETE: AtomicBool = AtomicBool::new(false);

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    /// Path to input file
    input_file: String,

    /// Path to output folder
    output_folder: String,

    /// Seconds between two checkpoints of the insertion of the stochastic
    /// fractures (checkpoint.bin and fractures_partial.vtp in the output
    /// folder; a complete run removes them)
    #[arg(long, value_name = "SECONDS", default_value_t = 60, value_parser = clap::value_parser!(u64).range(1..))]
    checkpoint_every: u64,

    /// Continue the stopped run of this checkpoint. The input files must be
    /// those of the stopped run.
    #[arg(long, value_name = "CHECKPOINT")]
    resume: Option<String>,

    /// With --resume: make the network from the fractures of the checkpoint,
    /// with no more insertion
    #[arg(long, requires = "resume")]
    finish: bool,
}

fn main() -> ExitCode {
    let cli = Cli::parse();

    init_tracing();

    if let Err(e) = ctrlc::set_handler(stop_on_ctrl_c) {
        warn!("Ctrl-C stops the run with no checkpoint: {}", e);
    }

    info!("Starting dfngen-rs");

    match run(&cli) {
        Ok(()) => {
            info!("dfngen-rs - Complete");
            ExitCode::SUCCESS
        }
        Err(e) => {
            error!("{}", e);
            ExitCode::FAILURE
        }
    }
}

/// The first Ctrl-C asks for a stop at the next accepted fracture; a second
/// one, or a Ctrl-C after the insertion, stops at once
fn stop_on_ctrl_c() {
    if IS_INSERTION_COMPLETE.load(Ordering::SeqCst)
        || IS_STOP_REQUESTED.swap(true, Ordering::SeqCst)
    {
        std::process::exit(STOPPED_EXIT_CODE);
    }
    warn!("Stop: a checkpoint is written at the next accepted fracture, then the run stops. Press Ctrl-C again to stop at once.");
}

/// Reads the input, generates the network and writes the output files.
fn run(cli: &Cli) -> Result<(), DfngenError> {
    let (mut input, family_collection) = read_input(&cli.input_file)?;

    // A restart uses the seed of its checkpoint
    if input.seed == 0 && cli.resume.is_none() {
        input.seed = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .map(|t| t.as_secs())?;
    }

    // The bar starts with the first accepted stochastic fracture
    let stop_condition = input.stop_condition;
    let mut progress = None;
    let mut checkpoints = CheckpointWriter::new(
        &cli.output_folder,
        Duration::from_secs(cli.checkpoint_every),
    );
    let mut on_progress = |update: &ProgressUpdate| {
        let bar = progress.get_or_insert_with(|| GenerationProgress::new(stop_condition));
        bar.update(update);

        let is_stop = IS_STOP_REQUESTED.load(Ordering::SeqCst);
        let is_last = update.fraction >= 1.;
        if is_stop || is_last || checkpoints.is_due() {
            match checkpoints.write(update) {
                Ok(path) if is_stop => {
                    bar.clear();
                    info!(
                        "Stopped after {} accepted stochastic fractures. Continue with: --resume {}",
                        update.accepted_count,
                        path.display()
                    );
                    std::process::exit(STOPPED_EXIT_CODE);
                }
                Ok(_) => {}
                Err(e) => error!("Cannot write the checkpoint: {}", e),
            }
        }
        if is_last {
            IS_INSERTION_COMPLETE.store(true, Ordering::SeqCst);
        }
        if is_stop {
            std::process::exit(STOPPED_EXIT_CODE);
        }
    };

    let network = match &cli.resume {
        None => generate(input, family_collection, &mut on_progress)?,
        Some(checkpoint_file) => {
            let checkpoint = read_checkpoint(checkpoint_file)?;
            let mode = if cli.finish {
                ResumeMode::Finish
            } else {
                ResumeMode::ContinueInsertion
            };
            let restart = Resume {
                checkpoint: &checkpoint,
                mode,
            };
            resume(input, family_collection, restart, &mut on_progress)?
        }
    };

    write_output_to_folder(&network, &cli.output_folder)?;
    remove_checkpoint_files(&cli.output_folder)?;

    Ok(())
}

/// Logs to stdout. The level is `info`, unless the environment variable `RUST_LOG` gives another.
fn init_tracing() {
    let level_filter = tracing_subscriber::EnvFilter::try_from_default_env()
        .unwrap_or_else(|_| tracing_subscriber::EnvFilter::new("info"));

    tracing_subscriber::registry()
        .with(level_filter)
        .with(
            tracing_subscriber::fmt::layer()
                .compact()
                .with_writer(|| LogWriter),
        )
        .init();
}
