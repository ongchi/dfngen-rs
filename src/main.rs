use std::process::ExitCode;
use std::time::{SystemTime, UNIX_EPOCH};

use clap::Parser;
use dfngen::{generate, read_input, write_output_to_folder, DfngenError};
use tracing::{error, info};
use tracing_subscriber::{layer::SubscriberExt, util::SubscriberInitExt};

use crate::progress::{GenerationProgress, LogWriter};

mod progress;

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    /// Path to input file
    input_file: String,

    /// Path to output folder
    output_folder: String,
}

fn main() -> ExitCode {
    let cli = Cli::parse();

    init_tracing();

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

/// Reads the input, generates the network and writes the output files.
fn run(cli: &Cli) -> Result<(), DfngenError> {
    let (mut input, family_collection) = read_input(&cli.input_file)?;

    if input.seed == 0 {
        input.seed = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .map(|t| t.as_secs())?;
    }

    // The bar starts with the first accepted stochastic fracture
    let stop_condition = input.stop_condition;
    let mut progress = None;
    let network = generate(input, family_collection, &mut |update| {
        progress
            .get_or_insert_with(|| GenerationProgress::new(stop_condition))
            .update(update);
    })?;

    write_output_to_folder(&network, &cli.output_folder)?;

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
