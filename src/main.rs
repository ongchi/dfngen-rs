use std::time::{SystemTime, UNIX_EPOCH};

use clap::Parser;
use dfngen::io::output::write_output;
use dfngen::{generate, read_input, DfngenError};
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

fn main() -> Result<(), DfngenError> {
    let cli = Cli::parse();

    init_tracing();

    info!("Starting dfngen-rs");

    let (mut input, frac_fam_col) = read_input(&cli.input_file).inspect_err(|e| error!("{}", e))?;

    if input.seed == 0 {
        input.seed = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .map(|t| t.as_secs())?;
    }

    // The bar starts with the first accepted stochastic fracture
    let (stop_condition, n_poly) = (input.stop_condition, input.n_poly);
    let mut progress = None;
    let mut network = generate(input, frac_fam_col, &mut |update| {
        progress
            .get_or_insert_with(|| GenerationProgress::new(stop_condition, n_poly))
            .update(update)
    })?;

    // Write all output files
    write_output(
        &network.input,
        &cli.output_folder,
        &network.dfn.accepted_poly,
        &network.dfn.rejected_poly,
        &mut network.dfn.intpts,
        &network.dfn.triple_points,
        &network.dfn.pstats,
        &network.final_fractures,
    )?;

    info!("dfngen-rs - Complete");

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
