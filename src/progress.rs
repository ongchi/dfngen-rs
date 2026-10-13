use std::io::{self, Write};
use std::sync::LazyLock;

use indicatif::{MultiProgress, ProgressBar, ProgressFinish, ProgressStyle};

use dfngen::fracture::fracture_family::FractureFamily;
use dfngen::io::input::StopCondition;
use dfngen::math_functions::p32_progress_fraction;

/// Number of bar steps for the P32 stop condition, where the number of fractures is not known
const P32_BAR_STEPS: u64 = 1000;

/// Draws the progress bars on stderr, only if stderr is a terminal.
/// The log writer and the progress bars use it together.
static TERMINAL_PROGRESS: LazyLock<MultiProgress> = LazyLock::new(MultiProgress::new);

/// Writes the log to stdout
///
/// The progress bar is hidden while a log line is written. Thus a log line does
/// not break the bar when stdout and stderr are the same terminal.
pub struct LogWriter;

impl Write for LogWriter {
    fn write(&mut self, buffer: &[u8]) -> io::Result<usize> {
        TERMINAL_PROGRESS.suspend(|| io::stdout().write(buffer))
    }

    fn flush(&mut self) -> io::Result<()> {
        io::stdout().flush()
    }
}

/// Progress bar for the insertion of the stochastic fractures
///
/// The bar is removed from the terminal when this value is dropped.
pub struct GenerationProgress {
    bar: ProgressBar,
    stop_condition: StopCondition,
}

impl GenerationProgress {
    /// # Arguments
    ///
    /// * `stop_condition` - Stop condition of the generation
    /// * `n_poly` - Number of fractures to accept with the `NPoly` stop condition
    pub fn new(stop_condition: StopCondition, n_poly: usize) -> Self {
        let (length, template) = match stop_condition {
            StopCondition::NPoly => (
                n_poly as u64,
                "[{elapsed_precise}] [{bar:40}] {percent:>3} %  {pos}/{len} fractures  ETA {eta_precise}",
            ),
            StopCondition::P32 => (
                P32_BAR_STEPS,
                "[{elapsed_precise}] [{bar:40}] {percent:>3} %  {msg} fractures  ETA {eta_precise}",
            ),
        };
        let style = ProgressStyle::with_template(template)
            .expect("the template is a constant and is valid")
            .progress_chars("#>-");
        let bar = ProgressBar::new(length)
            .with_style(style)
            .with_message("0")
            .with_finish(ProgressFinish::AndClear);

        Self {
            bar: TERMINAL_PROGRESS.add(bar),
            stop_condition,
        }
    }

    /// Shows the progress after a fracture was accepted
    ///
    /// # Arguments
    ///
    /// * `accepted_count` - Number of accepted fractures
    /// * `families` - Stochastic fracture families, with their current and target P32
    pub fn update(&self, accepted_count: usize, families: &[FractureFamily]) {
        match self.stop_condition {
            StopCondition::NPoly => self.bar.set_position(accepted_count as u64),
            StopCondition::P32 => {
                let current_and_target_p32: Vec<(f64, f64)> = families
                    .iter()
                    .map(|family| (family.current_p32, family.p32_target))
                    .collect();
                let fraction = p32_progress_fraction(&current_and_target_p32);

                self.bar.set_message(accepted_count.to_string());
                self.bar
                    .set_position((fraction * P32_BAR_STEPS as f64) as u64);
            }
        }
    }
}
