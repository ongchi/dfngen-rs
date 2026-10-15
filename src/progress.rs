use std::io::{self, Write};
use std::sync::LazyLock;

use indicatif::{MultiProgress, ProgressBar, ProgressFinish, ProgressStyle};

use dfngen::io::input::StopCondition;
use dfngen::ProgressUpdate;

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
    pub fn new(stop_condition: StopCondition) -> Self {
        let (length, template) = match stop_condition {
            StopCondition::NPoly(n_poly) => (
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
    /// The bar is removed from the terminal when the stop condition is complete.
    pub fn update(&self, update: &ProgressUpdate) {
        match self.stop_condition {
            StopCondition::NPoly(_) => self.bar.set_position(update.accepted_count as u64),
            StopCondition::P32 => {
                self.bar.set_message(update.accepted_count.to_string());
                self.bar
                    .set_position((update.fraction * P32_BAR_STEPS as f64) as u64);
            }
        }

        if update.fraction >= 1. {
            self.bar.finish_and_clear();
        }
    }
}
