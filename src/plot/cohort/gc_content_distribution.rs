//! GC content distribution plot for a cohort of samples.

use plotly::{
    common::{Line, LineShape, Mode, Title},
    layout::Axis,
    Layout, Scatter,
};
use tracing::error;

use crate::plot::command::{CohortPlot, FilepathResults};

/// Struct that represents a GC content distribution plot for a cohort.
pub struct GCContentDistributionPlot;

impl CohortPlot for GCContentDistributionPlot {
    fn name(&self) -> &'static str {
        "GC Content Distribution"
    }

    fn filename(&self) -> &'static str {
        "gc-content-distribution"
    }

    fn generate(&self, filepath_results: &[FilepathResults]) -> anyhow::Result<plotly::Plot> {
        let mut plot = plotly::Plot::new();

        for FilepathResults(filepath, result) in filepath_results {
            // (1) Check to make sure that the results files have the necessary
            // keys to plot the data. If they don't then we can ignore that
            // result and continue on.
            let gc_content = match &result.gc_content {
                Some(gc) => gc,
                None => {
                    error!(
                        "  [*] File {} has no GC content information! Skipping...",
                        filepath.display()
                    );
                    continue;
                }
            };

            // (2) Check to see if the GC content information is empty. This
            // just means whether or not the histogram has collected any
            // information. Often for files with shorter reads, all of the reads
            // will be rejected based on being too short, so the field will be
            // present, but the values will all be zero. If the sum of the
            // histogram is zero, then clearly the GC content did not complete
            // correctly. As such, we cannot process this sample (or else we'll
            // get a divide by zero error later).
            let histogram = &gc_content.histogram;
            let total = histogram.sum();

            if total == 0 {
                error!(
                    "  [*] File {} has GC content information, but it's empty! Skipping...",
                    filepath.display()
                );
                continue;
            }

            // (3) Prepare the X and Y points for this trace (line). We're going
            // to normalize the values so that they plot on the same scale by
            // dividing each value by the sum for this histogram.
            let mut x = Vec::new();
            let mut y = Vec::new();

            for (position, value) in histogram.values().iter().enumerate() {
                x.push(position);
                y.push(*value as f64 / total as f64);
            }

            // (4) Generalize the tracename based on the filename.
            // TODO: This could be done quite a bit better.
            let trace_name = filepath
                .file_name()
                .expect("the file to have a filename")
                .to_os_string()
                .into_string()
                .expect("the filename to be convertable to a string")
                .replace(".results.json", "");

            // (5) Add this trace to the plot.
            let trace = Scatter::new(x, y)
                .mode(Mode::Lines)
                .name(trace_name)
                .line(Line::new().shape(LineShape::Spline));

            plot.add_trace(trace);
        }

        // (6) Configure the graph for plotting and return.
        let layout = Layout::new()
            .title(Title::new(self.name()))
            .x_axis(
                Axis::new()
                    .title(Title::new("Percentage GC"))
                    .auto_range(true),
            )
            .y_axis(
                Axis::new()
                    .title(Title::new("Fraction of Records"))
                    .auto_range(true),
            );

        plot.set_layout(layout);
        Ok(plot)
    }
}
