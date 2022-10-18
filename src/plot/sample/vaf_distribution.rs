//! GC content distribution plot for a single sample.

use anyhow::bail;
use plotly::common::Line;
use plotly::common::LineShape;
use plotly::common::Mode;
use plotly::common::Title;
use plotly::layout::Axis;
use plotly::Layout;
use plotly::Scatter;

use crate::plot::command::FilepathResults;
use crate::plot::command::SamplePlot;

/// Struct that represents a GC content distribution plot for a single sample.
pub struct VariantAlleleFractionDistributionPlot;

impl SamplePlot for VariantAlleleFractionDistributionPlot {
    fn name(&self) -> &'static str {
        "Variant Allele Fraction Distribution"
    }

    fn description(&self) -> &'static str {
        "Shows the probability distribution from which VAFs are pulled."
    }

    fn filename(&self) -> &'static str {
        "vaf-distribution"
    }

    fn generate(
        &self,
        filepath_results: &FilepathResults,
        title: Title,
    ) -> anyhow::Result<plotly::Plot> {
        let mut plot = plotly::Plot::new();
        let FilepathResults(filepath, results) = filepath_results;

        // (1) Check to make sure that the results files have the necessary
        // keys to plot the data. If they don't then we need to fail as there
        // will be nothing to plot.
        let edits = match &results.edits {
            Some(e) => e,
            None => bail!("File {} has no Edits information!", filepath.display()),
        };

        // (2) Check to see if the Edits information is empty. This just means
        // whether or not the histogram has collected any information. If the
        // sum of the histogram is zero, then clearly the Edits did not complete
        // correctly. As such, we cannot process this sample get a divide by
        // zero error later).
        let histogram = &edits.vaf_histogram;
        let total = histogram.sum();

        if total == 0 {
            bail!(
                "File {} has  information, but it's empty!",
                filepath.display()
            );
        }

        // (3) Prepare the X and Y points for this trace (line). We're going
        // to normalize the values so that they plot on the same scale by
        // dividing each value by the sum for this histogram.
        let mut x = Vec::new();
        let mut y = Vec::new();

        let histogram = histogram.values();
        for (position, value) in histogram.iter().enumerate() {
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
            .mode(Mode::LinesMarkers)
            .name(trace_name)
            .line(Line::new().shape(LineShape::Spline));
        plot.add_trace(trace);

        // (6) Configure the graph for plotting and return.
        let layout = Layout::new()
            .title(title)
            .x_axis(
                Axis::new()
                    .title(Title::new("Variant Allele Fraction"))
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
