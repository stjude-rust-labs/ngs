use anyhow::bail;
use itertools::Itertools;
use plotly::{
    common::{ErrorData, ErrorType, Line, LineShape, Mode, Title},
    layout::Axis,
    Layout, Scatter,
};

use crate::commands::plot::{FilepathResults, SamplePlot};

pub struct QualityScoreDistributionPlot;

impl SamplePlot for QualityScoreDistributionPlot {
    fn name(&self) -> &'static str {
        "Quality Score Distribution"
    }

    fn filename(&self) -> &'static str {
        "quality-score-distribution"
    }

    fn generate(&self, filepath_results: &FilepathResults<'_>) -> anyhow::Result<plotly::Plot> {
        let mut plot = plotly::Plot::new();
        let FilepathResults(filepath, results) = filepath_results;

        // (1) Check to make sure that the results files have the necessary
        // keys to plot the data. If they don't then we need to fail as there
        // will be nothing to plot.
        let quality_scores = match results.quality_scores() {
            Some(qs) => qs,
            None => bail!(
                "File {} has no quality score information!",
                filepath.display()
            ),
        };

        // (2) Tally up the quality scores and configure the plot.
        let scores = &quality_scores.scores;
        let mut x = Vec::new();
        let mut y = Vec::new();
        let mut error_plus = Vec::new();
        let mut error_minus = Vec::new();
        let mut y_lim = 0.0;

        for position in scores.keys().sorted() {
            let bin = scores.get(position).unwrap();
            let median = bin.median().unwrap();
            let third_quartile = bin.third_quartile().unwrap();
            let first_quartile = bin.first_quartile().unwrap();

            x.push(*position);
            y.push(median);
            error_plus.push(third_quartile - median);
            error_minus.push(median - first_quartile);

            if third_quartile > y_lim {
                y_lim = third_quartile;
            }
        }

        // (3) Generalize the tracename based on the filename.
        // TODO: This could be done quite a bit better.
        let trace_name = filepath
            .file_name()
            .expect("the file to have a filename")
            .to_os_string()
            .into_string()
            .expect("the filename to be convertable to a string")
            .replace(".results.json", "");

        // (4) Add this trace to the plot.
        let trace = Scatter::new(x, y)
            .mode(Mode::LinesMarkers)
            .name(trace_name)
            .line(Line::new().shape(LineShape::Spline))
            .error_y(
                ErrorData::new(ErrorType::Data)
                    .array(error_plus)
                    .array_minus(error_minus),
            );
        plot.add_trace(trace);

        // (5) Configure the graph for plotting and return.
        let layout = Layout::new()
            .title(Title::new(self.name()))
            .x_axis(Axis::new().title(Title::new("Position")).auto_range(true))
            .y_axis(
                Axis::new()
                    .title(Title::new("Quality Score"))
                    .range(vec![0, y_lim as usize + 1]),
            );
        plot.set_layout(layout);

        Ok(plot)
    }
}
