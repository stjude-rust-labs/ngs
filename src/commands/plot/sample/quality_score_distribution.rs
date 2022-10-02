use anyhow::bail;
use itertools::Itertools;
use plotly::{
    common::{ErrorData, ErrorType, Line, LineShape, Mode, Title},
    layout::Axis,
    Layout, Scatter,
};

use crate::commands::plot::{Plot, PlotKind};

pub struct QualityScoreDistributionPlot;

impl Plot for QualityScoreDistributionPlot {
    fn name(&self) -> &'static str {
        "Quality Score Distribution"
    }

    fn kind(&self) -> crate::commands::plot::PlotKind {
        PlotKind::Sample
    }

    fn filename(&self) -> &'static str {
        "quality-score-distribution"
    }

    fn generate(&self, data: &crate::lib::qc::results::Results) -> anyhow::Result<plotly::Plot> {
        let mut plot = plotly::Plot::new();

        let quality_scores = match data.quality_scores() {
            Some(qs) => qs,
            None => bail!("Quality scores are required to plot {}", self.name()),
        };

        let scores = quality_scores.scores();
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

        let layout = Layout::new()
            .title(Title::new(self.name()))
            .x_axis(
                Axis::new()
                    .title(Title::new("Position"))
                    .range(vec![0, x.last().unwrap() + 1]),
            )
            .y_axis(
                Axis::new()
                    .title(Title::new("Quality Score"))
                    .range(vec![0, y_lim as usize + 1]),
            );

        let trace = Scatter::new(x, y)
            .mode(Mode::LinesMarkers)
            .name("Quality Scores")
            .line(Line::new().shape(LineShape::Spline))
            .error_y(
                ErrorData::new(ErrorType::Data)
                    .array(error_plus)
                    .array_minus(error_minus),
            );

        plot.set_layout(layout);
        plot.add_trace(trace);

        Ok(plot)
    }
}
