use anyhow::bail;
use plotly::{
    common::{ErrorData, ErrorType, Line, LineShape, Mode, Title},
    layout::Axis,
    Layout, Scatter,
};

use crate::commands::plot::{Plot, PlotKind};

pub struct GCContentDistributionPlot;

impl Plot for GCContentDistributionPlot {
    fn name(&self) -> &'static str {
        "GC Content Distribution"
    }

    fn kind(&self) -> crate::commands::plot::PlotKind {
        PlotKind::Sample
    }

    fn filename(&self) -> &'static str {
        "gc-content-distribution"
    }

    fn generate(&self, data: &crate::lib::qc::results::Results) -> anyhow::Result<plotly::Plot> {
        let mut plot = plotly::Plot::new();

        let gc_content = match data.gc_content() {
            Some(gc) => gc,
            None => bail!("GC Content Metrics are required to plot {}", self.name()),
        };

        let mut x = Vec::new();
        let mut y = Vec::new();

        let histogram = gc_content.histogram.values();
        for (position, value) in histogram.iter().enumerate() {
            x.push(position);
            y.push(*value);
        }

        let layout = Layout::new()
            .title(Title::new(self.name()))
            .x_axis(
                Axis::new()
                    .title(Title::new("Percentage GC"))
                    .range(vec![0, *x.last().unwrap()]),
            )
            .y_axis(
                Axis::new()
                    .title(Title::new("Number of Records"))
                    .auto_range(true),
            );

        let trace = Scatter::new(x, y)
            .mode(Mode::LinesMarkers)
            .name("GC Content Percentage")
            .line(Line::new().shape(LineShape::Spline));

        plot.set_layout(layout);
        plot.add_trace(trace);

        Ok(plot)
    }
}
