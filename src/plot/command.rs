//! Functionality related to the `ngs plot` command itself.

use std::path::PathBuf;

use anyhow::bail;
use clap::{command, Args, Subcommand};
use itertools::Itertools;

use crate::qc::results::Results;

use super::{
    cohort::{self, PlotCohortArgs},
    sample::{self, PlotSampleArgs},
};

//===============//
// Command setup //
//===============//

/// Command line arguments for `ngs plot`.
#[derive(Args)]
pub struct PlotArgs {
    /// The subcommand for `ngs plot`.
    #[command(subcommand)]
    pub subcommand: PlotSubcommand,
}

/// All possible subcommands for `ngs plot`.
#[derive(Subcommand)]
pub enum PlotSubcommand {
    /// Plots cohort-level information produced by the `ngs qc` command.
    Cohort(PlotCohortArgs),

    /// Plots sample-specific information produced by the `ngs qc` command.
    Sample(PlotSampleArgs),
}

/// Utility struct for passing around a combination of the filepath and
/// the loaded `ngs qc` results. This is useful when passing around this tuple
/// of arguments to various functions in the `ngs plot` subcommand(s).
pub struct FilepathResults(pub PathBuf, pub Results);

//===================//
// Sample Plot trait //
//===================//

/// A single-sample graph that can be plotted by `ngs plot`.
pub trait SamplePlot {
    /// The name of this plot.
    fn name(&self) -> &'static str;

    /// The filename to output for this plot.
    fn filename(&self) -> &'static str;

    /// Generates the plot given a loaded results file.
    fn generate(&self, filepath_results: &FilepathResults) -> anyhow::Result<plotly::Plot>;
}

//===================//
// Cohort Plot trait //
//===================//

/// A cohort-level graph that can be plotted by `ngs plot`.
pub trait CohortPlot {
    /// The name of this plot.
    fn name(&self) -> &'static str;

    /// The filename to output for this plot.
    fn filename(&self) -> &'static str;

    /// Generates the plot given a set of loaded results files.
    fn generate(&self, filepath_results: &[FilepathResults]) -> anyhow::Result<plotly::Plot>;
}

//==============//
// Sample plots //
//==============//

/// Gets all of the supported single-sample plots.
pub fn get_all_sample_plots(
    only_graph: Option<String>,
) -> anyhow::Result<Vec<Box<dyn SamplePlot>>> {
    let mut results: Vec<Box<dyn SamplePlot>> = vec![
        Box::new(sample::quality_score_distribution::QualityScoreDistributionPlot),
        Box::new(sample::gc_content_distribution::GCContentDistributionPlot),
    ];

    if let Some(only) = only_graph {
        results = results
            .into_iter()
            .filter(|x| x.name().eq_ignore_ascii_case(&only))
            .collect_vec();

        if results.is_empty() {
            bail!("No plots matched the specified `--only` flag: {}", only);
        }
    }

    Ok(results)
}

//==============//
// Cohort plots //
//==============//

/// Gets all of the supported cohort-level plots.
pub fn get_all_cohort_plots(
    only_graph: Option<String>,
) -> anyhow::Result<Vec<Box<dyn CohortPlot>>> {
    let mut results: Vec<Box<dyn CohortPlot>> = vec![Box::new(
        cohort::gc_content_distribution::GCContentDistributionPlot,
    )];

    if let Some(only) = only_graph {
        results = results
            .into_iter()
            .filter(|x| x.name().eq_ignore_ascii_case(&only))
            .collect_vec();

        if results.is_empty() {
            bail!("No plots matched the specified `--only` flag: {}", only);
        }
    }

    Ok(results)
}
