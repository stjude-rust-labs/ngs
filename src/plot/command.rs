use std::path::PathBuf;

use clap::{command, Args, Subcommand};

use crate::qc::results::Results;

use super::{
    cohort::{self, CohortArgs},
    sample::{self, SampleArgs},
};

//===============//
// Command setup //
//===============//

#[derive(Args)]
pub struct PlotArgs {
    #[command(subcommand)]
    pub subcommand: PlotSubcommand,
}

#[derive(Subcommand)]
pub enum PlotSubcommand {
    /// Plots cohort-level information produced by the `ngs qc` command.
    Cohort(CohortArgs),
    /// Plots sample-specific information produced by the `ngs qc` command.
    Sample(SampleArgs),
}

pub struct FilepathResults(pub PathBuf, pub Results);

//===================//
// Sample Plot trait //
//===================//

pub trait SamplePlot {
    fn name(&self) -> &'static str;
    fn filename(&self) -> &'static str;
    fn generate(&self, filepath_results: &FilepathResults) -> anyhow::Result<plotly::Plot>;
}

//===================//
// Cohort Plot trait //
//===================//

pub trait CohortPlot {
    fn name(&self) -> &'static str;
    fn filename(&self) -> &'static str;
    fn generate(&self, filepath_results: &[FilepathResults]) -> anyhow::Result<plotly::Plot>;
}

//================================//
// Get all supported sample plots //
//================================//

pub fn get_all_sample_plots() -> Vec<Box<dyn SamplePlot>> {
    vec![
        Box::new(sample::quality_score_distribution::QualityScoreDistributionPlot),
        Box::new(sample::gc_content_distribution::GCContentDistributionPlot),
    ]
}

//================================//
// Get all supported cohort plots //
//================================//

pub fn get_all_cohort_plots() -> Vec<Box<dyn CohortPlot>> {
    vec![Box::new(
        cohort::gc_content_distribution::GCContentDistributionPlot,
    )]
}
