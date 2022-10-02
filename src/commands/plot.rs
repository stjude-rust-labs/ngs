pub mod cohort;
pub mod sample;

use std::path::PathBuf;

use clap::Command;

use crate::{add_verbosity_args, lib::qc::results::Results};

//===============//
// Command setup //
//===============//

pub fn get_command<'a>() -> Command<'a> {
    Command::new("plot")
        .about("Produces plots for data generated by `ngs qc`.")
        .subcommand_required(true)
        .subcommand(add_verbosity_args(cohort::get_command()))
        .subcommand(add_verbosity_args(sample::get_command()))
}

pub struct FilepathResults<'a>(&'a PathBuf, Results);

//===================//
// Sample Plot trait //
//===================//

pub trait SamplePlot {
    fn name(&self) -> &'static str;
    fn filename(&self) -> &'static str;
    fn generate(&self, filepath_results: &FilepathResults<'_>) -> anyhow::Result<plotly::Plot>;
}

//===================//
// Cohort Plot trait //
//===================//

pub trait CohortPlot {
    fn name(&self) -> &'static str;
    fn filename(&self) -> &'static str;
    fn generate(&self, filepath_results: &[FilepathResults<'_>]) -> anyhow::Result<plotly::Plot>;
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
