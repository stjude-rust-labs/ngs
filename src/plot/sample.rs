//! Functionality related to the `ngs plot sample` subcommand.

pub mod gc_content_distribution;
pub mod quality_score_distribution;
pub mod vaf_distribution;

use std::path::PathBuf;

use anyhow::Context;
use clap::Args;
use plotly::common::Title;
use tracing::debug;
use tracing::info;

use crate::plot::command::get_all_sample_plots;
use crate::plot::command::FilepathResults;
use crate::qc::results::Results;

/// Clap arguments for the `ngs plot sample` subcommand.
#[derive(Args)]
pub struct PlotSampleArgs {
    /// `ngs qc` results file for which to generate the plot(s).
    #[arg(value_name = "JSON")]
    pub src: PathBuf,

    /// The directory to output all files within.
    #[arg(short, long, value_name = "PATH")]
    pub output_directory: Option<PathBuf>,

    /// If provided, the name of the sample being processed.
    #[arg(short, long = "sample")]
    pub sample_name: Option<String>,

    /// If provided, only prepares the plot specified (by plot name).
    #[arg(long = "only")]
    pub only_graph: Option<String>,
}

/// Main method for the `ngs plot sample` subcommand.
pub fn plot(args: PlotSampleArgs) -> anyhow::Result<()> {
    //========//
    // Source //
    //========//

    let results = Results::read(&args.src)
        .with_context(|| format!("invalid input file: {}", args.src.display()))?;

    let filepath_results = FilepathResults(args.src.clone(), results);
    debug!("  [*] Source: {}", args.src.display());

    //=============//
    // Sample name //
    //=============//

    let sample_name = args.sample_name;

    //==================//
    // Output Directory //
    //==================//

    let output_directory = if let Some(o) = args.output_directory {
        o
    } else {
        std::env::current_dir().expect("could not retrieve the current working directory.")
    };
    debug!("  [*] Output directory: {}", output_directory.display());

    let plots = get_all_sample_plots(args.only_graph)?;
    for p in plots {
        // (1) Determine what the title should be.
        let title = match sample_name {
            Some(ref sample) => Title::new(&format!("{} - {}", p.name(), sample)),
            None => Title::new(p.name()),
        };

        // (2) Generate the plot in question.
        let plot = p.generate(&filepath_results, title)?;

        // (3) Write the plot to the appropriate output file.
        let mut filename = output_directory.clone();
        filename.push(String::from(p.filename()) + ".sample.html");

        info!("  [*] Writing {} to {}", p.name(), filename.display());
        plot.write_html(filename);
    }

    Ok(())
}
