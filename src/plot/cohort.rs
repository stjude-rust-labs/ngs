//! Functionality related to the `ngs plot cohort` subcommand.

pub mod gc_content_distribution;

use anyhow::Context;
use std::path::PathBuf;
use tracing::{debug, info};

use clap::Args;

use crate::{
    plot::command::{get_all_cohort_plots, FilepathResults},
    qc::results::Results,
};

/// Clap arguments for the `ngs plot cohort` subcommand.
#[derive(Args)]
pub struct PlotCohortArgs {
    /// `ngs qc` results files for which to generate the plot(s).
    #[arg(required = true, value_name = "JSON")] // required implies one or more
    pub src: Vec<PathBuf>,

    /// The directory to output all files within.
    #[arg(short, long, value_name = "PATH")]
    pub output_directory: Option<PathBuf>,
}

/// Main method for the `ngs plot cohort` subcommand.
pub fn plot(args: PlotCohortArgs) -> anyhow::Result<()> {
    //========//
    // Source //
    //========//

    let mut filepath_results = Vec::new();

    for src in args.src {
        let results = Results::read(&src)
            .with_context(|| format!("invalid input file: {}", src.display()))?;
        filepath_results.push(FilepathResults(src.clone(), results));
        debug!("  [*] Source: {}", src.display());
    }

    //==================//
    // Output Directory //
    //==================//

    let output_directory = if let Some(o) = args.output_directory {
        o
    } else {
        std::env::current_dir().expect("could not retrieve the current working directory.")
    };
    debug!("  [*] Output directory: {}", output_directory.display());

    //================//
    // Generate Plots //
    //================//

    let plots = get_all_cohort_plots();
    for p in plots {
        let plot = p.generate(&filepath_results)?;

        let mut filename = output_directory.clone();
        filename.push(String::from(p.filename()) + ".cohort.html");

        info!("  [*] Writing {} to {}", p.name(), filename.display());
        plot.write_html(filename);
    }

    Ok(())
}
