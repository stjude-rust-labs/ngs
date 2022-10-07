pub mod gc_content_distribution;
pub mod quality_score_distribution;

use std::path::PathBuf;

use anyhow::Context;
use clap::Args;
use tracing::{debug, info};

use crate::{
    plot::command::{get_all_sample_plots, FilepathResults},
    qc::results::Results,
};

#[derive(Args)]
pub struct SampleArgs {
    /// `ngs qc` results file for which to generate the plot(s).
    #[arg(value_name = "JSON")]
    pub src: PathBuf,

    /// The directory to output all files within.
    #[arg(short, long, value_name = "PATH")]
    pub output_directory: Option<PathBuf>,
}

pub fn plot(args: SampleArgs) -> anyhow::Result<()> {
    //========//
    // Source //
    //========//

    let results = Results::read(&args.src)
        .with_context(|| format!("invalid input file: {}", args.src.display()))?;

    let filepath_results = FilepathResults(args.src.clone(), results);
    debug!("  [*] Source: {}", args.src.display());

    //==================//
    // Output Directory //
    //==================//

    let output_directory = if let Some(o) = args.output_directory {
        o
    } else {
        std::env::current_dir().expect("could not retrieve the current working directory.")
    };
    debug!("  [*] Output directory: {}", output_directory.display());

    let plots = get_all_sample_plots();
    for p in plots {
        let plot = p.generate(&filepath_results)?;

        let mut filename = output_directory.clone();
        filename.push(String::from(p.filename()) + ".sample.html");

        info!("  [*] Writing {} to {}", p.name(), filename.display());
        plot.write_html(filename);
    }

    Ok(())
}
