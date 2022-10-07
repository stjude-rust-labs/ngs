pub mod gc_content_distribution;
pub mod quality_score_distribution;

use std::path::PathBuf;

use anyhow::Context;
use clap::{value_parser, Arg, ArgMatches, Command};
use tracing::{debug, info};

use crate::{
    plot::command::{get_all_sample_plots, FilepathResults},
    qc::results::Results,
};

pub fn get_command<'a>() -> Command<'a> {
    Command::new("sample")
        .about("Plots sample-specific information produced by the `ngs qc` command.")
        .arg(
            Arg::new("src")
                .help("`ngs qc` results file for which to generate the plot(s)")
                .value_parser(value_parser!(PathBuf))
                .required(true),
        )
        .arg(
            Arg::new("output-directory")
                .long("--output-directory")
                .short('o')
                .help("The directory to output files to.")
                .value_parser(value_parser!(PathBuf))
                .required(false)
                .takes_value(true),
        )
}

pub fn plot(matches: &ArgMatches) -> anyhow::Result<()> {
    //========//
    // Source //
    //========//

    let src: &PathBuf = matches.get_one("src").expect("missing src filepath");
    let results = FilepathResults(
        src,
        Results::read(src).with_context(|| format!("invalid input file: {}", src.display()))?,
    );
    debug!("  [*] Source: {}", src.display());

    //==================//
    // Output Directory //
    //==================//

    let output_directory = if let Some(m) = matches.get_one::<PathBuf>("output-directory") {
        PathBuf::from(m)
    } else {
        std::env::current_dir().expect("Could not retrieve the current working directory.")
    };
    debug!("  [*] Output directory: {}", output_directory.display());

    let plots = get_all_sample_plots();
    for p in plots {
        let plot = p.generate(&results)?;

        let mut filename = output_directory.clone();
        filename.push(String::from(p.filename()) + ".sample.html");

        info!("  [*] Writing {} to {}", p.name(), filename.display());
        plot.write_html(filename);
    }

    Ok(())
}
