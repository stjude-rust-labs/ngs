pub mod gc_content_distribution;

use anyhow::Context;
use std::path::PathBuf;
use tracing::{debug, info};

use clap::{value_parser, Arg, ArgMatches, Command};

use crate::{
    commands::plot::{get_all_cohort_plots, FilepathResults},
    lib::qc::results::Results,
};

pub fn get_command<'a>() -> Command<'a> {
    Command::new("cohort")
        .about("Plots cohort-level information produced by the `ngs qc` command.")
        .arg(
            Arg::new("src")
                .help("`ngs qc` results files for which to generate the plot(s)")
                .value_parser(value_parser!(PathBuf))
                .multiple_values(true)
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

    let sources: Vec<&PathBuf> = matches
        .get_many("src")
        .expect("missing src filepath(s)")
        .collect();

    let mut results = Vec::new();

    for src in sources {
        results.push(FilepathResults(
            src,
            Results::read(src).with_context(|| format!("invalid input file: {}", src.display()))?,
        ));
        debug!("  [*] Source: {}", src.display());
    }

    //==================//
    // Output Directory //
    //==================//

    let output_directory = if let Some(m) = matches.get_one::<PathBuf>("output-directory") {
        PathBuf::from(m)
    } else {
        std::env::current_dir().expect("Could not retrieve the current working directory.")
    };
    debug!("  [*] Output directory: {}", output_directory.display());

    //================//
    // Generate Plots //
    //================//

    let plots = get_all_cohort_plots();
    for p in plots {
        let plot = p.generate(&results)?;

        let mut filename = output_directory.clone();
        filename.push(String::from(p.filename()) + ".cohort.html");

        info!("  [*] Writing {} to {}", p.name(), filename.display());
        plot.write_html(filename);
    }

    Ok(())
}
