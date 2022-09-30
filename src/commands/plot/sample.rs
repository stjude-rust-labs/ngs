pub mod quality_score_distribution;

use anyhow::Context;
use clap::{Arg, ArgMatches, Command};

use crate::lib::qc::results::Results;

use super::get_all_plots;

pub fn get_command<'a>() -> Command<'a> {
    Command::new("sample")
        .about("Plots sample-specific information produced by the `ngs qc` command.")
        .arg(
            Arg::new("src")
                .help("`ngs qc` results file for which to generate the plot(s)")
                .required(true),
        )
}

pub fn plot(matches: &ArgMatches) -> anyhow::Result<()> {
    let src: &String = matches.get_one("src").expect("missing src filepath");
    let results = Results::read(src).with_context(|| "invalid input file")?;
    for p in get_all_plots() {
        let plot = p.generate(&results)?;
        plot.show();
    }

    Ok(())
}
