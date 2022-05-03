use clap::{Arg, Command};
use std::io;

mod commands;
mod utils;

use git_testament::{git_testament, render_testament};
use tracing_subscriber;

git_testament!(TESTAMENT);

fn main() -> io::Result<()> {
    let version = render_testament!(TESTAMENT);

    let flagstat_cmd = Command::new("flagstat")
        .about("Generates stats from the marked flags in a SAM/BAM/CRAM file.")
        .arg(Arg::new("src").help("Source file.").index(1).required(true));

    let matches = Command::new("ngs")
        .version(version.as_str())
        .propagate_version(true)
        .subcommand_required(true)
        .subcommand(flagstat_cmd)
        .get_matches();

    let subscriber = tracing_subscriber::fmt::Subscriber::builder()
        .with_max_level(tracing::Level::DEBUG)
        .with_writer(std::io::stderr)
        .finish();
    let _ = tracing::subscriber::set_global_default(subscriber);

    if let Some(m) = matches.subcommand_matches("flagstat") {
        commands::flagstat(m)
    } else {
        unreachable!();
    }
}
