use clap::{arg, Command};
use std::io;

mod commands;
mod common;
mod derive;
mod errors;
mod generate;
mod qc;
mod utils;

use git_testament::{git_testament, render_testament};

git_testament!(TESTAMENT);

/// Utility method to add the verbosity arguments to any subcommands passed to Clap.
///
/// # Arguments
///
/// * `subcommand` — The Clap subcommand to add these arguments to.
fn add_verbosity_args(subcommand: Command) -> Command {
    subcommand
        .arg(arg!(-q --quiet "Only errors are printed to the stderr stream."))
        .arg(
            arg!(-v --verbose "All available information, including debug information, is \
                printed to stderr."),
        )
}

fn main() -> io::Result<()> {
    let version = render_testament!(TESTAMENT);

    let derive_cmd = commands::derive::get_command();
    let generate_cmd = commands::generate::get_command();
    let qc_cmd = commands::qc::get_command();

    let matches = Command::new("ngs")
        .version(version.as_str())
        .propagate_version(true)
        .subcommand_required(true)
        .subcommand(add_verbosity_args(derive_cmd))
        .subcommand(add_verbosity_args(generate_cmd))
        .subcommand(add_verbosity_args(qc_cmd))
        .get_matches();

    if let Some((name, subcommand)) = matches.subcommand() {
        let mut level = tracing::Level::INFO;
        if subcommand.is_present("quiet") {
            level = tracing::Level::ERROR;
        } else if subcommand.is_present("verbose") {
            level = tracing::Level::DEBUG;
        }

        let subscriber = tracing_subscriber::fmt::Subscriber::builder()
            .with_max_level(level)
            .with_writer(std::io::stderr)
            .finish();
        let _ = tracing::subscriber::set_global_default(subscriber);

        match name {
            "derive" => {
                if let Some(m) = subcommand.subcommand_matches("instrument") {
                    return commands::derive::instrument::derive(m);
                } else {
                    unreachable!();
                }
            }
            "generate" => return commands::generate(subcommand),
            "qc" => return commands::qc(subcommand),
            s => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("Unknown subcommand: {}", s),
                ))
            }
        }
    }

    unreachable!();
}
