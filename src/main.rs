use clap::{arg, Arg, ArgGroup, Command};
use std::io;

mod commands;
mod derive;
mod exitcodes;
mod generate;
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

    let derive_instrument_cmd = Command::new("instrument")
        .about("Derives the instrument used to produce the file (only Illumina is supported)")
        .arg(Arg::new("src").help("Source file.").index(1).required(true))
        .arg(
            Arg::new("first-n-reads")
                .short('n')
                .long("first-n-reads")
                .takes_value(true)
                .help(
                    "Only consider the first n reads in the file. If less \
                      than or equal to zero, the whole file will be read.",
                ),
        );

    let derive_cmd = Command::new("derive")
        .about("Forensic analysis tool useful for backwards computing information from next-generation sequencing data.")
        .subcommand_required(true)
        .subcommand(derive_instrument_cmd);

    let flagstat_cmd = Command::new("flagstat")
        .about("Generates stats from the marked flags in a SAM/BAM/CRAM file.")
        .arg(Arg::new("src").help("Source file.").index(1).required(true));

    let generate_cmd = Command::new("generate")
        .about("Generates a BAM file from a given reference genome.")
        .arg(
            Arg::new("reads-one-file")
                .long("--reads-one-file")
                .takes_value(true)
                .help("Destination for the first reads FASTQ file.")
                .required(true),
        )
        .arg(
            Arg::new("reads-two-file")
                .long("--reads-two-file")
                .takes_value(true)
                .help("Destination for the second reads FASTQ file.")
                .required(true),
        )
        .arg(
            Arg::new("reference")
                .long("--reference-fasta")
                .takes_value(true)
                .help("Reference FASTA to generate the data based off of.")
                .required(true),
        )
        .arg(Arg::new("error-rate")
                 .short('e')
                 .long("--error-rate")
                 .takes_value(true)
                 .default_value("0.0001")
                 .help("The error rate for the sequencer as a fraction between [0.0, 1.0] (per base).")
        )
        .arg(
            Arg::new("num-reads")
                .short('n')
                .long("--num-reads")
                .takes_value(true)
                .help("Specifies the exact number of read pairs to generate.")
                .conflicts_with("coverage"),
        )
        .arg(
            Arg::new("coverage")
                .short('c')
                .long("--coverage")
                .takes_value(true)
                .help("Dynamically calculate the number of reads needed for a particular mean coverage.")
                .conflicts_with("num-reads"),
        )
        .group(
            ArgGroup::new("reads-count")
                .arg("coverage")
                .arg("num-reads")
                .required(true),
        );

    let matches = Command::new("ngs")
        .version(version.as_str())
        .propagate_version(true)
        .subcommand_required(true)
        .subcommand(add_verbosity_args(derive_cmd))
        .subcommand(add_verbosity_args(flagstat_cmd))
        .subcommand(add_verbosity_args(generate_cmd))
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
            "flagstat" => return commands::flagstat(subcommand),
            "generate" => return commands::generate(subcommand),
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
