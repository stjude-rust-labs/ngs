use clap::{Arg, ArgGroup, Command};
use std::io;

mod commands;
mod derive;
mod exitcodes;
mod generate;
mod utils;

use git_testament::{git_testament, render_testament};

git_testament!(TESTAMENT);

fn main() -> io::Result<()> {
    let version = render_testament!(TESTAMENT);

    let derive_instrument_cmd = Command::new("instrument")
        .about("Derives the instrument used to produce the file (only Illumina is supported)")
        .arg(
            Arg::new("quiet")
                .short('q')
                .long("quiet")
                .help("Only errors are print to stderr."),
        )
        .arg(Arg::new("verbose").short('v').long("verbose").help(
            "All available information, including debug information, is \
                printed to stderr.",
        ))
        .arg(Arg::new("src").help("Source file.").index(1).required(true))
        .arg(
            Arg::new("first_n_reads")
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
        .subcommand(derive_cmd)
        .subcommand(flagstat_cmd)
        .subcommand(generate_cmd)
        .get_matches();

    let mut level = tracing::Level::INFO;
    if matches.value_of("quiet").is_some() {
        level = tracing::Level::ERROR;
    } else if matches.value_of("verbose").is_some() {
        level = tracing::Level::DEBUG;
    }

    let subscriber = tracing_subscriber::fmt::Subscriber::builder()
        .with_max_level(level)
        .with_writer(std::io::stderr)
        .finish();
    let _ = tracing::subscriber::set_global_default(subscriber);

    if let Some(derive) = matches.subcommand_matches("derive") {
        if let Some(m) = derive.subcommand_matches("instrument") {
            commands::derive::instrument::derive(m)
        } else {
            unreachable!();
        }
    } else if let Some(m) = matches.subcommand_matches("flagstat") {
        commands::flagstat(m)
    } else if let Some(m) = matches.subcommand_matches("generate") {
        commands::generate(m)
    } else {
        unreachable!();
    }
}
