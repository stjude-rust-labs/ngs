use std::path::PathBuf;

use anyhow::Context;
use clap::{value_parser, Arg, ArgGroup, ArgMatches, Command};
use indicatif::{ProgressBar, ProgressStyle};
use rand::prelude::*;
use tracing::info;

use crate::{
    generate::providers::{reference_provider::ReferenceGenomeSequenceProvider, SequenceProvider},
    utils::formats,
};

pub fn get_command() -> Command {
    Command::new("generate")
        .about("Generates a BAM file from a given reference genome.")
        .arg(
            Arg::new("reads-one-file")
                .long("--reads-one-file")
                .help("Destination for the first reads FASTQ file.")
                .value_parser(value_parser!(PathBuf))
                .required(true),
        )
        .arg(
            Arg::new("reads-two-file")
                .long("--reads-two-file")
                .help("Destination for the second reads FASTQ file.")
                .value_parser(value_parser!(PathBuf))
                .required(true),
        )
        .arg(
            Arg::new("reference-providers")
                .long("--reference-provider")
                .help("One or more reference FASTAs to generate the data based off of.")
                .required(true),
        )
        .arg(Arg::new("error-rate")
                 .short('e')
                 .long("--error-rate")
                 .default_value("0.0001")
                 .value_parser(value_parser!(f64))
                 .help("The error rate for the sequencer as a fraction between [0.0, 1.0] (per base).")
        )
        .arg(
            Arg::new("num-records")
                .short('n')
                .long("--num-records")
                .help("Specifies the number of records to generate.")
                .value_parser(value_parser!(usize))
                .conflicts_with("coverage"),
        )
        .arg(
            Arg::new("coverage")
                .short('c')
                .long("--coverage")
                .help("Dynamically calculate the number of reads needed for a particular mean coverage.")
                .value_parser(value_parser!(usize))
                .conflicts_with("num-reads"),
        )
        .group(
            ArgGroup::new("reads-count")
                .arg("coverage")
                .arg("num-reads")
                .required(true),
        )
}

pub fn generate(matches: &ArgMatches) -> anyhow::Result<()> {
    // (0) Parse arguments needed for subcommand.
    let reference_providers: Vec<_> = matches
        .get_many::<String>("reference-providers")
        .expect("missing reference providers")
        .map(|rp| rp.parse::<ReferenceGenomeSequenceProvider>().unwrap())
        .collect();

    let reads_one_file = matches
        .get_one::<PathBuf>("reads-one-file")
        .expect("missing reads one file");

    let reads_two_file = matches
        .get_one::<PathBuf>("reads-two-file")
        .expect("missing reads two file");

    info!("Starting generate command...");
    let mut writer_read_one = formats::fastq::writer(&reads_one_file).with_context(|| {
        format!(
            "Could not open reads one file: {}.",
            reads_one_file.display()
        )
    })?;
    let mut writer_read_two = formats::fastq::writer(&reads_two_file).with_context(|| {
        format!(
            "Could not open reads two file: {}.",
            reads_two_file.display()
        )
    })?;

    let mut total_reads: usize = 0;

    if let Some(num_reads) = matches.get_one::<usize>("num-reads") {
        total_reads = *num_reads
    } else if let Some(coverage) = matches.get_one::<usize>("coverage") {
        total_reads = reference_providers[0].reads_needed_for_coverage(*coverage);
    }

    // (2) Set up the output writers.
    info!("Generating {} reads...", total_reads);

    // (3) Set up the progress bar.
    let pb = ProgressBar::new(total_reads as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{prefix:.cyan.bold} {spinner:.green} [{elapsed_precise}] [{bar}] {pos}/{len} ({per_sec}, {eta})")
            .progress_chars("=> "),
    );
    pb.set_prefix("Generating");

    // (4) Generate the reads and write them to their respective files.
    let mut rng = ThreadRng::default();

    let mut i = 0;
    while i < total_reads {
        let selected_genome = reference_providers
            .choose_weighted(&mut rng, |x| x.weight())
            .unwrap();
        let read_pair = selected_genome.generate_read_pair(
            format!("ngs:{}", selected_genome.name()),
            (i + 1).to_string(),
        );
        writer_read_one
            .write_record(read_pair.get_forward_read())
            .with_context(|| "Could not write record to read one file.")?;
        writer_read_two
            .write_record(read_pair.get_reverse_read())
            .with_context(|| "Could not write record to read two file.")?;

        if i > 0 && i % 5_000 == 0 {
            pb.inc(5000);
        }

        i += 1;
    }

    pb.set_style(
        ProgressStyle::default_bar()
            .template("{prefix:.green.bold} {msg:.white.bold} [{elapsed_precise}] [{bar}] {pos}/{len} ({per_sec}, {eta})")
            .progress_chars("=> "),
    );
    pb.set_prefix("✓");
    pb.finish_with_message("Finished");
    Ok(())
}
