use std::{path::PathBuf};

use crate::{generate::provider::ReferenceGenomeSequenceProvider, formats};
use anyhow::Context;
use clap::{Arg, ArgGroup, ArgMatches, Command};
use generate::provider::SequenceProvider;
use indicatif::{ProgressBar, ProgressStyle};
use tracing::{info};

use crate::generate;

pub fn get_command<'a>() -> Command<'a> {
    Command::new("generate")
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
        )
}

pub fn generate(matches: &ArgMatches) -> anyhow::Result<()> {
    // (0) Parse arguments needed for subcommand.
    let reference = matches
        .value_of("reference")
        .map(PathBuf::from)
        .expect("missing reference");

    let reads_one_file = matches
        .value_of("reads-one-file")
        .map(PathBuf::from)
        .expect("missing reads one file");

    let reads_two_file = matches
        .value_of("reads-two-file")
        .map(PathBuf::from)
        .expect("missing reads two file");

    let error_rate = matches
        .value_of("error-rate")
        .map(|x| x.parse::<f64>())
        .expect("missing error rate")
        .expect("Could not parse error rate as a float.");
    let error_freq = (1.0 / error_rate) as usize;

    info!("Starting generate command...");
    let mut writer_read_one = formats::fastq::writer(&reads_one_file).with_context(|| format!("Could not open reads one file: {}.", reads_one_file.display()))?;
    let mut writer_read_two = formats::fastq::writer(&reads_two_file).with_context(|| format!("Could not open reads two file: {}.", reads_two_file.display()))?;

    // (1) Read in all sequences in reference FASTA and then precache all of the
    // results before showing the progress bar.
    info!("Loading reference genome...");
    let mut human = ReferenceGenomeSequenceProvider::new(reference, 150, -150..150, error_freq)?;

    let mut total_reads: usize = 0;
    if let Some(num_reads) = matches.value_of("num-reads") {
        total_reads = num_reads
            .parse()
            .expect("Could not parse number of reads from command line args.");
    } else if let Some(coverage) = matches.value_of("coverage") {
        let cov: usize = coverage
            .parse()
            .expect("Could not parse coverage as size from command line args.");
        total_reads = human.reads_needed_for_coverage(cov);
    }

    // (2) Set up the output writers.
    info!("Generating {} reads...", total_reads);

    // (3) Set up the progress bar.
    // let total_reads = 10;
    let pb = ProgressBar::new(total_reads as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{prefix:.cyan.bold} {spinner:.green} [{elapsed_precise}] [{bar}] {pos}/{len} ({per_sec}, {eta})")
            .progress_chars("=> "),
    );
    pb.set_prefix("Generating");

    // (4) Generate the reads and write them to their respective files.
    let mut i = 0;
    while i <= total_reads {
        let read_pair = human.generate_read_pair(format!("ngs:{}", i + 1));
        writer_read_one
            .write_record(read_pair.get_forward_read())
            .expect("Could not write record to read one file.");
        writer_read_two
            .write_record(read_pair.get_reverse_read())
            .expect("Could not write record to read two file.");

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
