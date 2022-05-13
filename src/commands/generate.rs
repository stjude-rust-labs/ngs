use std::{fs::File, io, path::PathBuf};

use crate::generate::provider::ReferenceGenomeSequenceProvider;
use clap::ArgMatches;
use flate2::{write::GzEncoder, Compression};
use generate::provider::SequenceProvider;
use indicatif::{ProgressBar, ProgressStyle};
use noodles_fastq as fastq;
use tracing::info;

use crate::generate;

pub fn generate(matches: &ArgMatches) -> io::Result<()> {
    // (0) Parse arguments needed for subcommand.
    let reference = matches
        .value_of("reference")
        .map(PathBuf::from)
        .expect("missing reference");
    info!("Starting generate command...");

    // (1) Read in all sequences in reference FASTA and then precache all of the
    // results before showing the progress bar.
    info!("Loading reference genome...");
    let mut human = ReferenceGenomeSequenceProvider::new(reference, 150, -150..150)?;
    human.precache();

    // (2) Set up the output writers.
    info!("Generating reads...");
    let mut writer_read_one = File::create("/Users/cmcleod/Desktop/Test_R1.fq.gz")
        .map(|x| GzEncoder::new(x, Compression::default()))
        .map(fastq::Writer::new)?;

    let mut writer_read_two = File::create("/Users/cmcleod/Desktop/Test_R2.fq.gz")
        .map(|x| GzEncoder::new(x, Compression::default()))
        .map(fastq::Writer::new)?;

    // (3) Set up the progress bar.
    let total_reads = 1_000_000;
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
