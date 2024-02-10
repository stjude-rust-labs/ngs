//! Functionality relating to the `ngs derive readlen` subcommand itself.

use std::collections::HashMap;
use std::path::PathBuf;

use anyhow::Context;
use clap::Args;
use num_format::Locale;
use num_format::ToFormattedString;
use tracing::info;

use crate::derive::readlen::compute;
use crate::utils::args::arg_in_range as cutoff_in_range;
use crate::utils::args::NumberOfRecords;
use crate::utils::display::RecordCounter;
use crate::utils::formats::bam::ParsedBAMFile;
use crate::utils::formats::utils::IndexCheck;

/// Clap arguments for the `ngs derive readlen` subcommand.
#[derive(Args)]
pub struct DeriveReadlenArgs {
    // Source BAM.
    #[arg(value_name = "BAM")]
    src: PathBuf,

    /// Only examine the first n records in the file.
    #[arg(short, long, value_name = "USIZE")]
    num_records: Option<usize>,

    /// Majority vote cutoff value as a fraction between [0.0, 1.0].
    #[arg(short, long, value_name = "F64", default_value = "0.7")]
    majority_vote_cutoff: f64,
}

/// Main function for the `ngs derive readlen` subcommand.
pub fn derive(args: DeriveReadlenArgs) -> anyhow::Result<()> {
    // (0) Parse arguments needed for subcommand.
    let majority_vote_cutoff = cutoff_in_range(args.majority_vote_cutoff, 0.0..=1.0)
        .with_context(|| "Majority vote cutoff is not within acceptable range")?;

    let mut read_lengths = HashMap::new();

    info!("Starting derive readlen subcommand.");

    let ParsedBAMFile {
        mut reader, header, ..
    } = crate::utils::formats::bam::open_and_parse(args.src, IndexCheck::None)?;

    // (1) Collect read lengths from reads within the
    // file. Support for sampling only a portion of the reads is provided.
    let num_records = NumberOfRecords::from(args.num_records);
    let mut counter = RecordCounter::default();

    for result in reader.records(&header.parsed) {
        let record = result?;
        let len = record.sequence().len();

        read_lengths.entry(len).and_modify(|e| *e += 1).or_insert(1);

        counter.inc();
        if counter.time_to_break(&num_records) {
            break;
        }
    }

    info!(
        "Processed {} records.",
        counter.get().to_formatted_string(&Locale::en)
    );

    // (2) Derive the consensus read length based on the read lengths gathered.
    let result = compute::predict(read_lengths, counter.get(), majority_vote_cutoff).unwrap();

    // (3) Print the output to stdout as JSON (more support for different output
    // types may be added in the future, but for now, only JSON).
    let output = serde_json::to_string_pretty(&result).unwrap();
    print!("{}", output);

    anyhow::Ok(())
}
