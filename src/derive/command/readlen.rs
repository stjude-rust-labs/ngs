//! Functionality relating to the `ngs derive readlen` subcommand itself.

use std::collections::HashMap;
use std::path::PathBuf;

use clap::Args;
use num_format::Locale;
use num_format::ToFormattedString;
use tracing::info;

use crate::derive::readlen::compute;
use crate::utils::args::NumberOfRecords;
use crate::utils::display::RecordCounter;
use crate::utils::formats::bam::ParsedBAMFile;
use crate::utils::formats::utils::IndexCheck;

/// Utility method to parse the Majority Vote Cutoff passed in on the command line and
/// ensure the cutoff is within the range [0.0, 1.0].
pub fn cutoff_in_range(cutoff_raw: &str) -> Result<f64, String> {
    let cutoff: f64 = cutoff_raw
        .parse()
        .map_err(|_| format!("{} isn't a float", cutoff_raw))?;

    match (0.0..=1.0).contains(&cutoff) {
        true => Ok(cutoff),
        false => Err(String::from(
            "Majority Vote Cutoff must be between 0.0 and 1.0",
        )),
    }
}

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
    #[arg(value_parser = cutoff_in_range)]
    majority_vote_cutoff: Option<f64>,
}

/// Main function for the `ngs derive readlen` subcommand.
pub fn derive(args: DeriveReadlenArgs) -> anyhow::Result<()> {
    let mut read_lengths = HashMap::new();

    info!("Starting derive readlen subcommand.");

    let ParsedBAMFile {
        mut reader, header, ..
    } = crate::utils::formats::bam::open_and_parse(args.src, IndexCheck::None)?;

    // (1) Collect read lengths from reads within the
    // file. Support for sampling only a portion of the reads is provided.
    let num_records = NumberOfRecords::from(args.num_records);
    let mut counter = RecordCounter::new();

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
    let result = compute::predict(
        read_lengths,
        counter.get(),
        args.majority_vote_cutoff.unwrap(),
    )
    .unwrap();

    // (3) Print the output to stdout as JSON (more support for different output
    // types may be added in the future, but for now, only JSON).
    let output = serde_json::to_string_pretty(&result).unwrap();
    print!("{}", output);

    anyhow::Ok(())
}
