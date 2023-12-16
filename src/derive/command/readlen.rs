//! Functionality relating to the `ngs derive readlen` subcommand itself.

use std::collections::HashMap;
use std::path::PathBuf;

use clap::Args;
use tracing::info;

use crate::derive::readlen::compute;
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
    let mut samples = 0;
    let mut sample_max = 0;

    if let Some(s) = args.num_records {
        sample_max = s;
    }

    for result in reader.records(&header.parsed) {
        let record = result?;
        let len = record.sequence().len();

        read_lengths.entry(len).and_modify(|e| *e += 1).or_insert(1);

        samples += 1;
        if sample_max > 0 && samples > sample_max {
            break;
        }
    }

    // (2) Derive the consensus read length based on the read lengths gathered.
    let result =
        compute::predict(read_lengths, samples, args.majority_vote_cutoff.unwrap()).unwrap();

    // (3) Print the output to stdout as JSON (more support for different output
    // types may be added in the future, but for now, only JSON).
    let output = serde_json::to_string_pretty(&result).unwrap();
    print!("{}", output);

    Ok(())
}
