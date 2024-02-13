//! Functionality relating to the `ngs derive instrument` subcommand itself.

use anyhow::bail;
use clap::Args;
use num_format::{Locale, ToFormattedString};
use std::collections::HashSet;
use std::path::PathBuf;
use tracing::info;

use crate::derive::instrument::compute;
use crate::derive::instrument::reads::IlluminaReadName;
use crate::utils::args::NumberOfRecords;
use crate::utils::display::RecordCounter;
use crate::utils::formats::bam::ParsedBAMFile;
use crate::utils::formats::utils::IndexCheck;

/// Clap arguments for the `ngs derive instrument` subcommand.
#[derive(Args)]
pub struct DeriveInstrumentArgs {
    // Source BAM.
    #[arg(value_name = "BAM")]
    src: PathBuf,

    /// Only examine the first n records in the file.
    #[arg(short, long, value_name = "USIZE")]
    num_records: Option<usize>,
}

/// Main function for the `ngs derive instrument` subcommand.
pub fn derive(args: DeriveInstrumentArgs) -> anyhow::Result<()> {
    let src = args.src;
    let mut instrument_names = HashSet::new();
    let mut flowcell_names = HashSet::new();
    let mut metrics = compute::RecordMetrics::default();

    info!("Starting derive instrument subcommand.");

    let ParsedBAMFile {
        mut reader, header, ..
    } = crate::utils::formats::bam::open_and_parse(src, IndexCheck::None)?;

    // (1) Collect instrument names and flowcell names from reads within the
    // file. Support for sampling only a portion of the reads is provided.
    let num_records = NumberOfRecords::from(args.num_records);
    let mut counter = RecordCounter::default();

    for result in reader.records(&header.parsed) {
        let record = result?;

        if let Some(read_name) = record.read_name() {
            let name: &str = read_name.as_ref();

            match name.parse::<IlluminaReadName>() {
                Ok(read) => {
                    instrument_names.insert(read.instrument_name);
                    metrics.found_instrument_name += 1;
                    if let Some(fc) = read.flowcell {
                        flowcell_names.insert(fc);
                        metrics.found_flowcell_name += 1;
                    }
                }
                Err(_) => {
                    bail!(
                        "Could not parse Illumina-formatted query names for read: {}",
                        name
                    );
                }
            }
        } else {
            metrics.bad_read_name += 1;
        }

        counter.inc();
        if counter.time_to_break(&num_records) {
            break;
        }
    }

    info!(
        "Processed {} records.",
        counter.get().to_formatted_string(&Locale::en)
    );
    metrics.total_records = counter.get();
    metrics.unique_instrument_names = instrument_names.clone();
    metrics.unique_flowcell_names = flowcell_names.clone();

    // (2) Derive the predict instrument results based on these detected
    // instrument names and flowcell names.
    let mut result = compute::predict(instrument_names, flowcell_names);
    result.records = metrics;

    // (3) Print the output to stdout as JSON (more support for different output
    // types may be added in the future, but for now, only JSON).
    let output = serde_json::to_string_pretty(&result).unwrap();
    println!("{}", output);

    Ok(())
}
