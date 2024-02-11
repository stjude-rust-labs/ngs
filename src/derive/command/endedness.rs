//! Functionality relating to the `ngs derive endedness` subcommand itself.

use std::collections::HashMap;
use std::collections::HashSet;
use std::path::PathBuf;
use std::sync::Arc;

use anyhow::Context;
use clap::Args;
use num_format::Locale;
use num_format::ToFormattedString;
use tracing::info;
use tracing::trace;

use crate::derive::endedness::compute;
use crate::derive::endedness::compute::OrderingFlagsCounts;
use crate::utils::args::arg_in_range as deviance_in_range;
use crate::utils::args::NumberOfRecords;
use crate::utils::display::RecordCounter;
use crate::utils::formats::bam::ParsedBAMFile;
use crate::utils::formats::utils::IndexCheck;
use crate::utils::read_groups::{get_read_group, validate_read_group_info, ReadGroupPtr};

/// Clap arguments for the `ngs derive endedness` subcommand.
#[derive(Args)]
pub struct DeriveEndednessArgs {
    /// Source BAM.
    #[arg(value_name = "BAM")]
    src: PathBuf,

    /// Only examine the first n records in the file.
    #[arg(short, long, value_name = "USIZE")]
    num_records: Option<usize>,

    /// Distance from 0.5 split between number of f+l- reads and f-l+ reads
    /// allowed to be called 'Paired-End'. Default of `0.0` only appropriate
    /// if the whole file is being processed.
    #[arg(long, value_name = "F64", default_value = "0.0")]
    paired_deviance: f64,

    /// Calculate and output Reads-Per-Template. This will produce a more
    /// sophisticated estimate for endedness, but uses substantially more memory.
    #[arg(long, default_value = "false")]
    calc_rpt: bool,

    /// Round RPT to the nearest INT before comparing to expected values.
    /// Appropriate if using `-n` > 0.
    #[arg(long, default_value = "false")]
    round_rpt: bool,
}

/// Main function for the `ngs derive endedness` subcommand.
pub fn derive(args: DeriveEndednessArgs) -> anyhow::Result<()> {
    // (0) Parse arguments needed for subcommand.
    let paired_deviance = deviance_in_range(args.paired_deviance, 0.0..=0.5)
        .with_context(|| "Paired deviance is not within acceptable range")?;

    info!("Starting derive endedness subcommand.");

    let mut found_rgs = HashSet::new();

    let mut ordering_flags: HashMap<ReadGroupPtr, OrderingFlagsCounts> = HashMap::new();

    // only used if args.calc_rpt is true
    let mut read_names: HashMap<String, Vec<ReadGroupPtr>> = HashMap::new();

    let ParsedBAMFile {
        mut reader, header, ..
    } = crate::utils::formats::bam::open_and_parse(args.src, IndexCheck::None)?;

    // (1) Collect ordering flags (and QNAMEs) from reads within the
    // file. Support for sampling only a portion of the reads is provided.
    let num_records = NumberOfRecords::from(args.num_records);
    let mut counter = RecordCounter::default();

    for result in reader.records(&header.parsed) {
        let record = result?;

        // Only count primary alignments and unmapped reads.
        if (record.flags().is_secondary() || record.flags().is_supplementary())
            && !record.flags().is_unmapped()
        {
            continue;
        }

        let read_group = get_read_group(&record, Some(&mut found_rgs));

        if args.calc_rpt {
            match record.read_name() {
                Some(rn) => {
                    let rn = rn.to_string();
                    let rg_vec = read_names.get_mut(&rn);

                    match rg_vec {
                        Some(rg_vec) => {
                            rg_vec.push(Arc::clone(&read_group));
                        }
                        None => {
                            read_names.insert(rn, vec![(Arc::clone(&read_group))]);
                        }
                    }
                }
                None => {
                    trace!("Could not parse a QNAME from a read in the file.");
                    trace!("Skipping this read and proceeding.");
                    continue;
                }
            }
        }

        if !record.flags().is_segmented() {
            ordering_flags.entry(read_group).or_default().unsegmented += 1;
        } else if record.flags().is_first_segment() && !record.flags().is_last_segment() {
            ordering_flags.entry(read_group).or_default().first += 1;
        } else if !record.flags().is_first_segment() && record.flags().is_last_segment() {
            ordering_flags.entry(read_group).or_default().last += 1;
        } else if record.flags().is_first_segment() && record.flags().is_last_segment() {
            ordering_flags.entry(read_group).or_default().both += 1;
        } else if !record.flags().is_first_segment() && !record.flags().is_last_segment() {
            ordering_flags.entry(read_group).or_default().neither += 1;
        } else {
            unreachable!();
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

    // (1.5) Validate the read group information.
    let rgs_in_header_not_records = validate_read_group_info(&found_rgs, &header.parsed);
    for rg_id in rgs_in_header_not_records {
        ordering_flags.insert(Arc::new(rg_id), OrderingFlagsCounts::new());
    }

    // (2) Derive the endedness based on the ordering flags gathered.
    let result = compute::predict(ordering_flags, read_names, paired_deviance, args.round_rpt);

    // (3) Print the output to stdout as JSON (more support for different output
    // types may be added in the future, but for now, only JSON).
    let output = serde_json::to_string_pretty(&result).unwrap();
    println!("{}", output);

    anyhow::Ok(())
}
