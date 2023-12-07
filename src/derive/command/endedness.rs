//! Functionality relating to the `ngs derive endedness` subcommand itself.

use std::collections::HashMap;
use std::collections::HashSet;
use std::path::PathBuf;
use std::rc::Rc;

use clap::Args;
use noodles::sam::record::data::field::Tag;
use radix_trie::Trie;
use tracing::info;
use tracing::trace;

use crate::derive::endedness::compute;
use crate::derive::endedness::compute::{OrderingFlagsCounts, OVERALL, UNKNOWN_READ_GROUP};
use crate::utils::formats::bam::ParsedBAMFile;
use crate::utils::formats::utils::IndexCheck;

/// Utility method to parse the Paired Deviance passed in on the command line and
/// ensure the value is within the range [0.0, 0.5].
pub fn deviance_in_range(deviance_raw: &str) -> Result<f64, String> {
    let deviance: f64 = deviance_raw
        .parse()
        .map_err(|_| format!("{} isn't a float", deviance_raw))?;

    match (0.0..=0.5).contains(&deviance) {
        true => Ok(deviance),
        false => Err(String::from("Paired Deviance must be between 0.0 and 0.5")),
    }
}

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
    #[arg(value_parser = deviance_in_range)]
    paired_deviance: Option<f64>,

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
    info!("Starting derive endedness subcommand.");

    let mut ordering_flags: HashMap<Rc<String>, OrderingFlagsCounts> = HashMap::new();
    ordering_flags.insert(Rc::new(OVERALL.to_string()), OrderingFlagsCounts::new());
    ordering_flags.insert(
        Rc::new(UNKNOWN_READ_GROUP.to_string()),
        OrderingFlagsCounts::new(),
    );

    // only used if args.calc_rpt is true
    let mut found_rgs = HashSet::new();
    let mut read_names = Trie::<String, Vec<Rc<String>>>::new();

    let ParsedBAMFile {
        mut reader, header, ..
    } = crate::utils::formats::bam::open_and_parse(args.src, IndexCheck::Full)?;

    // (1) Collect read lengths from reads within the
    // file. Support for sampling only a portion of the reads is provided.
    let mut samples = 0;
    let mut sample_max = 0;

    if let Some(s) = args.num_records {
        sample_max = s;
    }

    for result in reader.records(&header.parsed) {
        let record = result?;

        // Only count primary alignments and unmapped reads.
        if (record.flags().is_secondary() || record.flags().is_supplementary())
            && !record.flags().is_unmapped()
        {
            continue;
        }

        let read_group = match record.data().get(Tag::ReadGroup) {
            Some(rg) => Rc::new(rg.as_str().unwrap().to_owned()),
            None => Rc::new(UNKNOWN_READ_GROUP.to_string()),
        };

        if args.calc_rpt {
            found_rgs.insert(Rc::clone(&read_group));

            match record.read_name() {
                Some(rn) => {
                    let rg_vec = read_names.get_mut(&rn.to_string());

                    match rg_vec {
                        Some(rg_vec) => {
                            rg_vec.push(Rc::clone(&read_group));
                        }
                        None => {
                            let rg_vec = vec![(Rc::clone(&read_group))];
                            read_names.insert(rn.to_string(), rg_vec);
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

        if record.flags().is_first_segment() && !record.flags().is_last_segment() {
            ordering_flags
                .entry(Rc::new(OVERALL.to_string()))
                .and_modify(|e| {
                    e.first += 1;
                });

            ordering_flags
                .entry(read_group)
                .and_modify(|e| {
                    e.first += 1;
                })
                .or_insert(OrderingFlagsCounts {
                    first: 1,
                    last: 0,
                    both: 0,
                    neither: 0,
                });
        } else if !record.flags().is_first_segment() && record.flags().is_last_segment() {
            ordering_flags
                .entry(Rc::new(OVERALL.to_string()))
                .and_modify(|e| {
                    e.last += 1;
                });

            ordering_flags
                .entry(read_group)
                .and_modify(|e| {
                    e.last += 1;
                })
                .or_insert(OrderingFlagsCounts {
                    first: 0,
                    last: 1,
                    both: 0,
                    neither: 0,
                });
        } else if record.flags().is_first_segment() && record.flags().is_last_segment() {
            ordering_flags
                .entry(Rc::new(OVERALL.to_string()))
                .and_modify(|e| {
                    e.both += 1;
                });

            ordering_flags
                .entry(read_group)
                .and_modify(|e| {
                    e.both += 1;
                })
                .or_insert(OrderingFlagsCounts {
                    first: 0,
                    last: 0,
                    both: 1,
                    neither: 0,
                });
        } else if !record.flags().is_first_segment() && !record.flags().is_last_segment() {
            ordering_flags
                .entry(Rc::new(OVERALL.to_string()))
                .and_modify(|e| {
                    e.neither += 1;
                });

            ordering_flags
                .entry(read_group)
                .and_modify(|e| {
                    e.neither += 1;
                })
                .or_insert(OrderingFlagsCounts {
                    first: 0,
                    last: 0,
                    both: 0,
                    neither: 1,
                });
        } else {
            unreachable!();
        }

        if sample_max > 0 {
            samples += 1;
            if samples > sample_max {
                break;
            }
        }
    }

    // (2) Derive the consensus endedness based on the ordering flags gathered.
    let result = compute::predict(
        ordering_flags,
        read_names,
        args.paired_deviance.unwrap(),
        args.round_rpt,
    )
    .unwrap();

    // (3) Print the output to stdout as JSON (more support for different output
    // types may be added in the future, but for now, only JSON).
    let output = serde_json::to_string_pretty(&result).unwrap();
    print!("{}", output);

    anyhow::Ok(())
}
