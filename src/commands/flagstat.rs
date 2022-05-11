//! The `flagstat` subcommand intends to be a drop-in replacement for the `flagstat`
//! subcommand provided by [`samtools`]. You can learn more about the output for this
//! command by reviewing [the relevant documentation][samtools-flagstat-docs] for the
//! samtools project.
//!
//! This vast majority of this code was pulled directly from the excellent example in
//! the [`noodles-bam`] package at
//! [this commit](https://github.com/zaeleus/noodles/blob/eb129f0/noodles-bam/examples/bam_flagstat.rs).
//!
//! ## Example Usage
//!
//! ```bash
//! ngs flagstat FILE.bam
//! ```
//!
//! [`samtools`]: https://github.com/samtools/samtools
//! [`noodles-bam`]: https://github.com/zaeleus/noodles/tree/master/noodles-bam
//! [samtools-flagstat-docs]: https://www.htslib.org/doc/samtools-flagstat.html
//!
use clap::ArgMatches;
use futures::TryStreamExt;
use sam::reader::record::Fields;
use tokio::{fs::File, io};
use tracing::info;

use crate::utils;

use noodles_bam as bam;
use noodles_sam::{self as sam, AlignmentRecord};

#[derive(Debug, Default)]
struct Counts {
    read: u64,
    primary: u64,
    secondary: u64,
    supplementary: u64,
    duplicate: u64,
    primary_duplicate: u64,
    mapped: u64,
    primary_mapped: u64,
    paired: u64,
    read_1: u64,
    read_2: u64,
    proper_pair: u64,
    mate_mapped: u64,
    singleton: u64,
    mate_reference_sequence_id_mismatch: u64,
    mate_reference_sequence_id_mismatch_hq: u64,
}

fn count(counts: &mut Counts, record: &bam::Record) {
    let flags = record.flags();

    counts.read += 1;

    if !flags.is_unmapped() {
        counts.mapped += 1;
    }

    if flags.is_duplicate() {
        counts.duplicate += 1;
    }

    if flags.is_secondary() {
        counts.secondary += 1;
    } else if flags.is_supplementary() {
        counts.supplementary += 1;
    } else {
        counts.primary += 1;

        if !flags.is_unmapped() {
            counts.primary_mapped += 1;
        }

        if flags.is_duplicate() {
            counts.primary_duplicate += 1;
        }

        if flags.is_segmented() {
            counts.paired += 1;

            if flags.is_first_segment() {
                counts.read_1 += 1;
            }

            if flags.is_last_segment() {
                counts.read_2 += 1;
            }

            if !flags.is_unmapped() {
                if flags.is_properly_aligned() {
                    counts.proper_pair += 1;
                }

                if flags.is_mate_unmapped() {
                    counts.singleton += 1;
                } else {
                    counts.mate_mapped += 1;

                    if record.mate_reference_sequence_id() != record.reference_sequence_id() {
                        counts.mate_reference_sequence_id_mismatch += 1;

                        let mapq = record
                            .mapping_quality()
                            .map(u8::from)
                            .unwrap_or(sam::record::mapping_quality::MISSING);

                        if mapq >= 5 {
                            counts.mate_reference_sequence_id_mismatch_hq += 1;
                        }
                    }
                }
            }
        }
    }
}

fn print_stats(qc_pass_counts: &Counts, qc_fail_counts: &Counts) {
    println!(
        "{} + {} in total (QC-passed reads + QC-failed reads)",
        qc_pass_counts.read, qc_fail_counts.read
    );
    println!(
        "{} + {} primary",
        qc_pass_counts.primary, qc_fail_counts.primary
    );
    println!(
        "{} + {} secondary",
        qc_pass_counts.secondary, qc_fail_counts.secondary
    );
    println!(
        "{} + {} supplementary",
        qc_pass_counts.supplementary, qc_fail_counts.supplementary
    );
    println!(
        "{} + {} duplicates",
        qc_pass_counts.duplicate, qc_fail_counts.duplicate
    );
    println!(
        "{} + {} primary duplicates",
        qc_pass_counts.primary_duplicate, qc_fail_counts.primary_duplicate
    );
    println!(
        "{} + {} mapped ({} : {})",
        qc_pass_counts.mapped,
        qc_fail_counts.mapped,
        utils::display::PercentageFormat(qc_pass_counts.mapped, qc_pass_counts.read),
        utils::display::PercentageFormat(qc_fail_counts.mapped, qc_fail_counts.read)
    );
    println!(
        "{} + {} primary mapped ({} : {})",
        qc_pass_counts.primary_mapped,
        qc_fail_counts.primary_mapped,
        utils::display::PercentageFormat(qc_pass_counts.primary_mapped, qc_pass_counts.primary),
        utils::display::PercentageFormat(qc_fail_counts.primary_mapped, qc_fail_counts.primary)
    );
    println!(
        "{} + {} paired in sequencing",
        qc_pass_counts.paired, qc_fail_counts.paired
    );
    println!(
        "{} + {} read1",
        qc_pass_counts.read_1, qc_fail_counts.read_1
    );
    println!(
        "{} + {} read2",
        qc_pass_counts.read_2, qc_fail_counts.read_2
    );
    println!(
        "{} + {} properly paired ({} : {})",
        qc_pass_counts.proper_pair,
        qc_fail_counts.proper_pair,
        utils::display::PercentageFormat(qc_pass_counts.proper_pair, qc_pass_counts.paired),
        utils::display::PercentageFormat(qc_fail_counts.proper_pair, qc_fail_counts.paired)
    );
    println!(
        "{} + {} with itself and mate mapped",
        qc_pass_counts.mate_mapped, qc_fail_counts.mate_mapped
    );
    println!(
        "{} + {} singletons ({} : {})",
        qc_pass_counts.singleton,
        qc_fail_counts.singleton,
        utils::display::PercentageFormat(qc_pass_counts.singleton, qc_pass_counts.paired),
        utils::display::PercentageFormat(qc_fail_counts.singleton, qc_fail_counts.paired)
    );
    println!(
        "{} + {} with mate mapped to a different chr",
        qc_pass_counts.mate_reference_sequence_id_mismatch,
        qc_fail_counts.mate_reference_sequence_id_mismatch
    );
    println!(
        "{} + {} with mate mapped to a different chr (mapQ>=5)",
        qc_pass_counts.mate_reference_sequence_id_mismatch_hq,
        qc_fail_counts.mate_reference_sequence_id_mismatch_hq
    );
}

async fn app(
    src: &str,
    qc_pass_counts: &mut Counts,
    qc_fail_counts: &mut Counts,
) -> io::Result<()> {
    let mut reader = File::open(src).await.map(bam::AsyncReader::new)?;
    reader.read_header().await?;
    reader.read_reference_sequences().await?;

    let fields = Fields::FLAGS;
    let mut records = reader.records_with_fields(fields);
    while let Some(record) = records.try_next().await? {
        if record.flags().is_qc_fail() {
            count(qc_fail_counts, &record);
        } else {
            count(qc_pass_counts, &record);
        }
    }

    Ok(())
}

pub fn flagstat(matches: &ArgMatches) -> io::Result<()> {
    let src = matches.value_of("src").unwrap();
    info!("Starting flagstat...");

    let rt = tokio::runtime::Builder::new_current_thread()
        .enable_all()
        .build()
        .unwrap();

    let mut qc_pass_counts = Counts::default();
    let mut qc_fail_counts = Counts::default();

    let app = app(src, &mut qc_pass_counts, &mut qc_fail_counts);
    rt.block_on(app).unwrap();

    print_stats(&qc_pass_counts, &qc_fail_counts);
    Ok(())
}
