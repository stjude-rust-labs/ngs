//! Functionality relating to the `ngs derive junction_annotation` subcommand itself.

use std::collections::HashMap;
use std::path::PathBuf;

use anyhow::Context;
use clap::Args;
use num_format::Locale;
use num_format::ToFormattedString;
use tracing::debug;
use tracing::info;

use crate::derive::junction_annotation::compute;
use crate::derive::junction_annotation::results::JunctionAnnotationResults;
use crate::utils::display::RecordCounter;
use crate::utils::formats;
use crate::utils::formats::bam::ParsedBAMFile;
use crate::utils::formats::utils::IndexCheck;

/// Clap arguments for the `ngs derive junction_annotation` subcommand.
#[derive(Args)]
pub struct JunctionAnnotationArgs {
    /// Source BAM.
    #[arg(value_name = "BAM")]
    src: PathBuf,

    /// Features GFF file.
    #[arg(short = 'f', long, required = true, value_name = "PATH")]
    features_gff: PathBuf,

    /// Name of the exon region feature for the gene model used.
    #[arg(short, long, value_name = "STRING", default_value = "exon")]
    pub exon_feature_name: String,

    /// Minimum intron length to consider.
    /// An intron is defined as an `N` CIGAR operation of any length.
    #[arg(long, value_name = "USIZE", default_value = "50")]
    pub min_intron_length: usize,

    /// Add +- this amount to intron positions when looking up exon positions.
    #[arg(long, value_name = "U8", default_value = "0")]
    pub fuzzy_junction_match_range: u8,

    /// Minimum number of reads supporting a junction to be considered.
    #[arg(long, value_name = "U8", default_value = "2")]
    pub min_read_support: u8,

    /// Minumum mapping quality for a record to be considered.
    #[arg(short, long, value_name = "U8", default_value = "30")]
    pub min_mapq: u8,
}

/// Main function for the `ngs derive junction_annotation` subcommand.
pub fn derive(args: JunctionAnnotationArgs) -> anyhow::Result<()> {
    info!("Starting derive junction_annotation subcommand.");

    let mut exon_starts: HashMap<&str, Vec<usize>> = HashMap::new();
    let mut exon_ends: HashMap<&str, Vec<usize>> = HashMap::new();

    // (1) Parse the GFF file and collect all exon features.
    debug!("Reading all records in GFF.");
    let mut gff = formats::gff::open(&args.features_gff)
        .with_context(|| format!("opening GFF file: {}", args.features_gff.display()))?;

    let mut exon_records = Vec::new();
    for result in gff.records() {
        let record = result.unwrap();
        if record.ty() != args.exon_feature_name {
            continue;
        }
        exon_records.push(record);
    }

    debug!("Tabulating GFF exon features.");
    for record in &exon_records {
        let seq_name = record.reference_sequence_name();
        let start = record.start().into();
        let end = record.end().into();

        exon_starts.entry(seq_name).or_default().push(start);
        exon_ends.entry(seq_name).or_default().push(end);
    }

    debug!("Finalizing GFF features lookup.");
    for starts in exon_starts.values_mut() {
        starts.sort_unstable();
        starts.dedup();
    }
    for ends in exon_ends.values_mut() {
        ends.sort_unstable();
        ends.dedup();
    }

    debug!("Done reading GFF.");

    let mut counter = RecordCounter::new();
    let mut results = JunctionAnnotationResults::default();
    let params = compute::JunctionAnnotationParameters {
        min_intron_length: args.min_intron_length,
        fuzzy_junction_match_range: args.fuzzy_junction_match_range,
        min_read_support: args.min_read_support,
        min_mapq: args.min_mapq,
    };

    let ParsedBAMFile {
        mut reader, header, ..
    } = formats::bam::open_and_parse(args.src, IndexCheck::None)?;

    // (2) Process each record in the BAM file.
    for result in reader.records(&header.parsed) {
        let record = result?;
        compute::process(
            &record,
            &exon_starts,
            &exon_ends,
            &header.parsed,
            &params,
            &mut results,
        )?;
        counter.inc();
    }

    info!(
        "Processed {} records.",
        counter.get().to_formatted_string(&Locale::en)
    );

    // (3) Summarize found junctions.
    compute::summarize(&mut results, &params)?;

    // (3) Print the output to stdout as JSON (more support for different output
    // types may be added in the future, but for now, only JSON).
    let output = serde_json::to_string_pretty(&results).unwrap();
    print!("{}", output);

    Ok(())
}
