//! Functionality relating to the `ngs derive strandedness` subcommand itself.

use std::collections::HashMap;
use std::fs::File;
use std::path::PathBuf;

use anyhow::bail;
use anyhow::Context;
use clap::Args;
use noodles::bam;
use noodles::gff;
use noodles::sam;
use rust_lapper::{Interval, Lapper};
use tracing::debug;
use tracing::info;

use crate::derive::strandedness::compute;
use crate::derive::strandedness::compute::ParsedBAMFile;
use crate::derive::strandedness::compute::StrandednessFilters;
use crate::utils::formats;

/// Clap arguments for the `ngs derive strandedness` subcommand.
#[derive(Args)]
pub struct DeriveStrandednessArgs {
    /// Source BAM.
    #[arg(value_name = "BAM")]
    src: PathBuf,

    /// Features GFF file.
    #[arg(short = 'f', long, required = true, value_name = "PATH")]
    features_gff: PathBuf,

    /// When inconclusive, the test will repeat until this many tries have been reached.
    /// Evidence of previous attempts is saved and reused,
    /// leading to a larger sample size with multiple attempts.
    #[arg(long, value_name = "USIZE", default_value = "3")]
    max_tries: usize,

    /// Filter any genes that don't have at least `m` reads.
    #[arg(short = 'm', long, value_name = "USIZE", default_value = "10")]
    min_reads_per_gene: usize,

    /// How many genes to sample.
    #[arg(short = 'n', long, value_name = "USIZE", default_value = "1000")]
    num_genes: usize,

    /// Minimum mapping quality for a record to be considered.
    /// Set to 0 to disable this filter and allow reads _without_
    /// a mapping quality to be considered.
    #[arg(short = 'q', long, value_name = "U8", default_value = "30")]
    min_mapq: u8,

    /// Consider all genes, not just protein coding genes.
    #[arg(long)]
    all_genes: bool,

    /// Name of the gene region feature for the gene model used.
    #[arg(long, value_name = "STRING", default_value = "gene")]
    gene_feature_name: String,

    /// Name of the exon region feature for the gene model used.
    #[arg(long, value_name = "STRING", default_value = "exon")]
    exon_feature_name: String,

    /// Do not count supplementary alignments.
    #[arg(long)]
    no_supplementary: bool,

    /// Do count secondary alignments.
    #[arg(long)]
    count_secondary: bool,

    /// Do count duplicates.
    #[arg(long)]
    count_duplicates: bool,

    /// Do count QC failed reads.
    #[arg(long)]
    count_qc_failed: bool,

    /// At most, search this many times for genes that satisfy our search criteria.
    /// Default is 10 * --num-genes.
    #[arg(long, value_name = "USIZE")]
    max_iterations_per_try: Option<usize>,
}

/// Main function for the `ngs derive strandedness` subcommand.
pub fn derive(args: DeriveStrandednessArgs) -> anyhow::Result<()> {
    info!("Starting derive strandedness subcommand.");

    // (1) Parse the GFF file and collect all gene features.
    debug!("Reading all records in GFF.");
    let mut gff = formats::gff::open(&args.features_gff)
        .with_context(|| format!("opening GFF file: {}", args.features_gff.display()))?;

    let mut gene_records = Vec::new();
    let mut exon_records = Vec::new();
    let mut gene_metrics = compute::GeneRecordMetrics::default();
    let mut exon_metrics = compute::ExonRecordMetrics::default();
    for result in gff.records() {
        let record = result.unwrap();
        if record.ty() == args.gene_feature_name {
            // If --all-genes is set, keep the record.
            // Otherwise, check the gene type or biotype and keep the record if it's protein coding.
            // If the record does not have a gene type or biotype, discard it.
            let mut keep_record = false;
            if !args.all_genes {
                let mut gene_type_value = None;
                for entry in record.attributes().as_ref() {
                    gene_type_value = match entry.key() {
                        "gene_type" => Some(entry.value()),    // Gencode
                        "gene_biotype" => Some(entry.value()), // ENSEMBL
                        "biotype" => Some(entry.value()),      // also ENSEMBL
                        _ => gene_type_value,
                    };
                }
                if let Some(gene_type_value) = gene_type_value {
                    if gene_type_value.to_lowercase().contains("protein") {
                        keep_record = true;
                        gene_metrics.protein_coding += 1;
                    }
                }
            }
            gene_metrics.total += 1;
            if keep_record {
                gene_records.push(record);
            }
        } else if record.ty() == args.exon_feature_name {
            exon_metrics.total += 1;
            exon_records.push(record);
        }
    }

    debug!("Tabulating GFF gene and exon features.");

    if gene_records.is_empty() {
        bail!("No gene records matched criteria. Check your GFF file and `--gene-feature-name` and `--all-genes` options.");
    }
    if exon_records.is_empty() {
        bail!("No exon records matched criteria. Check your GFF file and `--exon-feature-name` option.");
    }

    let mut exon_intervals: HashMap<&str, Vec<Interval<usize, gff::record::Strand>>> =
        HashMap::new();
    for record in &exon_records {
        let seq_name = record.reference_sequence_name();
        let start: usize = record.start().into();
        let stop: usize = record.end().into();
        let strand = record.strand();

        exon_intervals.entry(seq_name).or_default().push(Interval {
            start,
            stop,
            val: strand,
        });
    }

    let mut exons: HashMap<&str, Lapper<usize, gff::record::Strand>> = HashMap::new();
    for (seq_name, intervals) in exon_intervals {
        exons.insert(seq_name, Lapper::new(intervals));
    }

    debug!("Done reading GFF.");

    let mut reader = File::open(&args.src)
        .map(bam::Reader::new)
        .with_context(|| format!("opening BAM file: {}", args.src.display()))?;
    let header = reader.read_header()?.parse()?;
    let index = bam::bai::read(&args.src.with_extension("bam.bai")).with_context(|| {
        format!(
            "reading BAM index: {}",
            args.src.with_extension("bam.bai").display()
        )
    })?;

    let parsed_bam = ParsedBAMFile {
        reader,
        header,
        index,
    };

    let filters = StrandednessFilters {
        min_reads_per_gene: args.min_reads_per_gene,
        min_mapq: args.min_mapq,
        count_qc_failed: args.count_qc_failed,
        no_supplementary: args.no_supplementary,
        count_secondary: args.count_secondary,
        count_duplicates: args.count_duplicates,
    };

    let max_iterations_per_try = args.max_iterations_per_try.unwrap_or(args.num_genes * 10);
    let max_iterations_per_try = match max_iterations_per_try > gene_records.len() {
        true => gene_records.len(),
        false => max_iterations_per_try,
    };

    for try_num in 1..=args.max_tries {
        info!("Starting try {} of {}", try_num, args.max_tries);

        compute::predict(
            &parsed_bam,
            &mut gene_records,
            &exons,
            max_iterations_per_try,
            args.num_genes,
            &filters,
            &mut gene_metrics,
            &mut exon_metrics,
        )?;
    }

    anyhow::Ok(())
}
