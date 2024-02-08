//! Module holding the logic for computing the strandedness.

use anyhow::bail;
use noodles::bam;
use noodles::core::{Position, Region};
use noodles::gff;
use noodles::sam;
use noodles::sam::record::data::field::Tag;
use rand::Rng;
use rust_lapper::Lapper;
use serde::Serialize;
use std::collections::HashMap;
use std::collections::HashSet;
use std::sync::Arc;

use crate::utils::read_groups::{validate_read_group_info, OVERALL, UNKNOWN_READ_GROUP};

const STRANDED_THRESHOLD: f64 = 0.80;
const UNSTRANDED_THRESHOLD: f64 = 0.40;

/// General gene metrics that are tallied as a part of the
/// strandedness subcommand.
#[derive(Clone, Default, Serialize, Debug)]
pub struct GeneRecordMetrics {
    /// The total number of genes found in the GFF.
    pub total: usize,

    /// The number of genes that were found to be protein coding.
    /// If --all-genes is set this will not be tallied.
    pub protein_coding: usize,

    /// The number of genes tested.
    pub tested: usize,

    /// The number of genes which were discarded due to having
    /// an unknown/invalid strand OR with exons on both strands.
    pub bad_strands: usize,

    /// The number of genes which were discarded due to not having
    /// enough reads.
    pub not_enough_reads: usize,
}

/// General exon metrics that are tallied as a part of the
/// strandedness subcommand.
#[derive(Clone, Default, Serialize, Debug)]
pub struct ExonRecordMetrics {
    /// The total number of exons found in the GFF.
    pub total: usize,

    /// The number of exons discarded due to having an unknown/invalid strand.
    pub bad_strand: usize,
}

/// General read record metrics that are tallied as a part of the
/// strandedness subcommand.
#[derive(Clone, Default, Serialize, Debug)]
pub struct ReadRecordMetrics {
    /// The number of records that have been filtered because of their flags.
    /// (i.e. they were qc_fail, duplicates, secondary, or supplementary)
    /// These conditions can be toggled on/off with CL flags
    pub ignored_flags: usize,

    /// The number of records that have been filtered because
    /// they failed the MAPQ filter.
    pub low_mapq: usize,

    /// The number of records whose MAPQ couldn't be parsed and were thus discarded.
    pub missing_mapq: usize,

    /// The number of records determined to be Paired-End.
    pub paired_end_reads: usize,

    /// The number of records determined to be Single-End.
    pub single_end_reads: usize,
}

/// Struct for tracking count results.
#[derive(Clone, Default)]
struct Counts {
    /// The number of reads that are evidence of Forward Strandedness.
    forward: usize,

    /// The number of reads that are evidence of Reverse Strandedness.
    reverse: usize,
}

/// Struct holding the per read group results for an `ngs derive strandedness`
/// subcommand call.
#[derive(Debug, Serialize)]
pub struct ReadGroupDerivedStrandednessResult {
    /// Name of the read group.
    pub read_group: String,

    /// Whether or not strandedness was determined for this read group.
    pub succeeded: bool,

    /// The strandedness of this read group or "Inconclusive".
    pub strandedness: String,

    /// The total number of reads in this read group.
    pub total: usize,

    /// The number of reads that are evidence of Forward Strandedness.
    pub forward: usize,

    /// The number of reads that are evidence of Reverse Strandedness.
    pub reverse: usize,

    /// The percent of evidence for Forward Strandedness.
    pub forward_pct: f64,

    /// The percent of evidence for Reverse Strandedness.
    pub reverse_pct: f64,
}

impl ReadGroupDerivedStrandednessResult {
    /// Creates a new [`ReadGroupDerivedStrandednessResult`].
    fn new(
        read_group: String,
        succeeded: bool,
        strandedness: String,
        forward: usize,
        reverse: usize,
    ) -> Self {
        ReadGroupDerivedStrandednessResult {
            read_group,
            succeeded,
            strandedness,
            total: forward + reverse,
            forward,
            reverse,
            forward_pct: (forward as f64 / (forward + reverse) as f64) * 100.0,
            reverse_pct: (reverse as f64 / (forward + reverse) as f64) * 100.0,
        }
    }
}

/// Struct holding the final results for an `ngs derive strandedness` subcommand
/// call.
#[derive(Debug, Serialize)]
pub struct DerivedStrandednessResult {
    /// Whether or not the `ngs derive strandedness` subcommand succeeded.
    pub succeeded: bool,

    /// The strandedness of this read group or "Inconclusive".
    pub strandedness: String,

    /// The total number of reads.
    pub total: usize,

    /// The number of reads that are evidence of Forward Strandedness.
    pub forward: usize,

    /// The number of reads that are evidence of Reverse Strandedness.
    pub reverse: usize,

    /// The percent of evidence for Forward Strandedness.
    pub forward_pct: f64,

    /// The percent of evidence for Reverse Strandedness.
    pub reverse_pct: f64,

    /// Vector of [`ReadGroupDerivedStrandednessResult`]s.
    /// One for each read group in the BAM,
    /// and potentially one for any reads with an unknown read group.
    pub read_groups: Vec<ReadGroupDerivedStrandednessResult>,

    /// General read record metrics that are tallied as a part of the
    /// strandedness subcommand.
    pub read_metrics: ReadRecordMetrics,

    /// General gene metrics that are tallied as a part of the
    /// strandedness subcommand.
    pub gene_metrics: GeneRecordMetrics,

    /// General exon metrics that are tallied as a part of the
    /// strandedness subcommand.
    pub exon_metrics: ExonRecordMetrics,
}

impl DerivedStrandednessResult {
    /// Creates a new [`DerivedStrandednessResult`].
    fn new(
        succeeded: bool,
        strandedness: String,
        forward: usize,
        reverse: usize,
        read_groups: Vec<ReadGroupDerivedStrandednessResult>,
        read_metrics: ReadRecordMetrics,
        gene_metrics: GeneRecordMetrics,
        exon_metrics: ExonRecordMetrics,
    ) -> Self {
        DerivedStrandednessResult {
            succeeded,
            strandedness,
            total: forward + reverse,
            forward,
            reverse,
            forward_pct: (forward as f64 / (forward + reverse) as f64) * 100.0,
            reverse_pct: (reverse as f64 / (forward + reverse) as f64) * 100.0,
            read_groups,
            read_metrics,
            gene_metrics,
            exon_metrics,
        }
    }
}

#[derive(Clone, Copy, Debug)]
enum Strand {
    Forward,
    Reverse,
}

impl From<sam::record::Flags> for Strand {
    fn from(flags: sam::record::Flags) -> Self {
        if flags.is_reverse_complemented() {
            Self::Reverse
        } else {
            Self::Forward
        }
    }
}

impl TryFrom<gff::record::Strand> for Strand {
    type Error = ();

    fn try_from(strand: gff::record::Strand) -> Result<Self, Self::Error> {
        match strand {
            gff::record::Strand::Forward => Ok(Self::Forward),
            gff::record::Strand::Reverse => Ok(Self::Reverse),
            _ => Err(()),
        }
    }
}

#[derive(Clone, Copy, Debug)]
enum SegmentOrder {
    First,
    Last,
}

impl TryFrom<sam::record::Flags> for SegmentOrder {
    type Error = ();

    fn try_from(flags: sam::record::Flags) -> Result<Self, Self::Error> {
        if !flags.is_segmented() {
            Err(())
        } else if flags.is_first_segment() && !flags.is_last_segment() {
            Ok(SegmentOrder::First)
        } else if flags.is_last_segment() && !flags.is_first_segment() {
            Ok(SegmentOrder::Last)
        } else {
            Err(())
        }
    }
}

/// Struct holding the parsed BAM file and its index.
pub struct ParsedBAMFile {
    pub reader: bam::Reader<noodles::bgzf::Reader<std::fs::File>>,
    pub header: sam::Header,
    pub index: bam::bai::Index,
}

/// Filters defining how to calculate strandedness.
pub struct StrandednessFilters {
    /// Minimum number of reads mapped to a gene to be considered
    /// for evidence of strandedness.
    pub min_reads_per_gene: usize,

    /// Minumum mapping quality for a record to be considered.
    /// 0 if MAPQ shouldn't be considered.
    pub min_mapq: u8,

    /// Allow qc failed reads to be counted.
    pub count_qc_failed: bool,

    /// Do not count supplementary alignments.
    pub no_supplementary: bool,

    /// Do count secondary alignments.
    pub count_secondary: bool,

    /// Do count duplicates.
    pub count_duplicates: bool,
}

fn disqualify_gene(
    gene: &gff::Record,
    exons: &HashMap<&str, Lapper<usize, gff::record::Strand>>,
) -> bool {
    let gene_strand = gene.strand();
    if gene_strand != gff::record::Strand::Forward && gene_strand != gff::record::Strand::Reverse {
        return true;
    }
    let mut all_on_same_strand = true;
    let mut at_least_one_exon = false;

    if let Some(intervals) = exons.get(gene.reference_sequence_name()) {
        for exon in intervals.find(gene.start().into(), gene.end().into()) {
            at_least_one_exon = true;
            if exon.val != gene_strand {
                all_on_same_strand = false;
                break;
            }
        }
    }

    if all_on_same_strand && at_least_one_exon {
        return false;
    }
    true
}

fn query_filtered_reads(
    parsed_bam: &ParsedBAMFile,
    gene: &gff::Record,
    filters: &StrandednessFilters,
    read_metrics: &mut ReadRecordMetrics,
) -> Vec<sam::alignment::Record> {
    let start = Position::from(gene.start());
    let end = Position::from(gene.end());
    let region = Region::new(gene.reference_sequence_name(), start..=end);

    let mut filtered_reads = Vec::new();

    let query = parsed_bam
        .reader
        .query(&parsed_bam.header, &parsed_bam.index, &region)
        .unwrap();
    for read in query {
        let read = read.unwrap();

        // (1) Parse the flags so we can see if the read should be discarded.
        let flags = read.flags();
        if (!filters.count_qc_failed && flags.is_qc_fail())
            || (filters.no_supplementary && flags.is_supplementary())
            || (!filters.count_secondary && flags.is_secondary())
            || (!filters.count_duplicates && flags.is_duplicate())
        {
            read_metrics.ignored_flags += 1;
            continue;
        }

        // (2) If the user is filtering by MAPQ, check if this read passes.
        if filters.min_mapq > 0 {
            match read.mapping_quality() {
                Some(mapq) => {
                    if mapq.get() < filters.min_mapq {
                        read_metrics.low_mapq += 1;
                        continue;
                    }
                }
                None => {
                    read_metrics.missing_mapq += 1;
                    continue;
                }
            }
        }

        filtered_reads.push(read);
    }

    if filtered_reads.len() < filters.min_reads_per_gene {
        filtered_reads.clear();
    }

    return filtered_reads;
}

fn classify_read(
    read: &sam::alignment::Record,
    gene_strand: &gff::record::Strand,
    all_counts: &mut HashMap<&str, Counts>,
    read_metrics: &mut ReadRecordMetrics,
) {
    let gene_strand = Strand::try_from(gene_strand).unwrap();

    let read_group = match read.data().get(Tag::ReadGroup) {
        Some(rg) => rg.as_str().unwrap_or_else(|| {
            tracing::warn!("Could not parse a RG tag from a read in the file.");
            UNKNOWN_READ_GROUP.as_str()
        }),
        None => UNKNOWN_READ_GROUP.as_str(),
    };

    let overall_counts = all_counts
        .entry(OVERALL.as_str())
        .or_insert(Counts::default());
    let rg_counts = all_counts.entry(read_group).or_insert(Counts::default());

    let read_strand = Strand::from(read.flags());
    if read.flags().is_segmented() {
        read_metrics.paired_end_reads += 1;

        let order = SegmentOrder::try_from(read.flags()).unwrap();

        match (order, read_strand, gene_strand) {
            (SegmentOrder::First, Strand::Forward, Strand::Forward)
            | (SegmentOrder::First, Strand::Reverse, Strand::Reverse)
            | (SegmentOrder::Last, Strand::Forward, Strand::Reverse)
            | (SegmentOrder::Last, Strand::Reverse, Strand::Forward) => {
                rg_counts.forward += 1;
                overall_counts.forward += 1;
            }
            (SegmentOrder::First, Strand::Forward, Strand::Reverse)
            | (SegmentOrder::First, Strand::Reverse, Strand::Forward)
            | (SegmentOrder::Last, Strand::Forward, Strand::Forward)
            | (SegmentOrder::Last, Strand::Reverse, Strand::Reverse) => {
                rg_counts.reverse += 1;
                overall_counts.reverse += 1;
            }
        }
    } else {
        read_metrics.single_end_reads += 1;

        match (read_strand, gene_strand) {
            (Strand::Forward, Strand::Forward) | (Strand::Reverse, Strand::Reverse) => {
                rg_counts.forward += 1;
                overall_counts.forward += 1;
            }
            (Strand::Forward, Strand::Reverse) | (Strand::Reverse, Strand::Forward) => {
                rg_counts.reverse += 1;
                overall_counts.reverse += 1;
            }
        }
    }
}

/// Method to predict the strandedness of a read group.
fn predict_strandedness(rg_name: &str, counts: &Counts) -> ReadGroupDerivedStrandednessResult {
    if counts.forward == 0 && counts.reverse == 0 {
        return ReadGroupDerivedStrandednessResult {
            read_group: rg_name.to_string(),
            succeeded: false,
            strandedness: "Inconclusive".to_string(),
            total: 0,
            forward: 0,
            reverse: 0,
            forward_pct: 0.0,
            reverse_pct: 0.0,
        };
    }
    let mut result = ReadGroupDerivedStrandednessResult::new(
        rg_name.to_string(),
        false,
        "Inconclusive".to_string(),
        counts.forward,
        counts.reverse,
    );

    if result.forward_pct > STRANDED_THRESHOLD {
        result.succeeded = true;
        result.strandedness = "Forward".to_string();
    } else if result.reverse_pct > STRANDED_THRESHOLD {
        result.succeeded = true;
        result.strandedness = "Reverse".to_string();
    } else if result.forward_pct > UNSTRANDED_THRESHOLD && result.reverse_pct > UNSTRANDED_THRESHOLD
    {
        result.succeeded = true;
        result.strandedness = "Unstranded".to_string();
    }

    return result;
}

/// Main method to evaluate the observed strand state and
/// return a result for the derived strandedness. This may fail, and the
/// resulting [`DerivedStrandednessResult`] should be evaluated accordingly.
pub fn predict(
    parsed_bam: &ParsedBAMFile,
    gene_records: &mut Vec<gff::Record>,
    exons: &HashMap<&str, Lapper<usize, gff::record::Strand>>,
    max_iterations_per_try: usize,
    num_genes: usize,
    filters: &StrandednessFilters,
    gene_metrics: &mut GeneRecordMetrics,
    exon_metrics: &mut ExonRecordMetrics,
    read_metrics: &mut ReadRecordMetrics,
) -> Result<DerivedStrandednessResult, anyhow::Error> {
    let mut rng = rand::thread_rng();
    let mut num_tested_genes: usize = 0; // Local to this attempt
    let mut all_counts: HashMap<&str, Counts> = HashMap::new();

    all_counts.insert(UNKNOWN_READ_GROUP.as_str(), Counts::default());
    all_counts.insert(OVERALL.as_str(), Counts::default());

    for _ in 0..max_iterations_per_try {
        if num_tested_genes > num_genes {
            break;
        }

        let cur_gene = gene_records.swap_remove(rng.gen_range(0..gene_records.len()));

        if disqualify_gene(&cur_gene, exons) {
            gene_metrics.bad_strands += 1;
            continue;
        }

        let mut enough_reads = false;
        for read in query_filtered_reads(parsed_bam, &cur_gene, filters, read_metrics) {
            enough_reads = true;

            classify_read(&read, &cur_gene.strand(), &mut all_counts, read_metrics);
        }
        if enough_reads {
            num_tested_genes += 1;
        } else {
            gene_metrics.not_enough_reads += 1;
        }
    }

    gene_metrics.tested += num_tested_genes; // Add to any other attempts

    // Overly complicated but IDK how to simplify this
    let found_rgs = all_counts
        .keys()
        .cloned()
        .map(|rg| rg.to_string())
        .collect::<Vec<_>>();
    let found_rgs_arc = found_rgs
        .iter()
        .map(|rg| Arc::new(rg.clone()))
        .collect::<HashSet<_>>();

    let rgs_in_header_not_found = validate_read_group_info(&found_rgs_arc, &parsed_bam.header);
    for rg in rgs_in_header_not_found {
        all_counts.insert(rg.as_str(), Counts::default());
    }

    let mut final_result = DerivedStrandednessResult::new(
        true,
        "Inconclusive".to_string(),
        0,
        0,
        Vec::new(),
        read_metrics.clone(),
        gene_metrics.clone(),
        exon_metrics.clone(),
    );

    for (rg, counts) in all_counts {
        if rg == UNKNOWN_READ_GROUP.as_str() && counts.forward == 0 && counts.reverse == 0 {
            continue;
        }
    }

    anyhow::Ok(final_result)
}
