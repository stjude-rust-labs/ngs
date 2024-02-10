//! Module holding the logic for computing the strandedness.

use noodles::bam;
use noodles::core::Region;
use noodles::gff;
use noodles::sam;
use noodles::sam::record::data::field::Tag;
use noodles::sam::record::MappingQuality;
use rand::Rng;
use rust_lapper::Lapper;
use std::collections::HashMap;
use std::collections::HashSet;
use std::sync::Arc;

use crate::derive::strandedness::results;
use crate::utils::alignment::filter_by_mapq;
use crate::utils::display::RecordCounter;
use crate::utils::read_groups::{validate_read_group_info, UNKNOWN_READ_GROUP};

const STRANDED_THRESHOLD: f64 = 80.0;
const UNSTRANDED_THRESHOLD: f64 = 40.0;

/// Struct for tracking count results.
#[derive(Clone, Default)]
pub struct Counts {
    /// The number of reads that are evidence of Forward Strandedness.
    forward: usize,

    /// The number of reads that are evidence of Reverse Strandedness.
    reverse: usize,
}

/// Struct for tracking possible strand orientations.
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

/// Struct for tracking the order of segments in a record.
#[derive(Clone, Copy, Debug)]
enum SegmentOrder {
    First,
    Last,
}

impl TryFrom<sam::record::Flags> for SegmentOrder {
    type Error = String;

    fn try_from(flags: sam::record::Flags) -> Result<Self, Self::Error> {
        if !flags.is_segmented() {
            Err("Expected segmented record.".to_string())
        } else if flags.is_first_segment() && !flags.is_last_segment() {
            Ok(SegmentOrder::First)
        } else if flags.is_last_segment() && !flags.is_first_segment() {
            Ok(SegmentOrder::Last)
        } else {
            Err("Expected first or last segment.".to_string())
        }
    }
}

/// Struct holding the parsed BAM file and its index.
pub struct ParsedBAMFile {
    /// The BAM reader.
    pub reader: bam::Reader<noodles::bgzf::Reader<std::fs::File>>,

    /// The BAM header.
    pub header: sam::Header,

    /// The BAM index.
    pub index: bam::bai::Index,
}

/// Struct holding the counts for all read groups.
/// Also holds the set of read groups found in the BAM.
pub struct AllReadGroupsCounts {
    /// The counts for all read groups.
    pub counts: HashMap<Arc<String>, Counts>,

    /// The set of read groups found in the BAM.
    pub found_rgs: HashSet<Arc<String>>,
}

/// Parameters defining how to calculate strandedness.
pub struct StrandednessParams {
    /// The number of genes to test for strandedness.
    pub num_genes: usize,

    /// The maximum number of iterations to try before giving up.
    pub max_genes_per_try: usize,

    /// Minimum number of reads mapped to a gene to be considered
    /// for evidence of strandedness.
    pub min_reads_per_gene: usize,

    /// Minumum mapping quality for a record to be considered.
    /// 0 if MAPQ shouldn't be considered.
    pub min_mapq: Option<MappingQuality>,

    /// Allow qc failed reads to be counted.
    pub count_qc_failed: bool,

    /// Do not count supplementary alignments.
    pub no_supplementary: bool,

    /// Do count secondary alignments.
    pub count_secondary: bool,

    /// Do count duplicates.
    pub count_duplicates: bool,
}

/// Function to disqualify a gene based on its strand and exons.
fn disqualify_gene(
    gene: &gff::Record,
    exons: &HashMap<&str, Lapper<usize, gff::record::Strand>>,
) -> bool {
    // gene_strand guaranteed to be Forward or Reverse by initialization code.
    let gene_strand = gene.strand();
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

/// Function to filter out records based on their flags.
fn filter_by_flags(record: &sam::alignment::Record, params: &StrandednessParams) -> bool {
    let flags = record.flags();
    if (!params.count_qc_failed && flags.is_qc_fail())
        || (params.no_supplementary && flags.is_supplementary())
        || (!params.count_secondary && flags.is_secondary())
        || (!params.count_duplicates && flags.is_duplicate())
    {
        return true;
    }
    false
}

/// Function to query the BAM file and filter the records based on the
/// parameters provided.
fn query_and_filter(
    parsed_bam: &mut ParsedBAMFile,
    gene: &gff::Record,
    params: &StrandednessParams,
    read_metrics: &mut results::ReadRecordMetrics,
) -> Vec<sam::alignment::Record> {
    let start = gene.start();
    let end = gene.end();
    let region = Region::new(gene.reference_sequence_name(), start..=end);

    let mut filtered_reads = Vec::new();

    let query = parsed_bam
        .reader
        .query(&parsed_bam.header, &parsed_bam.index, &region)
        .unwrap();
    for read in query {
        let read = read.unwrap();

        // (1) Filter by flags.
        if filter_by_flags(&read, params) {
            read_metrics.filtered_by_flags += 1;
            continue;
        }

        // (2) Filter by MAPQ.
        if filter_by_mapq(&read, params.min_mapq) {
            read_metrics.bad_mapq += 1;
            continue;
        }

        filtered_reads.push(read);
    }

    if filtered_reads.len() < params.min_reads_per_gene {
        filtered_reads.clear();
    }

    filtered_reads
}

/// Function to classify a read based on its strand and the strand of the gene.
fn classify_read(
    read: &sam::alignment::Record,
    gene_strand: &Strand,
    all_counts: &mut AllReadGroupsCounts,
    read_metrics: &mut results::ReadRecordMetrics,
) {
    let read_group = match read.data().get(Tag::ReadGroup) {
        Some(rg) => {
            let rg = rg.to_string();
            if !all_counts.found_rgs.contains(&rg) {
                all_counts.found_rgs.insert(Arc::new(rg.clone()));
            }
            Arc::clone(all_counts.found_rgs.get(&rg).unwrap())
        }
        None => Arc::clone(&UNKNOWN_READ_GROUP),
    };

    let rg_counts = all_counts.counts.entry(read_group).or_default();

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
            }
            (SegmentOrder::First, Strand::Forward, Strand::Reverse)
            | (SegmentOrder::First, Strand::Reverse, Strand::Forward)
            | (SegmentOrder::Last, Strand::Forward, Strand::Forward)
            | (SegmentOrder::Last, Strand::Reverse, Strand::Reverse) => {
                rg_counts.reverse += 1;
            }
        }
    } else {
        read_metrics.single_end_reads += 1;

        match (read_strand, gene_strand) {
            (Strand::Forward, Strand::Forward) | (Strand::Reverse, Strand::Reverse) => {
                rg_counts.forward += 1;
            }
            (Strand::Forward, Strand::Reverse) | (Strand::Reverse, Strand::Forward) => {
                rg_counts.reverse += 1;
            }
        }
    }
}

/// Method to predict the strandedness of a read group.
fn predict_strandedness(
    rg_name: &str,
    counts: &Counts,
) -> results::ReadGroupDerivedStrandednessResult {
    if counts.forward == 0 && counts.reverse == 0 {
        return results::ReadGroupDerivedStrandednessResult {
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
    let mut result = results::ReadGroupDerivedStrandednessResult::new(
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

    result
}

/// Main method to evaluate the observed strand state and
/// return a result for the derived strandedness. This may fail, and the
/// resulting [`DerivedStrandednessResult`] should be evaluated accordingly.
pub fn predict(
    parsed_bam: &mut ParsedBAMFile,
    gene_records: &mut Vec<gff::Record>,
    exons: &HashMap<&str, Lapper<usize, gff::record::Strand>>,
    all_counts: &mut AllReadGroupsCounts,
    params: &StrandednessParams,
    metrics: &mut results::RecordTracker,
) -> Result<results::DerivedStrandednessResult, anyhow::Error> {
    let mut rng = rand::thread_rng();
    let mut num_genes_considered: usize = 0; // Local to this attempt
    let mut counter = RecordCounter::new(Some(1_000));
    let genes_remaining = gene_records.len();

    let max_iters = if params.max_genes_per_try > genes_remaining {
        tracing::warn!(
            "The number of genes remaining ({}) is less than the maximum iterations per try ({}).",
            genes_remaining,
            params.max_genes_per_try,
        );
        genes_remaining
    } else {
        params.max_genes_per_try
    };

    for _ in 0..max_iters {
        if num_genes_considered >= params.num_genes {
            tracing::info!(
                "Reached the maximum number of considered genes ({}) for this try.",
                num_genes_considered,
            );
            break;
        }

        let cur_gene = gene_records.swap_remove(rng.gen_range(0..gene_records.len()));
        counter.inc();

        if disqualify_gene(&cur_gene, exons) {
            metrics.genes.mixed_strands += 1; // Tracked across attempts
            continue;
        }
        let cur_gene_strand = Strand::try_from(cur_gene.strand()).unwrap();

        let mut enough_reads = false;
        for read in query_and_filter(parsed_bam, &cur_gene, params, &mut metrics.reads) {
            enough_reads = true;

            classify_read(&read, &cur_gene_strand, all_counts, &mut metrics.reads);
        }
        if enough_reads {
            num_genes_considered += 1;
        } else {
            metrics.genes.not_enough_reads += 1; // Tracked across attempts
        }
    }
    if num_genes_considered < params.num_genes {
        tracing::warn!(
            "Reached the maximum number of iterations ({}) before considering the requested amount of genes ({}) for this try. Only considering {} genes.",
            max_iters,
            params.num_genes,
            num_genes_considered,
        );
    }

    metrics.genes.considered += num_genes_considered; // Add to any other attempts
    metrics.genes.evaluated += counter.get(); // Add to any other attempts

    // TODO: Should this be done in derive()? Will re-run for each attempt.
    // Might cause false positives?
    let rgs_in_header_not_found =
        validate_read_group_info(&all_counts.found_rgs, &parsed_bam.header);
    for rg in rgs_in_header_not_found {
        all_counts
            .counts
            .insert(Arc::new(rg.to_string()), Counts::default());
    }

    let mut overall_counts = Counts::default();
    let mut rg_results = Vec::new();
    for (rg, counts) in &all_counts.counts {
        overall_counts.forward += counts.forward;
        overall_counts.reverse += counts.reverse;

        let result = predict_strandedness(rg, counts);
        rg_results.push(result)
    }

    let overall_result = predict_strandedness("overall", &overall_counts);
    let final_result = results::DerivedStrandednessResult::new(
        overall_result.succeeded,
        overall_result.strandedness,
        overall_result.forward,
        overall_result.reverse,
        rg_results,
        metrics.clone(),
    );

    anyhow::Ok(final_result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_lapper::Interval;

    #[test]
    fn test_disqualify_gene() {
        let mut exons = HashMap::new();
        exons.insert(
            "chr1",
            Lapper::new(vec![
                Interval {
                    start: 1,
                    stop: 10,
                    val: gff::record::Strand::Forward,
                },
                Interval {
                    start: 11,
                    stop: 20,
                    val: gff::record::Strand::Reverse,
                },
            ]),
        );

        let gene = gff::Record::default();
        assert!(disqualify_gene(&gene, &exons));

        let mut exons = HashMap::new();
        exons.insert(
            "chr1",
            Lapper::new(vec![
                Interval {
                    start: 1,
                    stop: 10,
                    val: gff::record::Strand::Forward,
                },
                Interval {
                    start: 11,
                    stop: 20,
                    val: gff::record::Strand::Forward,
                },
            ]),
        );

        let s = "chr1\tNOODLES\tgene\t8\t13\t.\t+\t.\tgene_id=ndls0;gene_name=gene0";
        let record = s.parse::<gff::Record>().unwrap();
        assert!(!disqualify_gene(&record, &exons));
    }

    #[test]
    fn test_query_and_filter() { // TODO
    }

    #[test]
    fn test_classify_read() {
        // Set up
        let mut all_counts = AllReadGroupsCounts {
            counts: HashMap::new(),
            found_rgs: HashSet::new(),
        };
        let mut read_metrics = results::ReadRecordMetrics::default();
        let counts_key = Arc::new("rg1".to_string());
        let rg_tag = sam::record::data::field::Value::String("rg1".to_string());

        // Test Single-End read. Evidence for Forward Strandedness.
        let mut read = sam::alignment::Record::default();
        read.flags_mut().set(0x1.into(), false);
        read.data_mut().insert(Tag::ReadGroup, rg_tag.clone());
        classify_read(&read, &Strand::Forward, &mut all_counts, &mut read_metrics);
        assert_eq!(read_metrics.paired_end_reads, 0);
        assert_eq!(read_metrics.single_end_reads, 1);
        assert_eq!(read_metrics.filtered_by_flags, 0);
        assert_eq!(read_metrics.bad_mapq, 0);
        let counts = all_counts.counts.get(&counts_key).unwrap();
        assert_eq!(counts.forward, 1);
        assert_eq!(counts.reverse, 0);

        // Test Paired-End read. Evidence for Forward Strandedness.
        let mut read = sam::alignment::Record::default();
        read.flags_mut().set(0x1.into(), true);
        read.flags_mut().set(0x40.into(), true);
        read.data_mut().insert(Tag::ReadGroup, rg_tag.clone());
        classify_read(&read, &Strand::Forward, &mut all_counts, &mut read_metrics);
        assert_eq!(read_metrics.paired_end_reads, 1);
        assert_eq!(read_metrics.single_end_reads, 1);
        assert_eq!(read_metrics.filtered_by_flags, 0);
        assert_eq!(read_metrics.bad_mapq, 0);
        let counts = all_counts.counts.get(&counts_key).unwrap();
        assert_eq!(counts.forward, 2);
        assert_eq!(counts.reverse, 0);

        // Test Paired-End read. Evidence for Forward Strandedness.
        let mut read = sam::alignment::Record::default();
        read.flags_mut().set(0x1.into(), true);
        read.flags_mut().set(0x80.into(), true);
        read.data_mut().insert(Tag::ReadGroup, rg_tag.clone());
        classify_read(&read, &Strand::Reverse, &mut all_counts, &mut read_metrics);
        assert_eq!(read_metrics.paired_end_reads, 2);
        assert_eq!(read_metrics.single_end_reads, 1);
        assert_eq!(read_metrics.filtered_by_flags, 0);
        assert_eq!(read_metrics.bad_mapq, 0);
        let counts = all_counts.counts.get(&counts_key).unwrap();
        assert_eq!(counts.forward, 3);
        assert_eq!(counts.reverse, 0);

        // Test Paired-End read. Evidence for Reverse Strandedness.
        let mut read = sam::alignment::Record::default();
        read.flags_mut().set(0x1.into(), true);
        read.flags_mut().set(0x40.into(), true);
        read.data_mut().insert(Tag::ReadGroup, rg_tag.clone());
        classify_read(&read, &Strand::Reverse, &mut all_counts, &mut read_metrics);
        assert_eq!(read_metrics.paired_end_reads, 3);
        assert_eq!(read_metrics.single_end_reads, 1);
        assert_eq!(read_metrics.filtered_by_flags, 0);
        assert_eq!(read_metrics.bad_mapq, 0);
        let counts = all_counts.counts.get(&counts_key).unwrap();
        assert_eq!(counts.forward, 3);
        assert_eq!(counts.reverse, 1);
    }

    #[test]
    fn test_predict_strandedness() {
        let counts = Counts {
            forward: 10,
            reverse: 90,
        };
        let result = predict_strandedness("rg1", &counts);
        assert!(result.succeeded);
        assert_eq!(result.strandedness, "Reverse");
        assert_eq!(result.forward, 10);
        assert_eq!(result.reverse, 90);
        assert_eq!(result.forward_pct, 10.0);
        assert_eq!(result.reverse_pct, 90.0);

        let counts = Counts {
            forward: 50,
            reverse: 50,
        };
        let result = predict_strandedness("rg1", &counts);
        assert!(result.succeeded);
        assert_eq!(result.strandedness, "Unstranded");
        assert_eq!(result.forward, 50);
        assert_eq!(result.reverse, 50);
        assert_eq!(result.forward_pct, 50.0);
        assert_eq!(result.reverse_pct, 50.0);

        let counts = Counts {
            forward: 90,
            reverse: 10,
        };
        let result = predict_strandedness("rg1", &counts);
        assert!(result.succeeded);
        assert_eq!(result.strandedness, "Forward");
        assert_eq!(result.forward, 90);
        assert_eq!(result.reverse, 10);
        assert_eq!(result.forward_pct, 90.0);
        assert_eq!(result.reverse_pct, 10.0);

        let counts = Counts {
            forward: 0,
            reverse: 0,
        };
        let result = predict_strandedness("rg1", &counts);
        assert!(!result.succeeded);
        assert_eq!(result.strandedness, "Inconclusive");
        assert_eq!(result.forward, 0);
        assert_eq!(result.reverse, 0);
        assert_eq!(result.forward_pct, 0.0);
        assert_eq!(result.reverse_pct, 0.0);
    }
}
