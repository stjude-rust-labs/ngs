//! Module holding the logic for annotating junctions.

use anyhow::bail;
use noodles::sam::alignment::Record;
use noodles::sam::record::cigar::op::Kind;
use noodles::sam::Header;
use std::collections::HashMap;
use std::num::NonZeroUsize;

use crate::derive::junction_annotation::results::JunctionAnnotationResults;

/// Parameters defining how to annotate found junctions
pub struct JunctionAnnotationParameters {
    /// Minimum intron length to consider.
    pub min_intron_length: usize,

    /// Add +- this amount to intron positions when looking up exon positions.
    pub fuzzy_junction_match_range: u8,

    /// Minimum number of reads supporting a junction to be considered.
    pub min_read_support: u8,

    /// Minumum mapping quality for a record to be considered.
    /// 0 if MAPQ shouldn't be considered.
    pub min_mapq: u8,

    /// Do not count supplementary alignments.
    pub no_supplementary: bool,

    /// Do count secondary alignments.
    pub count_secondary: bool,

    /// Do count duplicates.
    pub count_duplicates: bool,
}

/// Main function to annotate junctions one record at a time.
pub fn process(
    record: &Record,
    exon_starts: &HashMap<&str, Vec<usize>>,
    exon_ends: &HashMap<&str, Vec<usize>>,
    header: &Header,
    params: &JunctionAnnotationParameters,
    results: &mut JunctionAnnotationResults,
) -> anyhow::Result<()> {
    // (1) Parse the read name.
    let read_name = match record.read_name() {
        Some(name) => name,
        _ => bail!("Could not parse read name"),
    };

    // (2) Parse the flags so we can see if the read should be ignored.
    let flags = record.flags();

    if flags.is_unmapped()
        || (params.no_supplementary && flags.is_supplementary())
        || (!params.count_secondary && flags.is_secondary())
        || (!params.count_duplicates && flags.is_duplicate())
    {
        results.records.ignored_flags += 1;
        return Ok(());
    }

    // (3) Parse the CIGAR string from the record.
    // We only care about reads with introns, so if there are no introns
    // we can skip this read.
    let cigar = record.cigar();
    if !cigar.iter().any(|op| matches!(op.kind(), Kind::Skip)) {
        results.records.not_spliced += 1;
        return Ok(());
    }

    // (4) If the user is filtering by MAPQ, check if this read passes.
    // Log if the read is filtered out for a too low MAPQ or a missing MAPQ.
    if params.min_mapq > 0 {
        match record.mapping_quality() {
            Some(mapq) => {
                if mapq.get() < params.min_mapq {
                    results.records.low_mapq += 1;
                    return Ok(());
                }
            }
            None => {
                results.records.missing_mapq += 1;
                return Ok(());
            }
        }
    }

    // (5) Parse the reference sequence id from the record.
    let id = match record.reference_sequence_id() {
        Some(id) => id,
        _ => {
            bail!(
                "Could not parse reference sequence id for read: {}",
                read_name
            )
        }
    };

    // (6) Map the parsed reference sequence id to a reference sequence name.
    let seq_name = match header
        .reference_sequences()
        .get_index(id)
        .map(|(name, _)| Some(name))
    {
        Some(Some(name)) => name.as_str(),
        _ => {
            bail!(
                "Could not map reference sequence id to header for read: {}",
                read_name
            )
        }
    };

    // (7) Check if there will be annotations for this reference sequence.
    let mut ref_is_annotated = true;
    if !exon_starts.contains_key(&seq_name) || !exon_ends.contains_key(&seq_name) {
        ref_is_annotated = false;
    }

    // (8) Calculate the start position of this read. This will
    // later be used to find the position of any introns.
    let start = match record.alignment_start() {
        Some(s) => usize::from(s),
        _ => bail!("Could not parse record's start position."),
    };

    // (9) Find introns
    let mut cur_pos = start;
    for op in cigar.iter() {
        match op.kind() {
            // Operations that increment the reference position.
            Kind::Match | Kind::Deletion | Kind::SequenceMatch | Kind::SequenceMismatch => {
                cur_pos += op.len();
            }
            // This is an intron.
            Kind::Skip => {
                // Do this check later, for better metric reporting.
                // if op.len() < params.min_intron_length {
                //     continue;
                // }

                let intron_start = cur_pos;
                let intron_end = cur_pos + op.len();
                // Update cur_pos to the end of the intron
                // in case there are multiple introns in the read.
                cur_pos = intron_end;

                // If the reference sequence is not annotated, we can skip
                // the lookup of exon positions, and directly insert the
                // intron into the unannotated_reference HashMap.
                if !ref_is_annotated {
                    results
                        .junction_annotations
                        .unannotated_reference
                        .entry(seq_name.to_string())
                        .or_default()
                        .entry((
                            NonZeroUsize::new(intron_start).unwrap(),
                            NonZeroUsize::new(intron_end).unwrap(),
                        ))
                        .and_modify(|e| *e += 1)
                        .or_insert(1);
                    continue;
                }

                let exon_starts = match exon_starts.get(&seq_name) {
                    Some(starts) => starts,
                    _ => bail!("Could not find exon starts for contig: {}", seq_name),
                };
                let exon_ends = match exon_ends.get(&seq_name) {
                    Some(ends) => ends,
                    _ => bail!("Could not find exon ends for contig: {}", seq_name),
                };

                let mut intron_start_known = false;
                let mut intron_end_known = false;
                // To allow collapsing fuzzy junctions,
                // we need to store the reference positions of the exon boundaries.
                // We initialize these values to the position of the found intron.
                let mut ref_intron_start = intron_start;
                let mut ref_intron_end = intron_end;
                for exon_end in exon_ends.iter() {
                    if intron_start >= (exon_end - params.fuzzy_junction_match_range as usize)
                        && intron_start <= (exon_end + params.fuzzy_junction_match_range as usize)
                    {
                        intron_start_known = true;
                        ref_intron_start = *exon_end;
                        break;
                    }
                }
                for exon_start in exon_starts.iter() {
                    if intron_end >= (exon_start - params.fuzzy_junction_match_range as usize)
                        && intron_end <= (exon_start + params.fuzzy_junction_match_range as usize)
                    {
                        intron_end_known = true;
                        ref_intron_end = *exon_start;
                        break;
                    }
                }

                match (intron_start_known, intron_end_known) {
                    (true, true) => {
                        // We found both ends of the intron.
                        // This is a Known Junction.
                        results
                            .junction_annotations
                            .known
                            .entry(seq_name.to_string())
                            .or_default()
                            .entry((
                                NonZeroUsize::new(ref_intron_start).unwrap(),
                                NonZeroUsize::new(ref_intron_end).unwrap(),
                            ))
                            .and_modify(|e| *e += 1)
                            .or_insert(1);
                    }
                    (true, false) | (false, true) => {
                        // We found one end of the intron,
                        // but not the other.
                        // This is a Partial Novel Junction.
                        results
                            .junction_annotations
                            .partial_novel
                            .entry(seq_name.to_string())
                            .or_default()
                            .entry((
                                NonZeroUsize::new(ref_intron_start).unwrap(),
                                NonZeroUsize::new(ref_intron_end).unwrap(),
                            ))
                            .and_modify(|e| *e += 1)
                            .or_insert(1);
                    }
                    (false, false) => {
                        // We found neither end of the intron.
                        // This is a Complete Novel Junction.
                        results
                            .junction_annotations
                            .complete_novel
                            .entry(seq_name.to_string())
                            .or_default()
                            .entry((
                                NonZeroUsize::new(ref_intron_start).unwrap(),
                                NonZeroUsize::new(ref_intron_end).unwrap(),
                            ))
                            .and_modify(|e| *e += 1)
                            .or_insert(1);
                    }
                }
            }
            // Operations that do not increment the reference position.
            _ => {}
        }
    }

    results.records.processed += 1;
    Ok(())
}

/// Main function to summarize the results of the junction_annotation subcommand.
pub fn summarize(
    results: &mut JunctionAnnotationResults,
    params: &JunctionAnnotationParameters,
) -> anyhow::Result<()> {
    // Filter out junctions that are too short or don't have enough read support.
    let mut num_junctions_too_short: usize = 0;
    let mut num_not_enough_support: usize = 0;
    for (_, v) in results.junction_annotations.known.iter_mut() {
        v.retain(|(start, end), count| {
            if end.get() - start.get() < params.min_intron_length {
                num_junctions_too_short += 1;
                false
            } else if *count < params.min_read_support as usize {
                num_not_enough_support += 1;
                false
            } else {
                true
            }
        });
    }
    for (_, v) in results.junction_annotations.partial_novel.iter_mut() {
        v.retain(|(start, end), count| {
            if end.get() - start.get() < params.min_intron_length {
                num_junctions_too_short += 1;
                false
            } else if *count < params.min_read_support as usize {
                num_not_enough_support += 1;
                false
            } else {
                true
            }
        });
    }
    for (_, v) in results.junction_annotations.complete_novel.iter_mut() {
        v.retain(|(start, end), count| {
            if end.get() - start.get() < params.min_intron_length {
                num_junctions_too_short += 1;
                false
            } else if *count < params.min_read_support as usize {
                num_not_enough_support += 1;
                false
            } else {
                true
            }
        });
    }
    for (_, v) in results
        .junction_annotations
        .unannotated_reference
        .iter_mut()
    {
        v.retain(|(start, end), count| {
            if end.get() - start.get() < params.min_intron_length {
                num_junctions_too_short += 1;
                false
            } else if *count < params.min_read_support as usize {
                num_not_enough_support += 1;
                false
            } else {
                true
            }
        });
    }
    results.summary.intron_too_short = num_junctions_too_short;
    results.summary.junctions_with_not_enough_read_support = num_not_enough_support;

    // Tally up observed junctions and spliced reads.
    results.summary.known_junctions = results
        .junction_annotations
        .known
        .values()
        .map(|v| v.len())
        .sum();
    results.summary.known_spliced_reads = results
        .junction_annotations
        .known
        .values()
        .map(|v| v.values().sum::<usize>())
        .sum();
    results.summary.partial_novel_junctions = results
        .junction_annotations
        .partial_novel
        .values()
        .map(|v| v.len())
        .sum();
    results.summary.partial_novel_spliced_reads = results
        .junction_annotations
        .partial_novel
        .values()
        .map(|v| v.values().sum::<usize>())
        .sum();
    results.summary.complete_novel_junctions = results
        .junction_annotations
        .complete_novel
        .values()
        .map(|v| v.len())
        .sum();
    results.summary.complete_novel_spliced_reads = results
        .junction_annotations
        .complete_novel
        .values()
        .map(|v| v.values().sum::<usize>())
        .sum();
    results.summary.unannotated_reference_junctions = results
        .junction_annotations
        .unannotated_reference
        .values()
        .map(|v| v.len())
        .sum();
    results.summary.unannotated_reference_spliced_reads = results
        .junction_annotations
        .unannotated_reference
        .values()
        .map(|v| v.values().sum::<usize>())
        .sum();

    // Tally up total junctions and spliced reads.
    results.summary.total_junctions = results.summary.known_junctions
        + results.summary.partial_novel_junctions
        + results.summary.complete_novel_junctions
        + results.summary.unannotated_reference_junctions;
    results.summary.total_spliced_reads = results.summary.known_spliced_reads
        + results.summary.partial_novel_spliced_reads
        + results.summary.complete_novel_spliced_reads
        + results.summary.unannotated_reference_spliced_reads;

    Ok(())
}
