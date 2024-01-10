//! Module holding the logic for annotating junctions.

use anyhow::bail;
use anyhow::Ok;
use noodles::sam::alignment::Record;
use noodles::sam::record::cigar::op::Kind;
use noodles::sam::Header;
use std::collections::HashMap;
use std::collections::HashSet;
use std::num::NonZeroUsize;

use crate::derive::junction_annotation::results::JunctionAnnotationResults;

/// Parameters defining how to annotate found junctions
pub struct JunctionAnnotationParameters {
    /// Minimum intron length to consider.
    pub min_intron_length: usize,

    /// Minimum number of reads supporting a junction to be considered.
    pub min_read_support: usize,

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
    exon_starts: &HashMap<&str, HashSet<usize>>,
    exon_ends: &HashMap<&str, HashSet<usize>>,
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
                if exon_ends.contains(&intron_start) {
                    intron_start_known = true;
                }
                if exon_starts.contains(&intron_end) {
                    intron_end_known = true;
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
                                NonZeroUsize::new(intron_start).unwrap(),
                                NonZeroUsize::new(intron_end).unwrap(),
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
                                NonZeroUsize::new(intron_start).unwrap(),
                                NonZeroUsize::new(intron_end).unwrap(),
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
                                NonZeroUsize::new(intron_start).unwrap(),
                                NonZeroUsize::new(intron_end).unwrap(),
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
pub fn summarize(results: &mut JunctionAnnotationResults, params: &JunctionAnnotationParameters) {
    // Filter out junctions that are too short or don't have enough read support.
    let mut num_rejected: usize = 0;
    let mut num_junctions_too_short: usize = 0;
    let mut num_not_enough_support: usize = 0;
    for (_, v) in results.junction_annotations.known.iter_mut() {
        v.retain(|(start, end), count| {
            if end.get() - start.get() < params.min_intron_length {
                num_junctions_too_short += 1;
                if *count < params.min_read_support {
                    num_not_enough_support += 1;
                }
                num_rejected += 1;
                false
            } else if *count < params.min_read_support {
                num_not_enough_support += 1;
                if end.get() - start.get() < params.min_intron_length {
                    num_junctions_too_short += 1;
                }
                num_rejected += 1;
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
                if *count < params.min_read_support {
                    num_not_enough_support += 1;
                }
                num_rejected += 1;
                false
            } else if *count < params.min_read_support {
                num_not_enough_support += 1;
                if end.get() - start.get() < params.min_intron_length {
                    num_junctions_too_short += 1;
                }
                num_rejected += 1;
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
                if *count < params.min_read_support {
                    num_not_enough_support += 1;
                }
                num_rejected += 1;
                false
            } else if *count < params.min_read_support {
                num_not_enough_support += 1;
                if end.get() - start.get() < params.min_intron_length {
                    num_junctions_too_short += 1;
                }
                num_rejected += 1;
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
                if *count < params.min_read_support {
                    num_not_enough_support += 1;
                }
                num_rejected += 1;
                false
            } else if *count < params.min_read_support {
                num_not_enough_support += 1;
                if end.get() - start.get() < params.min_intron_length {
                    num_junctions_too_short += 1;
                }
                num_rejected += 1;
                false
            } else {
                true
            }
        });
    }
    results.summary.total_rejected_junctions = num_rejected;
    results.summary.intron_too_short = num_junctions_too_short;
    results.summary.junctions_with_not_enough_read_support = num_not_enough_support;

    // Tally up observed junctions and spliced reads.
    results.summary.known_junctions = results
        .junction_annotations
        .known
        .values()
        .map(|v| v.len())
        .sum();
    results.summary.known_junctions_read_support = results
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
    results.summary.partial_novel_junctions_read_support = results
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
    results.summary.complete_novel_junctions_read_support = results
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
    results.summary.unannotated_reference_junctions_read_support = results
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
    results.summary.total_junctions_read_support = results.summary.known_junctions_read_support
        + results.summary.partial_novel_junctions_read_support
        + results.summary.complete_novel_junctions_read_support
        + results.summary.unannotated_reference_junctions_read_support;

    // Calculate percentages.
    let total_junctions = results.summary.total_junctions as f64
        - results.summary.unannotated_reference_junctions as f64;
    results.summary.known_junctions_percent =
        results.summary.known_junctions as f64 / total_junctions * 100.0;
    results.summary.partial_novel_junctions_percent =
        results.summary.partial_novel_junctions as f64 / total_junctions * 100.0;
    results.summary.complete_novel_junctions_percent =
        results.summary.complete_novel_junctions as f64 / total_junctions * 100.0;

    // Calculate average read support.
    results.summary.average_junction_read_support = results.summary.total_junctions_read_support
        as f64
        / results.summary.total_junctions as f64;
    results.summary.average_known_junction_read_support =
        results.summary.known_junctions_read_support as f64
            / results.summary.known_junctions as f64;
    results.summary.average_partial_novel_junction_read_support =
        results.summary.partial_novel_junctions_read_support as f64
            / results.summary.partial_novel_junctions as f64;
    results.summary.average_complete_novel_junction_read_support =
        results.summary.complete_novel_junctions_read_support as f64
            / results.summary.complete_novel_junctions as f64;
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles::core::Position;
    use noodles::sam::header::record::value::map;
    use noodles::sam::header::record::value::map::header::Version;
    use noodles::sam::header::record::value::map::{Map, ReferenceSequence};
    use noodles::sam::record::MappingQuality;
    use noodles::sam::record::ReadName;

    #[test]
    fn test_process_and_summarize() {
        // Setup
        let mut results = JunctionAnnotationResults::default();
        let params = JunctionAnnotationParameters {
            min_intron_length: 10,
            min_read_support: 2,
            min_mapq: 30,
            no_supplementary: false,
            count_secondary: false,
            count_duplicates: false,
        };
        let header = Header::builder()
            .set_header(Map::<map::Header>::new(Version::new(1, 6)))
            .add_reference_sequence(
                "sq1".parse().unwrap(),
                Map::<ReferenceSequence>::new(NonZeroUsize::try_from(800).unwrap()),
            )
            .add_reference_sequence(
                "sq1_random".parse().unwrap(), // unannotated
                Map::<ReferenceSequence>::new(NonZeroUsize::try_from(400).unwrap()),
            )
            .build();
        let exon_starts: HashMap<&str, HashSet<usize>> =
            HashMap::from([("sq1", HashSet::from([1, 11, 21, 31, 41, 51, 61, 71]))]);
        let exon_ends = exon_starts
            .iter()
            .map(|(k, v)| (*k, v.iter().map(|e| e + 10).collect()))
            .collect::<HashMap<&str, HashSet<usize>>>();

        // Test known junction
        let mut record = Record::default();
        let r1_name: ReadName = "known1".parse().unwrap();
        *record.read_name_mut() = Some(r1_name);
        *record.flags_mut() = 0.into();
        *record.reference_sequence_id_mut() = Some(0);
        *record.alignment_start_mut() = Position::new(1);
        *record.cigar_mut() = "10M10N10M".parse().unwrap();
        *record.mapping_quality_mut() = MappingQuality::new(60);
        process(
            &record,
            &exon_starts,
            &exon_ends,
            &header,
            &params,
            &mut results,
        )
        .unwrap();
        assert_eq!(results.records.processed, 1);
        assert_eq!(results.records.ignored_flags, 0);
        assert_eq!(results.records.not_spliced, 0);
        assert_eq!(results.records.low_mapq, 0);
        assert_eq!(results.records.missing_mapq, 0);

        // Test that unmapped gets ignored
        let mut record = Record::default();
        let r2_name: ReadName = "unmapped".parse().unwrap();
        *record.read_name_mut() = Some(r2_name);
        *record.flags_mut() = 0x4.into();
        *record.reference_sequence_id_mut() = Some(0);
        *record.alignment_start_mut() = Position::new(1);
        *record.cigar_mut() = "10M10N10M".parse().unwrap();
        *record.mapping_quality_mut() = MappingQuality::new(255);
        process(
            &record,
            &exon_starts,
            &exon_ends,
            &header,
            &params,
            &mut results,
        )
        .unwrap();
        assert_eq!(results.records.processed, 1);
        assert_eq!(results.records.ignored_flags, 1);
        assert_eq!(results.records.not_spliced, 0);
        assert_eq!(results.records.low_mapq, 0);
        assert_eq!(results.records.missing_mapq, 0);

        // Test partial novel junction
        let mut record = Record::default();
        let r3_name: ReadName = "partial1".parse().unwrap();
        *record.read_name_mut() = Some(r3_name);
        *record.flags_mut() = 0x0.into();
        *record.reference_sequence_id_mut() = Some(0);
        *record.alignment_start_mut() = Position::new(1);
        *record.cigar_mut() = "10M12N10M".parse().unwrap();
        *record.mapping_quality_mut() = MappingQuality::new(60);
        process(
            &record,
            &exon_starts,
            &exon_ends,
            &header,
            &params,
            &mut results,
        )
        .unwrap();
        assert_eq!(results.records.processed, 2);
        assert_eq!(results.records.ignored_flags, 1);
        assert_eq!(results.records.not_spliced, 0);
        assert_eq!(results.records.low_mapq, 0);
        assert_eq!(results.records.missing_mapq, 0);

        // Test partial novel junction (again for more read support)
        let mut record = Record::default();
        let r3_name: ReadName = "partial2".parse().unwrap();
        *record.read_name_mut() = Some(r3_name);
        *record.flags_mut() = 0x0.into();
        *record.reference_sequence_id_mut() = Some(0);
        *record.alignment_start_mut() = Position::new(1);
        *record.cigar_mut() = "10M12N10M".parse().unwrap();
        *record.mapping_quality_mut() = MappingQuality::new(60);
        process(
            &record,
            &exon_starts,
            &exon_ends,
            &header,
            &params,
            &mut results,
        )
        .unwrap();
        assert_eq!(results.records.processed, 3);
        assert_eq!(results.records.ignored_flags, 1);
        assert_eq!(results.records.not_spliced, 0);
        assert_eq!(results.records.low_mapq, 0);
        assert_eq!(results.records.missing_mapq, 0);

        // Test that supplementary alignments get counted
        let mut record = Record::default();
        let r4_name: ReadName = "supplementary_and_known2".parse().unwrap();
        *record.read_name_mut() = Some(r4_name);
        *record.flags_mut() = 0x800.into();
        *record.reference_sequence_id_mut() = Some(0);
        *record.alignment_start_mut() = Position::new(1);
        *record.cigar_mut() = "10M10N10M".parse().unwrap();
        *record.mapping_quality_mut() = MappingQuality::new(60);
        process(
            &record,
            &exon_starts,
            &exon_ends,
            &header,
            &params,
            &mut results,
        )
        .unwrap();
        assert_eq!(results.records.processed, 4);
        assert_eq!(results.records.ignored_flags, 1);
        assert_eq!(results.records.not_spliced, 0);
        assert_eq!(results.records.low_mapq, 0);
        assert_eq!(results.records.missing_mapq, 0);

        // Test that secondary alignments get ignored
        let mut record = Record::default();
        let r5_name: ReadName = "secondary".parse().unwrap();
        *record.read_name_mut() = Some(r5_name);
        *record.flags_mut() = 0x100.into();
        *record.reference_sequence_id_mut() = Some(0);
        *record.alignment_start_mut() = Position::new(1);
        *record.cigar_mut() = "10M10N10M".parse().unwrap();
        *record.mapping_quality_mut() = MappingQuality::new(60);
        process(
            &record,
            &exon_starts,
            &exon_ends,
            &header,
            &params,
            &mut results,
        )
        .unwrap();
        assert_eq!(results.records.processed, 4);
        assert_eq!(results.records.ignored_flags, 2);
        assert_eq!(results.records.not_spliced, 0);
        assert_eq!(results.records.low_mapq, 0);
        assert_eq!(results.records.missing_mapq, 0);

        // Test complete novel junction
        let mut record = Record::default();
        let r6_name: ReadName = "novel1".parse().unwrap();
        *record.read_name_mut() = Some(r6_name);
        *record.flags_mut() = 0x0.into();
        *record.reference_sequence_id_mut() = Some(0);
        *record.alignment_start_mut() = Position::new(1);
        *record.cigar_mut() = "8M15N8M".parse().unwrap();
        *record.mapping_quality_mut() = MappingQuality::new(60);
        process(
            &record,
            &exon_starts,
            &exon_ends,
            &header,
            &params,
            &mut results,
        )
        .unwrap();
        assert_eq!(results.records.processed, 5);
        assert_eq!(results.records.ignored_flags, 2);
        assert_eq!(results.records.not_spliced, 0);
        assert_eq!(results.records.low_mapq, 0);
        assert_eq!(results.records.missing_mapq, 0);

        // Test complete novel junction (again for more read support)
        let mut record = Record::default();
        let r6_name: ReadName = "novel2".parse().unwrap();
        *record.read_name_mut() = Some(r6_name);
        *record.flags_mut() = 0x0.into();
        *record.reference_sequence_id_mut() = Some(0);
        *record.alignment_start_mut() = Position::new(1);
        *record.cigar_mut() = "8M15N8M".parse().unwrap();
        *record.mapping_quality_mut() = MappingQuality::new(60);
        process(
            &record,
            &exon_starts,
            &exon_ends,
            &header,
            &params,
            &mut results,
        )
        .unwrap();
        assert_eq!(results.records.processed, 6);
        assert_eq!(results.records.ignored_flags, 2);
        assert_eq!(results.records.not_spliced, 0);
        assert_eq!(results.records.low_mapq, 0);
        assert_eq!(results.records.missing_mapq, 0);

        // Test fails MAPQ filter
        let mut record = Record::default();
        let r7_name: ReadName = "low_mapq".parse().unwrap();
        *record.read_name_mut() = Some(r7_name);
        *record.flags_mut() = 0x0.into();
        *record.reference_sequence_id_mut() = Some(0);
        *record.alignment_start_mut() = Position::new(1);
        *record.cigar_mut() = "10M10N10M".parse().unwrap();
        *record.mapping_quality_mut() = MappingQuality::new(20);
        process(
            &record,
            &exon_starts,
            &exon_ends,
            &header,
            &params,
            &mut results,
        )
        .unwrap();
        assert_eq!(results.records.processed, 6);
        assert_eq!(results.records.ignored_flags, 2);
        assert_eq!(results.records.not_spliced, 0);
        assert_eq!(results.records.low_mapq, 1);
        assert_eq!(results.records.missing_mapq, 0);

        // Test missing MAPQ
        let mut record = Record::default();
        let r8_name: ReadName = "missing_mapq".parse().unwrap();
        *record.read_name_mut() = Some(r8_name);
        *record.flags_mut() = 0x0.into();
        *record.reference_sequence_id_mut() = Some(0);
        *record.alignment_start_mut() = Position::new(1);
        *record.cigar_mut() = "10M10N10M".parse().unwrap();
        *record.mapping_quality_mut() = MappingQuality::new(255);
        process(
            &record,
            &exon_starts,
            &exon_ends,
            &header,
            &params,
            &mut results,
        )
        .unwrap();
        assert_eq!(results.records.processed, 6);
        assert_eq!(results.records.ignored_flags, 2);
        assert_eq!(results.records.not_spliced, 0);
        assert_eq!(results.records.low_mapq, 1);
        assert_eq!(results.records.missing_mapq, 1);

        // Test that intron is too short
        let mut record = Record::default();
        let r9_name: ReadName = "short".parse().unwrap();
        *record.read_name_mut() = Some(r9_name);
        *record.flags_mut() = 0x0.into();
        *record.reference_sequence_id_mut() = Some(0);
        *record.alignment_start_mut() = Position::new(1);
        *record.cigar_mut() = "5M5N5M".parse().unwrap();
        *record.mapping_quality_mut() = MappingQuality::new(60);
        process(
            &record,
            &exon_starts,
            &exon_ends,
            &header,
            &params,
            &mut results,
        )
        .unwrap();
        assert_eq!(results.records.processed, 7); // Still gets processed, will be filtered later
        assert_eq!(results.records.ignored_flags, 2);
        assert_eq!(results.records.not_spliced, 0);
        assert_eq!(results.records.low_mapq, 1);
        assert_eq!(results.records.missing_mapq, 1);

        // That that reads not spliced are ignored
        let mut record = Record::default();
        let r10_name: ReadName = "not_spliced".parse().unwrap();
        *record.read_name_mut() = Some(r10_name);
        *record.flags_mut() = 0x0.into();
        *record.reference_sequence_id_mut() = Some(0);
        *record.alignment_start_mut() = Position::new(1);
        *record.cigar_mut() = "10M".parse().unwrap();
        *record.mapping_quality_mut() = MappingQuality::new(60);
        process(
            &record,
            &exon_starts,
            &exon_ends,
            &header,
            &params,
            &mut results,
        )
        .unwrap();
        assert_eq!(results.records.processed, 7);
        assert_eq!(results.records.ignored_flags, 2);
        assert_eq!(results.records.not_spliced, 1);
        assert_eq!(results.records.low_mapq, 1);
        assert_eq!(results.records.missing_mapq, 1);

        // Test unannoted reference
        let mut record = Record::default();
        let r11_name: ReadName = "unannotated1".parse().unwrap();
        *record.read_name_mut() = Some(r11_name);
        *record.flags_mut() = 0x0.into();
        *record.reference_sequence_id_mut() = Some(1);
        *record.alignment_start_mut() = Position::new(1);
        *record.cigar_mut() = "10M10N10M".parse().unwrap();
        *record.mapping_quality_mut() = MappingQuality::new(60);
        process(
            &record,
            &exon_starts,
            &exon_ends,
            &header,
            &params,
            &mut results,
        )
        .unwrap();
        assert_eq!(results.records.processed, 8);
        assert_eq!(results.records.ignored_flags, 2);
        assert_eq!(results.records.not_spliced, 1);
        assert_eq!(results.records.low_mapq, 1);
        assert_eq!(results.records.missing_mapq, 1);

        // Test unannoted reference (again for more read support)
        let mut record = Record::default();
        let r11_name: ReadName = "unannotated2".parse().unwrap();
        *record.read_name_mut() = Some(r11_name);
        *record.flags_mut() = 0x0.into();
        *record.reference_sequence_id_mut() = Some(1);
        *record.alignment_start_mut() = Position::new(1);
        *record.cigar_mut() = "10M10N10M".parse().unwrap();
        *record.mapping_quality_mut() = MappingQuality::new(60);
        process(
            &record,
            &exon_starts,
            &exon_ends,
            &header,
            &params,
            &mut results,
        )
        .unwrap();
        assert_eq!(results.records.processed, 9);
        assert_eq!(results.records.ignored_flags, 2);
        assert_eq!(results.records.not_spliced, 1);
        assert_eq!(results.records.low_mapq, 1);
        assert_eq!(results.records.missing_mapq, 1);

        // Test summarize
        summarize(&mut results, &params);

        assert_eq!(results.summary.total_rejected_junctions, 1);
        assert_eq!(results.summary.intron_too_short, 1);
        assert_eq!(results.summary.junctions_with_not_enough_read_support, 1);
        assert_eq!(results.summary.known_junctions, 1);
        assert_eq!(results.summary.known_junctions_read_support, 2);
        assert_eq!(results.summary.partial_novel_junctions, 1);
        assert_eq!(results.summary.partial_novel_junctions_read_support, 2);
        assert_eq!(results.summary.complete_novel_junctions, 1);
        assert_eq!(results.summary.complete_novel_junctions_read_support, 2);
        assert_eq!(results.summary.unannotated_reference_junctions, 1);
        assert_eq!(
            results.summary.unannotated_reference_junctions_read_support,
            2
        );
        assert_eq!(results.summary.total_junctions, 4);
        assert_eq!(results.summary.total_junctions_read_support, 8);
        assert_eq!(results.summary.known_junctions_percent, 33.33333333333333);
        assert_eq!(
            results.summary.partial_novel_junctions_percent,
            33.33333333333333
        );
        assert_eq!(
            results.summary.complete_novel_junctions_percent,
            33.33333333333333
        );
        assert_eq!(results.summary.average_junction_read_support, 2.0);
        assert_eq!(results.summary.average_known_junction_read_support, 2.0);
        assert_eq!(
            results.summary.average_partial_novel_junction_read_support,
            2.0
        );
        assert_eq!(
            results.summary.average_complete_novel_junction_read_support,
            2.0
        );
    }
}
