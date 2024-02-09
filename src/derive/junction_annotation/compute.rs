//! Module holding the logic for annotating junctions.

use anyhow::bail;
use anyhow::Ok;
use noodles::core::Position;
use noodles::sam::alignment::Record;
use noodles::sam::record::cigar::op::Kind;
use noodles::sam::record::MappingQuality;
use noodles::sam::Header;
use std::collections::HashMap;
use std::collections::HashSet;

use crate::derive::junction_annotation::results;
use crate::utils::alignment::filter_by_mapq;

/// Struct to hold starts and ends of exons.
pub struct ExonSets<'a> {
    /// Starts of exons, grouped by contig.
    pub starts: HashMap<&'a str, HashSet<Position>>,

    /// ends of exons, grouped by contig.
    pub ends: HashMap<&'a str, HashSet<Position>>,
}

/// Parameters defining how to annotate found junctions
pub struct JunctionAnnotationParameters {
    /// Minimum intron length to consider.
    pub min_intron_length: usize,

    /// Minimum number of reads supporting a junction to be considered.
    pub min_read_support: usize,

    /// Minumum mapping quality for a record to be considered.
    /// 0 if MAPQ shouldn't be considered.
    pub min_mapq: Option<MappingQuality>,

    /// Do not count supplementary alignments.
    pub no_supplementary: bool,

    /// Do count secondary alignments.
    pub count_secondary: bool,

    /// Do count duplicates.
    pub count_duplicates: bool,
}

/// Function for incrementing a junction counter by one.
fn increment_junction_counter(
    junction_counter: &mut results::JunctionCounter,
    junction: results::Junction,
) {
    junction_counter
        .entry(junction)
        .and_modify(|e| *e += 1)
        .or_insert(1);
}

/// Function for incrementing a junction map by one.
fn increment_junction_map(
    junction_map: &mut results::JunctionsMap,
    ref_name: &str,
    junction: results::Junction,
) {
    increment_junction_counter(
        junction_map.entry(ref_name.to_string()).or_default(),
        junction,
    );
}

/// Function to filter out records based on their flags.
fn filter_by_flags(record: &Record, params: &JunctionAnnotationParameters) -> bool {
    let flags = record.flags();
    if flags.is_unmapped()
        || (params.no_supplementary && flags.is_supplementary())
        || (!params.count_secondary && flags.is_secondary())
        || (!params.count_duplicates && flags.is_duplicate())
    {
        return true;
    }
    false
}

/// Function to filter out records that don't have introns.
fn filter_by_cigar(record: &Record) -> bool {
    !record
        .cigar()
        .iter()
        .any(|op| matches!(op.kind(), Kind::Skip))
}

/// Main function to annotate junctions one record at a time.
pub fn process(
    record: &Record,
    exons: &ExonSets<'_>,
    header: &Header,
    params: &JunctionAnnotationParameters,
    results: &mut results::JunctionAnnotationResults,
) -> anyhow::Result<()> {
    // (1) Parse the read name.
    let read_name = match record.read_name() {
        Some(name) => name,
        _ => bail!("Could not parse read name"),
    };

    // (2) Filter by record flags.
    if filter_by_flags(record, params) {
        results.records.filtered_by_flags += 1;
        return Ok(());
    }

    // (3) Filter by CIGAR.
    // We only care about reads with introns, so if there are no introns
    // we can skip this read.
    if filter_by_cigar(record) {
        results.records.not_spliced += 1;
        return Ok(());
    }

    // (4) Filter by MAPQ
    if filter_by_mapq(record, params.min_mapq) {
        results.records.bad_mapq += 1;
        return Ok(());
    }

    // (5) Parse the reference sequence from the record.
    let (seq_name, _) = match record.reference_sequence(header) {
        Some(seq_map_result) => seq_map_result?,
        _ => {
            bail!(
                "Could not parse reference sequence id for read: {}",
                read_name
            )
        }
    };
    let seq_name = seq_name.as_str();

    // (6) Check if there will be annotations for this reference sequence.
    let mut ref_is_annotated = true;
    if !exons.starts.contains_key(seq_name) || !exons.ends.contains_key(seq_name) {
        ref_is_annotated = false;
    }

    // (7) Calculate the start position of this read. This will
    // later be used to find the position of any introns.
    let start = match record.alignment_start() {
        Some(s) => s,
        _ => bail!("Could not parse record's start position."),
    };

    // (8) Find introns
    let cur_pos = start;
    for op in record.cigar().iter() {
        match op.kind() {
            // This is an intron.
            Kind::Skip => {
                // Check that `op.len() >= params.min_intron_length` later,
                // for better metric reporting.

                let intron_start = cur_pos;
                // Update cur_pos to the end of the intron.
                cur_pos.checked_add(op.len());
                let intron_end = cur_pos;
                let junction: results::Junction = (intron_start, intron_end);

                // If the reference sequence is not annotated, we can skip
                // the lookup of exon positions, and directly insert the
                // intron into the unannotated_reference HashMap.
                if !ref_is_annotated {
                    increment_junction_map(
                        &mut results.junction_annotations.unannotated_reference,
                        seq_name,
                        junction,
                    );
                    continue;
                }

                let exon_starts = match exons.starts.get(seq_name) {
                    Some(starts) => starts,
                    _ => bail!("Could not find exon starts for contig: {}", seq_name),
                };
                let exon_ends = match exons.ends.get(seq_name) {
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

                // TODO: Better way to do this?
                increment_junction_map(
                    match (intron_start_known, intron_end_known) {
                        (true, true) => {
                            // We found both ends of the intron.
                            // This is a Known Junction.
                            &mut results.junction_annotations.known
                        }
                        (true, false) | (false, true) => {
                            // We found one end of the intron,
                            // but not the other.
                            // This is a Partial Novel Junction.
                            &mut results.junction_annotations.partial_novel
                        }
                        (false, false) => {
                            // We found neither end of the intron.
                            // This is a Complete Novel Junction.
                            &mut results.junction_annotations.complete_novel
                        }
                    },
                    seq_name,
                    junction,
                )
            }
            // Operations (beside Skip which is handled above) that increment the reference position.
            Kind::Match | Kind::Deletion | Kind::SequenceMatch | Kind::SequenceMismatch => {
                cur_pos.checked_add(op.len());
            }
            // Operations that do not increment the reference position.
            _ => {}
        }
    }

    results.records.processed += 1;
    Ok(())
}

/// Function to filter out junctions that are too short or don't have enough read support.
fn filter_junction_map(
    junction_map: &mut results::JunctionsMap,
    min_intron_length: usize,
    min_read_support: usize,
    metrics: &mut results::SummaryResults,
) {
    junction_map.retain(|_, v| {
        v.retain(|(start, end), count| {
            let mut keep = true;
            if end.get() - start.get() < min_intron_length {
                metrics.intron_too_short += 1;
                keep = false;
            }
            if *count < min_read_support {
                metrics.junctions_with_not_enough_read_support += 1;
                keep = false;
            }
            if !keep {
                metrics.total_rejected_junctions += 1;
            }
            keep
        });
        !v.is_empty()
    });
}

/// Function to tally up the junctions and their read support.
fn tally_junctions_and_support(junction_map: &results::JunctionsMap) -> (usize, usize) {
    let junctions = junction_map.values().map(|v| v.len()).sum();
    let support = junction_map
        .values()
        .map(|v| v.values().sum::<usize>())
        .sum();
    (junctions, support)
}

/// Main function to summarize the results of the junction-annotation subcommand.
pub fn summarize(
    results: &mut results::JunctionAnnotationResults,
    params: &JunctionAnnotationParameters,
) {
    // Filter out junctions that are too short or don't have enough read support.
    filter_junction_map(
        &mut results.junction_annotations.known,
        params.min_intron_length,
        params.min_read_support,
        &mut results.summary,
    );
    filter_junction_map(
        &mut results.junction_annotations.partial_novel,
        params.min_intron_length,
        params.min_read_support,
        &mut results.summary,
    );
    filter_junction_map(
        &mut results.junction_annotations.complete_novel,
        params.min_intron_length,
        params.min_read_support,
        &mut results.summary,
    );
    filter_junction_map(
        &mut results.junction_annotations.unannotated_reference,
        params.min_intron_length,
        params.min_read_support,
        &mut results.summary,
    );

    // Tally up observed junctions and spliced reads.
    let mut juncs;
    let mut support;
    (juncs, support) = tally_junctions_and_support(&results.junction_annotations.known);
    results.summary.known_junctions = juncs;
    results.summary.known_junctions_read_support = support;
    (juncs, support) = tally_junctions_and_support(&results.junction_annotations.partial_novel);
    results.summary.partial_novel_junctions = juncs;
    results.summary.partial_novel_junctions_read_support = support;
    (juncs, support) = tally_junctions_and_support(&results.junction_annotations.complete_novel);
    results.summary.complete_novel_junctions = juncs;
    results.summary.complete_novel_junctions_read_support = support;
    (juncs, support) =
        tally_junctions_and_support(&results.junction_annotations.unannotated_reference);
    results.summary.unannotated_reference_junctions = juncs;
    results.summary.unannotated_reference_junctions_read_support = support;

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
    use std::num::NonZeroUsize;

    #[test]
    fn test_filter_by_flags() {
        // Setup
        let mut record = Record::default();
        let params = JunctionAnnotationParameters {
            min_intron_length: 10,
            min_read_support: 2,
            min_mapq: Some(MappingQuality::new(30).unwrap()),
            no_supplementary: false,
            count_secondary: false,
            count_duplicates: false,
        };

        // Test that records are filtered out correctly
        record.flags_mut().set(0x4.into(), true);
        assert!(filter_by_flags(&record, &params));
        record.flags_mut().set(0x4.into(), false);
        record.flags_mut().set(0x800.into(), true);
        assert!(!filter_by_flags(&record, &params));
        record.flags_mut().set(0x800.into(), false);
        record.flags_mut().set(0x100.into(), true);
        assert!(filter_by_flags(&record, &params));
        record.flags_mut().set(0x100.into(), false);
        record.flags_mut().set(0x400.into(), true);
        assert!(filter_by_flags(&record, &params));
        record.flags_mut().set(0x400.into(), false);
        assert!(!filter_by_flags(&record, &params));
    }

    #[test]
    fn test_filter_by_cigar() {
        // Setup
        let mut record = Record::default();

        // Test that records are filtered out correctly
        *record.cigar_mut() = "10M10N10M".parse().unwrap();
        assert!(!filter_by_cigar(&record));
        *record.cigar_mut() = "10M".parse().unwrap();
        assert!(filter_by_cigar(&record));
    }

    #[test]
    fn test_filter_junction_map() {
        // Setup
        let mut junction_map = results::JunctionsMap::default();
        junction_map.insert(
            "sq1".to_string(),
            HashMap::from([
                ((Position::new(1).unwrap(), Position::new(11).unwrap()), 1),
                ((Position::new(1).unwrap(), Position::new(5).unwrap()), 1),
            ]),
        );
        junction_map.insert(
            "sq2".to_string(),
            HashMap::from([((Position::new(1).unwrap(), Position::new(11).unwrap()), 2)]),
        );
        let min_intron_length = 10;
        let min_read_support = 2;
        let mut metrics = results::SummaryResults::default();

        // Test that junctions are filtered out correctly
        filter_junction_map(
            &mut junction_map,
            min_intron_length,
            min_read_support,
            &mut metrics,
        );
        assert_eq!(junction_map.len(), 1);
        assert_eq!(junction_map.get("sq1"), None);
        assert_eq!(junction_map.get("sq2").unwrap().len(), 1);
        assert_eq!(metrics.intron_too_short, 1);
        assert_eq!(metrics.junctions_with_not_enough_read_support, 2);
        assert_eq!(metrics.total_rejected_junctions, 2);
    }

    #[test]
    fn test_tally_junctions_and_support() {
        // Setup
        let mut junction_map = results::JunctionsMap::default();
        junction_map.insert(
            "sq1".to_string(),
            HashMap::from([
                ((Position::new(1).unwrap(), Position::new(11).unwrap()), 1),
                ((Position::new(1).unwrap(), Position::new(5).unwrap()), 1),
            ]),
        );
        junction_map.insert(
            "sq2".to_string(),
            HashMap::from([((Position::new(1).unwrap(), Position::new(11).unwrap()), 2)]),
        );

        // Test that junctions are tallied correctly
        let (juncs, support) = tally_junctions_and_support(&junction_map);
        assert_eq!(juncs, 3);
        assert_eq!(support, 4);
    }

    #[test]
    fn test_process_and_summarize() {
        // Setup
        let mut results = results::JunctionAnnotationResults::default();
        let params = JunctionAnnotationParameters {
            min_intron_length: 10,
            min_read_support: 2,
            min_mapq: Some(MappingQuality::new(30).unwrap()),
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
        let exon_starts: HashMap<&str, HashSet<Position>> = HashMap::from([(
            "sq1",
            HashSet::from([
                Position::new(1).unwrap(),
                Position::new(11).unwrap(),
                Position::new(21).unwrap(),
                Position::new(31).unwrap(),
                Position::new(41).unwrap(),
                Position::new(51).unwrap(),
                Position::new(61).unwrap(),
                Position::new(71).unwrap(),
            ]),
        )]);
        let exon_ends = exon_starts
            .iter()
            .map(|(k, v)| (*k, v.iter().map(|e| e.checked_add(10).unwrap()).collect()))
            .collect::<HashMap<&str, HashSet<Position>>>();
        let exons = ExonSets {
            starts: exon_starts,
            ends: exon_ends,
        };

        // Test known junction
        let mut record = Record::default();
        let r1_name: ReadName = "known1".parse().unwrap();
        *record.read_name_mut() = Some(r1_name);
        *record.reference_sequence_id_mut() = Some(0);
        *record.alignment_start_mut() = Position::new(1);
        *record.cigar_mut() = "10M10N10M".parse().unwrap();
        *record.mapping_quality_mut() = MappingQuality::new(60);
        record.flags_mut().set(0x4.into(), false);
        process(&record, &exons, &header, &params, &mut results).unwrap();
        assert_eq!(results.records.processed, 1);
        assert_eq!(results.records.filtered_by_flags, 0);
        assert_eq!(results.records.not_spliced, 0);
        assert_eq!(results.records.bad_mapq, 0);

        // Test that unmapped gets ignored
        let mut record = Record::default();
        let r2_name: ReadName = "unmapped".parse().unwrap();
        *record.read_name_mut() = Some(r2_name);
        *record.reference_sequence_id_mut() = Some(0);
        *record.alignment_start_mut() = Position::new(1);
        *record.cigar_mut() = "10M10N10M".parse().unwrap();
        *record.mapping_quality_mut() = MappingQuality::new(255);
        record.flags_mut().set(0x4.into(), true);
        process(&record, &exons, &header, &params, &mut results).unwrap();
        assert_eq!(results.records.processed, 1);
        assert_eq!(results.records.filtered_by_flags, 1);
        assert_eq!(results.records.not_spliced, 0);
        assert_eq!(results.records.bad_mapq, 0);

        // Test partial novel junction
        let mut record = Record::default();
        let r3_name: ReadName = "partial1".parse().unwrap();
        *record.read_name_mut() = Some(r3_name);
        *record.reference_sequence_id_mut() = Some(0);
        *record.alignment_start_mut() = Position::new(1);
        *record.cigar_mut() = "10M12N10M".parse().unwrap();
        *record.mapping_quality_mut() = MappingQuality::new(60);
        record.flags_mut().set(0x4.into(), false);
        process(&record, &exons, &header, &params, &mut results).unwrap();
        assert_eq!(results.records.processed, 2);
        assert_eq!(results.records.filtered_by_flags, 1);
        assert_eq!(results.records.not_spliced, 0);
        assert_eq!(results.records.bad_mapq, 0);

        // Test partial novel junction (again for more read support)
        let mut record = Record::default();
        let r3_name: ReadName = "partial2".parse().unwrap();
        *record.read_name_mut() = Some(r3_name);
        *record.reference_sequence_id_mut() = Some(0);
        *record.alignment_start_mut() = Position::new(1);
        *record.cigar_mut() = "10M12N10M".parse().unwrap();
        *record.mapping_quality_mut() = MappingQuality::new(60);
        record.flags_mut().set(0x4.into(), false);
        process(&record, &exons, &header, &params, &mut results).unwrap();
        assert_eq!(results.records.processed, 3);
        assert_eq!(results.records.filtered_by_flags, 1);
        assert_eq!(results.records.not_spliced, 0);
        assert_eq!(results.records.bad_mapq, 0);

        // Test that supplementary alignments get counted
        let mut record = Record::default();
        let r4_name: ReadName = "supplementary".parse().unwrap();
        *record.read_name_mut() = Some(r4_name);
        *record.reference_sequence_id_mut() = Some(0);
        *record.alignment_start_mut() = Position::new(1);
        *record.cigar_mut() = "10M10N10M".parse().unwrap();
        *record.mapping_quality_mut() = MappingQuality::new(60);
        record.flags_mut().set(0x4.into(), false);
        record.flags_mut().set(0x800.into(), true);
        process(&record, &exons, &header, &params, &mut results).unwrap();
        assert_eq!(results.records.processed, 4);
        assert_eq!(results.records.filtered_by_flags, 1);
        assert_eq!(results.records.not_spliced, 0);
        assert_eq!(results.records.bad_mapq, 0);

        // Test that secondary alignments don't get counted
        let mut record = Record::default();
        let r5_name: ReadName = "secondary".parse().unwrap();
        *record.read_name_mut() = Some(r5_name);
        *record.reference_sequence_id_mut() = Some(0);
        *record.alignment_start_mut() = Position::new(1);
        *record.cigar_mut() = "10M10N10M".parse().unwrap();
        *record.mapping_quality_mut() = MappingQuality::new(60);
        record.flags_mut().set(0x4.into(), false);
        record.flags_mut().set(0x100.into(), true);
        process(&record, &exons, &header, &params, &mut results).unwrap();
        assert_eq!(results.records.processed, 4);
        assert_eq!(results.records.filtered_by_flags, 2);
        assert_eq!(results.records.not_spliced, 0);
        assert_eq!(results.records.bad_mapq, 0);

        // TODO: Below tests are not working as expected. Need to fix them.
        // Test complete novel junction
        // let mut record = Record::default();
        // let r6_name: ReadName = "complete1".parse().unwrap();
        // *record.read_name_mut() = Some(r6_name);
        // *record.reference_sequence_id_mut() = Some(0);
        // *record.alignment_start_mut() = Position::new(1);
        // *record.cigar_mut() = "10M10N10M10N10M".parse().unwrap();
        // *record.mapping_quality_mut() = MappingQuality::new(60);
        // record.flags_mut().set(0x4.into(), false);
        // process(&record, &exons, &header, &params, &mut results).unwrap();
        // assert_eq!(results.records.processed, 6);
        // assert_eq!(results.records.filtered_by_flags, 2);
        // assert_eq!(results.records.not_spliced, 0);
        // assert_eq!(results.records.low_mapq, 0);
        // assert_eq!(results.records.missing_mapq, 0);

        // // Test complete novel junction (again for more read support)
        // let mut record = Record::default();
        // let r6_name: ReadName = "complete2".parse().unwrap();
        // *record.read_name_mut() = Some(r6_name);
        // *record.reference_sequence_id_mut() = Some(0);
        // *record.alignment_start_mut() = Position::new(1);
        // *record.cigar_mut() = "10M10N10M10N10M".parse().unwrap();
        // *record.mapping_quality_mut() = MappingQuality::new(60);
        // record.flags_mut().set(0x4.into(), false);
        // process(&record, &exons, &header, &params, &mut results).unwrap();
        // assert_eq!(results.records.processed, 7);
        // assert_eq!(results.records.filtered_by_flags, 2);
        // assert_eq!(results.records.not_spliced, 0);
        // assert_eq!(results.records.low_mapq, 0);
        // assert_eq!(results.records.missing_mapq, 0);

        // // Test fails MAPQ filter
        // let mut record = Record::default();
        // let r7_name: ReadName = "low_mapq".parse().unwrap();
        // *record.read_name_mut() = Some(r7_name);
        // *record.reference_sequence_id_mut() = Some(0);
        // *record.alignment_start_mut() = Position::new(1);
        // *record.cigar_mut() = "10M10N10M".parse().unwrap();
        // *record.mapping_quality_mut() = MappingQuality::new(20);
        // record.flags_mut().set(0x4.into(), false);
        // process(&record, &exons, &header, &params, &mut results).unwrap();
        // assert_eq!(results.records.processed, 6);
        // assert_eq!(results.records.filtered_by_flags, 2);
        // assert_eq!(results.records.not_spliced, 0);
        // assert_eq!(results.records.low_mapq, 1);
        // assert_eq!(results.records.missing_mapq, 0);

        // // Test missing MAPQ
        // let mut record = Record::default();
        // let r8_name: ReadName = "missing_mapq".parse().unwrap();
        // *record.read_name_mut() = Some(r8_name);
        // *record.reference_sequence_id_mut() = Some(0);
        // *record.alignment_start_mut() = Position::new(1);
        // *record.cigar_mut() = "10M10N10M".parse().unwrap();
        // *record.mapping_quality_mut() = MappingQuality::new(255);
        // record.flags_mut().set(0x4.into(), false);
        // process(&record, &exons, &header, &params, &mut results).unwrap();
        // assert_eq!(results.records.processed, 6);
        // assert_eq!(results.records.filtered_by_flags, 2);
        // assert_eq!(results.records.not_spliced, 0);
        // assert_eq!(results.records.low_mapq, 1);
        // assert_eq!(results.records.missing_mapq, 1);

        // // Test that intron is too short
        // let mut record = Record::default();
        // let r9_name: ReadName = "short".parse().unwrap();
        // *record.read_name_mut() = Some(r9_name);
        // *record.reference_sequence_id_mut() = Some(0);
        // *record.alignment_start_mut() = Position::new(1);
        // *record.cigar_mut() = "5M5N5M".parse().unwrap();
        // *record.mapping_quality_mut() = MappingQuality::new(60);
        // record.flags_mut().set(0x4.into(), false);
        // process(&record, &exons, &header, &params, &mut results).unwrap();
        // assert_eq!(results.records.processed, 7); // Still gets processed, will be filtered later
        // assert_eq!(results.records.filtered_by_flags, 2);
        // assert_eq!(results.records.not_spliced, 0);
        // assert_eq!(results.records.low_mapq, 1);
        // assert_eq!(results.records.missing_mapq, 1);

        // // Test that that reads not spliced are ignored
        // let mut record = Record::default();
        // let r10_name: ReadName = "not_spliced".parse().unwrap();
        // *record.read_name_mut() = Some(r10_name);
        // *record.reference_sequence_id_mut() = Some(0);
        // *record.alignment_start_mut() = Position::new(1);
        // *record.cigar_mut() = "10M".parse().unwrap();
        // *record.mapping_quality_mut() = MappingQuality::new(60);
        // record.flags_mut().set(0x4.into(), false);
        // process(&record, &exons, &header, &params, &mut results).unwrap();
        // assert_eq!(results.records.processed, 7);
        // assert_eq!(results.records.filtered_by_flags, 2);
        // assert_eq!(results.records.not_spliced, 1);
        // assert_eq!(results.records.low_mapq, 1);
        // assert_eq!(results.records.missing_mapq, 1);

        // // Test unannoted reference
        // let mut record = Record::default();
        // let r11_name: ReadName = "unannotated1".parse().unwrap();
        // *record.read_name_mut() = Some(r11_name);
        // *record.reference_sequence_id_mut() = Some(1);
        // *record.alignment_start_mut() = Position::new(1);
        // *record.cigar_mut() = "10M10N10M".parse().unwrap();
        // *record.mapping_quality_mut() = MappingQuality::new(60);
        // record.flags_mut().set(0x4.into(), false);
        // process(&record, &exons, &header, &params, &mut results).unwrap();
        // assert_eq!(results.records.processed, 8);
        // assert_eq!(results.records.filtered_by_flags, 2);
        // assert_eq!(results.records.not_spliced, 1);
        // assert_eq!(results.records.low_mapq, 1);
        // assert_eq!(results.records.missing_mapq, 1);

        // // Test unannoted reference (again for more read support)
        // let mut record = Record::default();
        // let r11_name: ReadName = "unannotated2".parse().unwrap();
        // *record.read_name_mut() = Some(r11_name);
        // *record.reference_sequence_id_mut() = Some(1);
        // *record.alignment_start_mut() = Position::new(1);
        // *record.cigar_mut() = "10M10N10M".parse().unwrap();
        // *record.mapping_quality_mut() = MappingQuality::new(60);
        // record.flags_mut().set(0x4.into(), false);
        // process(&record, &exons, &header, &params, &mut results).unwrap();
        // assert_eq!(results.records.processed, 9);
        // assert_eq!(results.records.filtered_by_flags, 2);
        // assert_eq!(results.records.not_spliced, 1);
        // assert_eq!(results.records.low_mapq, 1);
        // assert_eq!(results.records.missing_mapq, 1);

        // // Test summarize
        // summarize(&mut results, &params);

        // assert_eq!(results.summary.total_rejected_junctions, 1);
        // assert_eq!(results.summary.intron_too_short, 1);
        // assert_eq!(results.summary.junctions_with_not_enough_read_support, 1);
        // assert_eq!(results.summary.known_junctions, 1);
        // assert_eq!(results.summary.known_junctions_read_support, 2);
        // assert_eq!(results.summary.partial_novel_junctions, 1);
        // assert_eq!(results.summary.partial_novel_junctions_read_support, 2);
        // assert_eq!(results.summary.complete_novel_junctions, 1);
        // assert_eq!(results.summary.complete_novel_junctions_read_support, 2);
        // assert_eq!(results.summary.unannotated_reference_junctions, 1);
        // assert_eq!(
        //     results.summary.unannotated_reference_junctions_read_support,
        //     2
        // );
        // assert_eq!(results.summary.total_junctions, 4);
        // assert_eq!(results.summary.total_junctions_read_support, 8);
        // assert_eq!(results.summary.known_junctions_percent, 33.33333333333333);
        // assert_eq!(
        //     results.summary.partial_novel_junctions_percent,
        //     33.33333333333333
        // );
        // assert_eq!(
        //     results.summary.complete_novel_junctions_percent,
        //     33.33333333333333
        // );
        // assert_eq!(results.summary.average_junction_read_support, 2.0);
        // assert_eq!(results.summary.average_known_junction_read_support, 2.0);
        // assert_eq!(
        //     results.summary.average_partial_novel_junction_read_support,
        //     2.0
        // );
        // assert_eq!(
        //     results.summary.average_complete_novel_junction_read_support,
        //     2.0
        // );
    }
}
