//! Results related to the `ngs derive junction_annotation` subcommand.

use serde::ser::SerializeStruct;
use serde::Serialize;
use serde::Serializer;
use std::collections::HashMap;
use std::num::NonZeroUsize;

/// Lists of annotated junctions.
#[derive(Clone, Default)]
pub struct JunctionAnnotations {
    /// Known junctions. The outer key is the referece name, and the value is another
    /// HashMap. The inner key is the (start, end) coordinates of the junction,
    /// and the value is the number of spliced reads that support the junction.
    pub known: HashMap<String, HashMap<(NonZeroUsize, NonZeroUsize), usize>>,

    /// Partially novel junctions. The outer key is the referece name, and the value is another
    /// HashMap. The inner key is the (start, end) coordinates of the junction,
    /// and the value is the number of spliced reads that support the junction.
    pub partial_novel: HashMap<String, HashMap<(NonZeroUsize, NonZeroUsize), usize>>,

    /// Complete novel junctions. The outer key is the referece name, and the value is another
    /// HashMap. The inner key is the (start, end) coordinates of the junction,
    /// and the value is the number of spliced reads that support the junction.
    pub complete_novel: HashMap<String, HashMap<(NonZeroUsize, NonZeroUsize), usize>>,

    /// Junctions on reference sequences for which junction annotations were not found.
    /// The outer key is the referece name, and the value is another
    /// HashMap. The inner key is the (start, end) coordinates of the junction,
    /// and the value is the number of spliced reads that support the junction.
    pub unannotated_reference: HashMap<String, HashMap<(NonZeroUsize, NonZeroUsize), usize>>,
}

impl Serialize for JunctionAnnotations {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let mut known = Vec::new();
        for (ref_name, junctions) in &self.known {
            for ((start, end), count) in junctions {
                known.push((ref_name, start, end, count));
            }
        }

        let mut partial_novel = Vec::new();
        for (ref_name, junctions) in &self.partial_novel {
            for ((start, end), count) in junctions {
                partial_novel.push((ref_name, start, end, count));
            }
        }

        let mut complete_novel = Vec::new();
        for (ref_name, junctions) in &self.complete_novel {
            for ((start, end), count) in junctions {
                complete_novel.push((ref_name, start, end, count));
            }
        }

        let mut unannotated_reference = Vec::new();
        for (ref_name, junctions) in &self.unannotated_reference {
            for ((start, end), count) in junctions {
                unannotated_reference.push((ref_name, start, end, count));
            }
        }

        let mut s = serializer.serialize_struct("JunctionAnnotations", 4)?;
        s.serialize_field("known", &known)?;
        s.serialize_field("partial_novel", &partial_novel)?;
        s.serialize_field("complete_novel", &complete_novel)?;
        s.serialize_field("unannotated_reference", &unannotated_reference)?;
        s.end()
    }
}

/// General record metrics that are tallied as a part of the
/// junction_annotation subcommand.
#[derive(Clone, Default, Serialize)]
pub struct RecordMetrics {
    /// The number of records that have been fully processed.
    pub processed: usize,

    /// The number of records that have been ignored because of their flags.
    /// (i.e. they were unmapped, duplicates, secondary, or supplementary)
    /// The last 3 conditions can be toggled on/off with CL flags
    pub ignored_flags: usize,

    /// The number of records that have been ignored because they were not
    /// spliced.
    pub not_spliced: usize,

    /// The number of records with junctions that have been ignored because
    /// they failed the MAPQ filter.
    pub low_mapq: usize,

    /// The number of records whose MAPQ couldn't be parsed and were thus ignored.
    pub missing_mapq: usize,
}

/// Summary statistics for the junction_annotation subcommand.
#[derive(Clone, Default, Serialize)]
pub struct SummaryResults {
    /// The total number of junctions observed in the file.
    pub total_junctions: usize,

    /// The total number of spliced reads observed in the file.
    pub total_spliced_reads: usize,

    /// The total number of known junctions observed in the file.
    pub known_junctions: usize,

    ///The total number of partially novel junctions observed in the file.
    pub partial_novel_junctions: usize,

    /// The total number of complete novel junctions observed in the file.
    pub complete_novel_junctions: usize,

    /// The total number of junctions on reference sequences for which junction
    /// annotations were not found.
    pub unannotated_reference_junctions: usize,

    /// The total number of known spliced reads observed in the file.
    pub known_spliced_reads: usize,

    /// The total number of partially novel spliced reads observed in the file.
    pub partial_novel_spliced_reads: usize,

    /// The total number of complete novel spliced reads observed in the file.
    pub complete_novel_spliced_reads: usize,

    /// The total number of spliced reads on reference sequences for which
    /// junction annotations were not found.
    pub unannotated_reference_spliced_reads: usize,

    /// The total number of junctions which were discarded due to lack of
    /// read support.
    pub junctions_with_not_enough_read_support: usize,

    /// The number of junctions that have been ignored because
    /// they failed the min_intron_length filter.
    pub intron_too_short: usize,
}

/// Main Results struct. This struct aggregates all of the minor metrics structs
/// outlined in this file so they can be kept track of as a unit.
#[derive(Clone, Default, Serialize)]
pub struct JunctionAnnotationResults {
    /// Lists of annotated junctions.
    pub junction_annotations: JunctionAnnotations,

    /// General record metrics.
    pub records: RecordMetrics,

    /// Summary statistics for the junction_annotation subcommand.
    pub summary: SummaryResults,
}
