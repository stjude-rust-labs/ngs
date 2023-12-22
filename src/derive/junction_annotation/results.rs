//! Results related to the `ngs derive junction_annotation` subcommand.

use serde::Deserialize;
use serde::Serialize;
use std::collections::HashMap;
use std::num::NonZeroUsize;

/// Lists of annotated junctions.
#[derive(Clone, Default, Serialize, Deserialize)]
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

/// General record metrics that are tallied as a part of the
/// junction_annotation subcommand.
#[derive(Clone, Default, Serialize, Deserialize)]
pub struct RecordMetrics {
    /// The number of records that have been fully processed.
    /// Should equal Metrics::SummaryMetrics::total_spliced_reads.
    pub processed: usize,

    /// The number of records that couldn't be parsed.
    pub couldnt_parse: usize,

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
}

/// Summary statistics for the junction_annotation subcommand.
#[derive(Clone, Default, Serialize, Deserialize)]
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
#[derive(Clone, Default, Serialize, Deserialize)]
pub struct JunctionAnnotationResults {
    /// Lists of annotated junctions.
    pub junction_annotations: JunctionAnnotations,

    /// General record metrics.
    pub records: RecordMetrics,

    /// Summary statistics for the junction_annotation subcommand.
    pub summary: SummaryResults,
}
