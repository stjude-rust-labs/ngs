//! Metrics related to the GC content quality control facet.

use serde::{Deserialize, Serialize};

use crate::utils::histogram::SimpleHistogram;

/// Metrics related to nucleobase counting.
#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct NucleobaseMetrics {
    /// Total number of nucleobases read as a 'G' or a 'C'.
    pub total_gc_count: usize,

    /// Total number of nucleobases read as an 'A' or a 'T'.
    pub total_at_count: usize,

    /// Total number of nucleobases which were not read as an 'A', a 'C', a 'G',
    /// or a 'T'.
    pub total_other_count: usize,
}

/// General metrics related to record counting.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct RecordMetrics {
    /// Number of records that have been processed by this struct.
    pub processed: usize,

    /// Number of records that were ignored because of their flags.
    pub ignored_flags: usize,

    /// Number of records that were ignored because they were too short.
    pub ignored_too_short: usize,
}

/// Summary statistics for the GC content control facet.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SummaryMetrics {
    /// Mean GC content for the given sample.
    pub gc_content_pct: f64,

    /// Percentage of records that were ignored because of flags.
    pub ignored_flags_pct: f64,

    /// Percentage of records that were ignored because they were too short.
    pub ignored_too_short_pct: f64,
}

/// Primary struct used to compile stats regarding GC content.
///
/// Within this struct, the histogram represents the number of reads which have
/// 0% GC content all the way up to 100% GC content. The other fields are for
/// counting the number of nucleobases which are G/C, the number of nucleobases
/// that are A/T, or the number of nucleobases that fall into other.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct GCContentMetrics {
    /// Histogram that represents the number of reads which have 0% GC content
    /// all the way up to 100% GC content.
    pub histogram: SimpleHistogram,

    /// Struct holding all of the nucleobase metrics.
    pub nucleobases: NucleobaseMetrics,

    /// Struct containing all of the status of processed/ignored records.
    pub records: RecordMetrics,

    /// Summary statistics for the GC content control facet.
    pub summary: Option<SummaryMetrics>,
}
