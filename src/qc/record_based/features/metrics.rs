//! Metrics related to the Features quality control facet.

use serde::Deserialize;
use serde::Serialize;

/// Metrics related to the tallying of records in exonic translation regions
/// (five prime UTR regions, three prime UTR regions, coding sequence regions).
#[derive(Clone, Default, Serialize, Deserialize)]
pub struct ExonicTranslationRegionMetrics {
    /// Count for the number of records falling within a five prime UTR region that have
    /// been detected within this file.
    pub utr_five_prime_count: usize,

    /// Count for the number of records falling within a three prime UTR region
    /// that have been detected within this file.
    pub utr_three_prime_count: usize,

    /// Count for the number of records falling within a coding sequence region
    /// that have been detected within this file.
    pub coding_sequence_count: usize,
}

/// Metrics related to the tallying of records in gene regions (intronic,
/// exonic, intergenic).
#[derive(Clone, Default, Serialize, Deserialize)]
pub struct GeneRegionMetrics {
    /// Count for the number of records falling within an intergenic region that
    /// have been detected within this file.
    pub intergenic_count: usize,

    /// Count for the number of records falling within an exonic region that
    /// have been detected within this file.
    pub exonic_count: usize,

    /// Count for the number of records falling within an intronic region that
    /// have been detected within this file.
    pub intronic_count: usize,
}

/// General record metrics that are tallied as a part of the Features quality
/// control facet.
#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct RecordMetrics {
    /// The number of records that have been processed.
    pub processed: usize,

    /// The number of records that have been ignored because of bad flags.
    pub ignored_flags: usize,

    /// The number of records that have been ignored because they were not
    /// aligned to a primary chromosome.
    pub ignored_nonprimary_chromosome: usize,
}

/// Summary statistics for the Features quality control facet.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SummaryMetrics {
    /// The percentage of records that were ignored because of bad flags.
    pub ignored_flags_pct: f64,
    /// The percentage of records that were ignored because they were not
    /// aligned to a primary chromosome.
    pub ignored_nonprimary_chromosome_pct: f64,
}

/// Main metrics struct. This struct aggregates all of the minor metrics structs
/// outlined in this file so they can be kept track of as a unit.
#[derive(Clone, Default, Serialize, Deserialize)]
pub struct Metrics {
    /// Exonic translation region metrics.
    pub exonic_translation_regions: ExonicTranslationRegionMetrics,

    /// Gene region metrics.
    pub gene_regions: GeneRegionMetrics,

    /// General record metrics.
    pub records: RecordMetrics,

    /// Summary statistics for the Features quality control facet.
    pub summary: Option<SummaryMetrics>,
}
