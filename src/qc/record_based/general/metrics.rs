//! Metrics related to the General quality control facet.

use std::collections::HashMap;

use serde::Deserialize;
use serde::Serialize;

/// Metrics related to tallying read designations (primary, secondary,
/// supplementary).
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct ReadDesignationMetrics {
    /// Number of records that are marked as primary.
    pub primary: usize,

    /// Number of records that are marked as secondary.
    pub secondary: usize,

    /// Number of records that are marked as supplementary.
    pub supplementary: usize,
}

/// General purpose record metrics.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct RecordMetrics {
    /// The total number of records within the file.
    pub total: usize,

    /// The total number of records marked as unmapped (`0x4`) within the file.
    pub unmapped: usize,

    /// The number of records marked as duplicate (`0x400`) within the file.
    pub duplicate: usize,

    /// The number of _primary_, _secondary_, and _supplementary_ records in the
    /// file respectively.
    ///
    /// * If a read is marked as secondary (`0x100`), then the read is counted
    ///   as secondary.
    /// * Else, if a read is marked as supplementary (`0x800`), then the read is
    ///   counted as supplementary.
    /// * Else, the read is counted as primary.
    pub designation: ReadDesignationMetrics,

    /// The number of records that are counted as primary and and are marked as
    /// mapped (`!0x4`).
    pub primary_mapped: usize,

    /// The number of records that counted as primary and are marked as
    /// duplicate (`0x400`).
    pub primary_duplicate: usize,

    /// The number of records that are designated as primary and marked as
    /// segmented (`0x01`).
    pub paired: usize,

    /// The number of records that are designated as primary, marked as
    /// segmented (`0x01`), and marked as being the first record within a
    /// segment (`0x40`).
    pub read_1: usize,

    /// The number of records that are designated as primary, marked as
    /// segmented (`0x01`), and marked as being the last record within a segment
    /// (`0x80`).
    pub read_2: usize,

    /// The number of records that are designated as primary, marked as segmented
    /// (`0x01`), marked as mapped (`!0x04`), and properly aligned (`0x2`).
    pub proper_pair: usize,

    /// The number of records that are designated as primary, marked as
    /// segmented (`0x01`), marked as mapped (`!0x04`), and marked as mate is
    /// unmapped (`0x08`).
    pub singleton: usize,

    /// The number of records that are designated as primary, marked as
    /// segmented (`0x01`), marked as mapped (`!0x04`), and marked as mate is
    /// mapped (`!0x08`).
    pub mate_mapped: usize,

    /// The number of records that are designated as primary, marked as
    /// segmented (`0x01`), marked as mapped (`!0x04`), marked as mate is mapped
    /// (`!0x08`), but the sequence id that the mate is matched to is different
    /// that the record being examined.
    pub mate_reference_sequence_id_mismatch: usize,

    /// The number of records that are designated as primary, marked as
    /// segmented (`0x01`), marked as mapped (`!0x04`), marked as mate is mapped
    /// (`!0x08`), but the sequence id that the mate is matched to is different
    /// that the record being examined **and** the mapping quality of the
    /// current record is greater than 5.
    pub mate_reference_sequence_id_mismatch_hq: usize,
}

/// Metrics related the to the CIGAR string.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct CigarMetrics {
    /// Count of each of the various CIGAR operations present within all read
    /// ones.
    pub read_one_cigar_ops: HashMap<String, usize>,

    /// Count of each of the various CIGAR operations present within all read
    /// twos.
    pub read_two_cigar_ops: HashMap<String, usize>,
}

/// Summary metrics related to the General quality control facet.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct SummaryMetrics {
    /// Percentage of records that were marked as duplicate.
    pub duplication_pct: f64,

    /// Percentage of records that were marked as mapped.
    pub mapped_pct: f64,

    /// Percentage of records that qualified as
    /// `mate_reference_sequence_id_mismatch`.
    pub mate_reference_sequence_id_mismatch_pct: f64,

    /// Percentage of records that qualified as
    /// `mate_reference_sequence_id_mismatch_hq`.
    pub mate_reference_sequence_id_mismatch_hq_pct: f64,
}

/// Aggregate struct for all minor metrics structs for the General quality
/// control facet.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct GeneralMetrics {
    /// General metrics pertaining the the records counted within the file.
    pub records: RecordMetrics,

    /// Metrics related to CIGAR string pileups for read ones and read twos.
    pub cigar: CigarMetrics,

    /// Summary statistics for the General quality control facet.
    pub summary: Option<SummaryMetrics>,
}
