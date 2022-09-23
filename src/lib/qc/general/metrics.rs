use std::collections::HashMap;

use serde::Serialize;

#[derive(Clone, Debug, Default, Serialize)]
pub struct ReadDesignation {
    pub primary: usize,
    pub secondary: usize,
    pub supplementary: usize,
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct RecordMetrics {
    // Flagstat-like metrics
    pub total: usize,
    pub unmapped: usize,
    pub duplicate: usize,
    pub designation: ReadDesignation,
    pub primary_mapped: usize,
    pub primary_duplicate: usize,
    pub paired: usize,
    pub read_1: usize,
    pub read_2: usize,
    pub proper_pair: usize,
    pub singleton: usize,
    pub mate_mapped: usize,
    pub mate_reference_sequence_id_mismatch: usize,
    pub mate_reference_sequence_id_mismatch_hq: usize,

    //
    pub read_one_cigar_ops: HashMap<String, usize>,
    pub read_two_cigar_ops: HashMap<String, usize>,
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct SummaryMetrics {
    pub duplication_pct: f64,
    pub unmapped_pct: f64,
    pub mate_reference_sequence_id_mismatch_pct: f64,
    pub mate_reference_sequence_id_mismatch_hq_pct: f64,
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct GeneralMetricsFacet {
    pub records: RecordMetrics,
    pub summary: Option<SummaryMetrics>,
}
