use serde::Serialize;

#[derive(Debug, Default, Serialize)]
pub struct ReadDesignation {
    pub primary: usize,
    pub secondary: usize,
    pub supplementary: usize,
}

#[derive(Debug, Default, Serialize)]
pub struct RecordMetrics {
    pub total: usize,
    pub unmapped: usize,
    pub duplicate: usize,
    pub designation: ReadDesignation,

    /**
     * Other Flagstat Metrics
     */
    pub primary_mapped: usize,
    pub primary_duplicate: usize,
    pub paired: usize,
    pub read_1: usize,
    pub read_2: usize,
    pub proper_pair: usize,
    pub singleton: usize,
    pub mate_mapped: usize,
}

#[derive(Debug, Default, Serialize)]
pub struct SummaryMetrics {
    pub duplication_pct: f64,
    pub unmapped_pct: f64,
}

#[derive(Debug, Default, Serialize)]
pub struct GeneralMetricsFacet {
    pub records: RecordMetrics,
    pub summary: Option<SummaryMetrics>,
}
