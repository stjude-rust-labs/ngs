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
