use serde::{Deserialize, Serialize};

#[derive(Clone, Default, Serialize, Deserialize)]
pub struct ExonicTranslationRegionMetrics {
    pub utr_five_prime_count: usize,
    pub utr_three_prime_count: usize,
    pub coding_sequence_count: usize,
}

#[derive(Clone, Default, Serialize, Deserialize)]
pub struct GeneRegionMetrics {
    pub intergenic_count: usize,
    pub exonic_count: usize,
    pub intronic_count: usize,
}

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct RecordMetrics {
    pub processed: usize,
    pub ignored_flags: usize,
    pub ignored_nonprimary_chromosome: usize,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SummaryMetrics {
    pub ignored_flags_pct: f64,
    pub ignored_nonprimary_chromosome_pct: f64,
}

#[derive(Clone, Default, Serialize, Deserialize)]
pub struct Metrics {
    pub exonic_translation_regions: ExonicTranslationRegionMetrics,
    pub gene_regions: GeneRegionMetrics,
    pub records: RecordMetrics,
    pub summary: Option<SummaryMetrics>,
}
