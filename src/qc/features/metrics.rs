use serde::Serialize;

#[derive(Serialize)]
pub struct ExonicTranslationRegionMetrics {
    pub utr_five_prime_count: usize,
    pub utr_three_prime_count: usize,
    pub coding_sequence_count: usize,
}

#[derive(Serialize)]
pub struct GeneRegionMetrics {
    pub intergenic_count: usize,
    pub exonic_count: usize,
    pub intronic_count: usize,
}

#[derive(Debug, Serialize)]
pub struct RecordMetrics {
    pub processed: usize,
    pub ignored_flags: usize,
    pub ignored_nonprimary_chromosome: usize,
}

#[derive(Debug, Serialize)]
pub struct SummaryMetrics {
    pub ignored_flags_pct: f64,
    pub ignored_nonprimary_chromosome_pct: f64,
}

#[derive(Serialize)]
pub struct Metrics {
    pub exonic_translation_regions: ExonicTranslationRegionMetrics,
    pub gene_regions: GeneRegionMetrics,
    pub records: RecordMetrics,
    pub summary: Option<SummaryMetrics>,
}

impl Metrics {
    pub fn new() -> Self {
        Metrics {
            exonic_translation_regions: ExonicTranslationRegionMetrics {
                utr_five_prime_count: 0,
                utr_three_prime_count: 0,
                coding_sequence_count: 0,
            },
            gene_regions: GeneRegionMetrics {
                intergenic_count: 0,
                exonic_count: 0,
                intronic_count: 0,
            },
            records: RecordMetrics {
                processed: 0,
                ignored_flags: 0,
                ignored_nonprimary_chromosome: 0,
            },
            summary: None,
        }
    }
}
