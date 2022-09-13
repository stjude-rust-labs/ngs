use serde::Serialize;

#[derive(Serialize)]
pub struct ExonicTranslationMetrics {
    utr_five_prime_count: usize,
    utr_three_prime_count: usize,
    coding_sequence_count: usize,
}

impl ExonicTranslationMetrics {
    pub fn inc_utr_five_prime_count(&mut self) {
        self.utr_five_prime_count += 1;
    }

    pub fn inc_utr_three_prime_count(&mut self) {
        self.utr_three_prime_count += 1;
    }

    pub fn inc_coding_sequence_count(&mut self) {
        self.coding_sequence_count += 1;
    }
}

#[derive(Serialize)]
pub struct GeneRegionMetrics {
    intergenic_count: usize,
    exonic_count: usize,
    intronic_count: usize,
}

impl GeneRegionMetrics {
    pub fn inc_intergenic_count(&mut self) {
        self.intergenic_count += 1;
    }

    pub fn inc_exonic_count(&mut self) {
        self.exonic_count += 1;
    }

    pub fn inc_intronic_count(&mut self) {
        self.intronic_count += 1;
    }
}

#[derive(Serialize)]
pub struct Metrics {
    pub exonic_translations: ExonicTranslationMetrics,
    pub gene_regions: GeneRegionMetrics,
}

impl Metrics {
    pub fn new() -> Self {
        Metrics {
            exonic_translations: ExonicTranslationMetrics {
                utr_five_prime_count: 0,
                utr_three_prime_count: 0,
                coding_sequence_count: 0,
            },
            gene_regions: GeneRegionMetrics {
                intergenic_count: 0,
                exonic_count: 0,
                intronic_count: 0,
            },
        }
    }
}
