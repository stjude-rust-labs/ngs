use std::sync::atomic::{AtomicUsize, Ordering};

use serde::{Deserialize, Serialize};

//===================================//
// Exonic Translation Region Metrics //
//===================================//

/// An atomic version of the exonic translation region metrics struct. Useful
/// when tallying these statistics in a multi-threaded fashion.
pub struct AtomicExonicTranslationRegionMetrics {
    pub utr_five_prime_count: AtomicUsize,
    pub utr_three_prime_count: AtomicUsize,
    pub coding_sequence_count: AtomicUsize,
}

/// An equivalent, non-atomic version of the exonic translation region metrics
/// struct. Useful for reading and writing to files.
#[derive(Serialize, Deserialize)]
pub struct ExonicTranslationRegionMetrics {
    utr_five_prime_count: usize,
    utr_three_prime_count: usize,
    coding_sequence_count: usize,
}

impl AtomicExonicTranslationRegionMetrics {
    /// The associated conversion from the atomic exonic translation region metrics
    /// struct to the non-atomic version. This will generally be called onceh
    /// tallying as completed and the results are ready to be written to a file.
    pub fn clone_nonatomic(&self) -> ExonicTranslationRegionMetrics {
        ExonicTranslationRegionMetrics {
            utr_five_prime_count: self.utr_five_prime_count.load(Ordering::SeqCst),
            utr_three_prime_count: self.utr_three_prime_count.load(Ordering::SeqCst),
            coding_sequence_count: self.coding_sequence_count.load(Ordering::SeqCst),
        }
    }
}

//=====================//
// Gene Region Metrics //
//=====================//

/// An atomic version of the gene region metrics struct. Useful when tallying
/// these statistics in a multi-threaded fashion.
pub struct AtomicGeneRegionMetrics {
    pub intergenic_count: AtomicUsize,
    pub exonic_count: AtomicUsize,
    pub intronic_count: AtomicUsize,
}

/// An equivalent, non-atomic version of the gene region metrics struct. Useful
/// for reading and writing to files.
#[derive(Debug, Serialize, Deserialize)]
pub struct GeneRegionMetrics {
    intergenic_count: usize,
    exonic_count: usize,
    intronic_count: usize,
}

impl AtomicGeneRegionMetrics {
    /// The associated conversion from the atomic gene region metrics struct to the
    /// non-atomic version. This will generally be called onceh tallying as
    /// completed and the results are ready to be written to a file.
    pub fn clone_nonatomic(&self) -> GeneRegionMetrics {
        GeneRegionMetrics {
            intergenic_count: self.intergenic_count.load(Ordering::SeqCst),
            exonic_count: self.exonic_count.load(Ordering::SeqCst),
            intronic_count: self.intronic_count.load(Ordering::SeqCst),
        }
    }
}

//================//
// Record Metrics //
//================//

/// An atomic version of the record metrics struct. Useful when tallying these
/// statistics in a multi-threaded fashion.
pub struct AtomicRecordMetrics {
    pub processed: AtomicUsize,
    pub ignored_flags: AtomicUsize,
    pub ignored_nonprimary_chromosome: AtomicUsize,
}

/// An equivalent, non-atomic version of the record metrics struct. Useful for
/// reading and writing to files.
#[derive(Debug, Serialize, Deserialize)]
pub struct RecordMetrics {
    processed: usize,
    ignored_flags: usize,
    ignored_nonprimary_chromosome: usize,
}

impl AtomicRecordMetrics {
    /// The associated conversion from the atomic gene region metrics struct to the
    /// non-atomic version. This will generally be called onceh tallying as
    /// completed and the results are ready to be written to a file.
    pub fn clone_nonatomic(&self) -> RecordMetrics {
        RecordMetrics {
            processed: self.processed.load(Ordering::SeqCst),
            ignored_flags: self.ignored_flags.load(Ordering::SeqCst),
            ignored_nonprimary_chromosome: self
                .ignored_nonprimary_chromosome
                .load(Ordering::SeqCst),
        }
    }
}

//=================//
// Summary Metrics //
//=================//

/// A struct which hold summary data for this quality control facet. An atomic
/// version of the struct is not needed, as these are calculated after all of
/// the results have been tallied (in a single threaded environment).
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SummaryMetrics {
    pub ignored_flags_pct: f64,
    pub ignored_nonprimary_chromosome_pct: f64,
}

//==============//
// Main Metrics //
//==============//

/// An atomic version of the main metrics struct. This struct contains all
/// atomic versions of the underlying structs. This will be the main vehicle
/// within which metrics are tallied.
pub struct AtomicMetrics {
    pub exonic_translation_regions: AtomicExonicTranslationRegionMetrics,
    pub gene_regions: AtomicGeneRegionMetrics,
    pub records: AtomicRecordMetrics,
    pub summary: Option<SummaryMetrics>,
}

/// An equivalent, non-atomic struct containing all metrics for this quality
/// check facet. Useful when reading from or writing to files.
#[derive(Serialize, Deserialize)]
pub struct Metrics {
    pub exonic_translation_regions: ExonicTranslationRegionMetrics,
    pub gene_regions: GeneRegionMetrics,
    pub records: RecordMetrics,
    pub summary: Option<SummaryMetrics>,
}

impl AtomicMetrics {
    /// Creates a new, empty struct for tallying metrics within.
    pub fn new() -> Self {
        AtomicMetrics {
            exonic_translation_regions: AtomicExonicTranslationRegionMetrics {
                utr_five_prime_count: AtomicUsize::new(0),
                utr_three_prime_count: AtomicUsize::new(0),
                coding_sequence_count: AtomicUsize::new(0),
            },
            gene_regions: AtomicGeneRegionMetrics {
                intergenic_count: AtomicUsize::new(0),
                exonic_count: AtomicUsize::new(0),
                intronic_count: AtomicUsize::new(0),
            },
            records: AtomicRecordMetrics {
                processed: AtomicUsize::new(0),
                ignored_flags: AtomicUsize::new(0),
                ignored_nonprimary_chromosome: AtomicUsize::new(0),
            },
            summary: None,
        }
    }

    /// The associated conversion from the main metrics struct to the non-atomic
    /// version. This will generally be called onceh tallying as completed and
    /// the results are ready to be written to a file.
    pub fn clone_nonatomic(&self) -> Metrics {
        Metrics {
            exonic_translation_regions: self.exonic_translation_regions.clone_nonatomic(),
            gene_regions: self.gene_regions.clone_nonatomic(),
            records: self.records.clone_nonatomic(),
            summary: self.summary.clone(),
        }
    }
}
