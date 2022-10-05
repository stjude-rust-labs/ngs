use std::sync::atomic::{AtomicUsize, Ordering};

use serde::{Deserialize, Serialize};

use crate::lib::utils::histogram::{Histogram, MutexBackedHistogram};

//====================//
// Nucleobase Metrics //
//====================//

/// An atomic version of the nucleobase metrics struct. Useful when tallying these
/// statistics in a multi-threaded fashion.
pub struct AtomicNucleobaseMetrics {
    // Total number of nucleobases read as a 'G' or a 'C'.
    pub total_gc_count: AtomicUsize,

    // Total number of nucleobases read as an 'A' or a 'T'.
    pub total_at_count: AtomicUsize,

    // Total number of nucleobases which were not read as an 'A', a 'C', a 'G',
    // or a 'T'.
    pub total_other_count: AtomicUsize,
}

/// An equivalent, non-atomic version of the nucleobase metrics struct. Useful for
/// reading and writing to files.
#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct NucleobaseMetrics {
    // Total number of nucleobases read as a 'G' or a 'C'.
    total_gc_count: usize,

    // Total number of nucleobases read as an 'A' or a 'T'.
    total_at_count: usize,

    // Total number of nucleobases which were not read as an 'A', a 'C', a 'G',
    // or a 'T'.
    total_other_count: usize,
}

impl AtomicNucleobaseMetrics {
    /// The associated conversion from the nucleobase metrics struct to the
    /// non-atomic version. This will generally be called onceh tallying as
    /// completed and the results are ready to be written to a file.
    pub fn clone_nonatomic(&self) -> NucleobaseMetrics {
        NucleobaseMetrics {
            total_gc_count: self.total_gc_count.load(Ordering::SeqCst),
            total_at_count: self.total_at_count.load(Ordering::SeqCst),
            total_other_count: self.total_other_count.load(Ordering::SeqCst),
        }
    }
}

//================//
// Record Metrics //
//================//

/// An atomic version of the record metrics struct. Useful when tallying these
/// statistics in a multi-threaded fashion.
pub struct AtomicRecordMetrics {
    // Number of records that have been processed by this struct.
    pub processed: AtomicUsize,

    // Number of records that were ignored because of their flags.
    pub ignored_flags: AtomicUsize,

    // Number of records that were ignored because they were too short.
    pub ignored_too_short: AtomicUsize,
}

/// An equivalent, non-atomic version of the record metrics struct. Useful for
/// reading and writing to files.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct RecordMetrics {
    // Number of records that have been processed by this struct.
    processed: usize,

    // Number of records that were ignored because of their flags.
    ignored_flags: usize,

    // Number of records that were ignored because they were too short.
    ignored_too_short: usize,
}

impl AtomicRecordMetrics {
    /// The associated conversion from the record metrics struct to the
    /// non-atomic version. This will generally be called onceh tallying as
    /// completed and the results are ready to be written to a file.
    pub fn clone_nonatomic(&self) -> RecordMetrics {
        RecordMetrics {
            processed: self.processed.load(Ordering::SeqCst),
            ignored_flags: self.ignored_flags.load(Ordering::SeqCst),
            ignored_too_short: self.ignored_too_short.load(Ordering::SeqCst),
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
    // Mean GC content for the given sample.
    pub gc_content_pct: f64,

    // Percentage of records that were ignored because of flags.
    pub ignored_flags_pct: f64,

    // Percentage of records that were ignored because they were too short.
    pub ignored_too_short_pct: f64,
}

//==============//
// Main Metrics //
//==============//

pub struct AtomicGCContentMetrics {
    // Histogram that represents the number of reads which have 0% GC content
    // all the way up to 100% GC content.
    pub histogram: MutexBackedHistogram,

    // Struct holding all of the nucleobase metrics.
    pub nucleobases: AtomicNucleobaseMetrics,

    // Struct containing all of the status of processed/ignored records.
    pub records: AtomicRecordMetrics,

    pub summary: Option<SummaryMetrics>,
}

impl Default for AtomicGCContentMetrics {
    fn default() -> Self {
        Self {
            histogram: MutexBackedHistogram::zero_based_with_capacity(100),
            nucleobases: AtomicNucleobaseMetrics {
                total_at_count: AtomicUsize::new(0),
                total_gc_count: AtomicUsize::new(0),
                total_other_count: AtomicUsize::new(0),
            },
            records: AtomicRecordMetrics {
                processed: AtomicUsize::new(0),
                ignored_flags: AtomicUsize::new(0),
                ignored_too_short: AtomicUsize::new(0),
            },
            summary: None,
        }
    }
}

/// Primary struct used to compile stats regarding GC content. Within this
/// struct, the histogram represents the number of reads which have 0% GC
/// content all the way up to 100% GC content. The other fields are for counting
/// the number of nucleobases which are G/C, the number of nucleobases that are
/// A/T, or the number of nucleobases that fall into other.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct GCContentMetrics {
    // Histogram that represents the number of reads which have 0% GC content
    // all the way up to 100% GC content.
    pub histogram: Histogram,

    // Struct holding all of the nucleobase metrics.
    nucleobases: NucleobaseMetrics,

    // Struct containing all of the status of processed/ignored records.
    records: RecordMetrics,

    summary: Option<SummaryMetrics>,
}

impl AtomicGCContentMetrics {
    pub fn clone_nonatomic(&self) -> GCContentMetrics {
        GCContentMetrics {
            histogram: self.histogram.clone_nonmutex(),
            nucleobases: self.nucleobases.clone_nonatomic(),
            records: self.records.clone_nonatomic(),
            summary: self.summary.clone(),
        }
    }
}
