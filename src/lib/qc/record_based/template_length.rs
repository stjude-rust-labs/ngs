//! Functionality related to computing template lenght and related metrics.

use std::sync::atomic::{AtomicUsize, Ordering};

use noodles_sam::alignment::Record;
use serde::{Deserialize, Serialize};

use crate::lib::{
    qc::{results, ComputationalLoad, RecordBasedQualityCheckFacet},
    utils::histogram::{Histogram, MutexBackedHistogram},
};

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SummaryMetrics {
    template_length_unknown_pct: f64,
    template_length_out_of_range_pct: f64,
}

//================//
// Record Metrics //
//================//

#[derive(Debug)]
pub struct AtomicRecordMetrics {
    // Number of records that were processed (and, as such, had template lengths
    // that fell within our histogram's range).
    pub processed: AtomicUsize,

    // Number of records that were ignored (and, as such, had template lengths
    // that fell outside of our histogram's range).
    pub ignored: AtomicUsize,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct RecordMetrics {
    // Number of records that were processed (and, as such, had template lengths
    // that fell within our histogram's range).
    processed: usize,

    // Number of records that were ignored (and, as such, had template lengths
    // that fell outside of our histogram's range).
    ignored: usize,
}

impl AtomicRecordMetrics {
    pub fn clone_nonatomic(&self) -> RecordMetrics {
        RecordMetrics {
            processed: self.processed.load(Ordering::SeqCst),
            ignored: self.ignored.load(Ordering::SeqCst),
        }
    }
}

//=========================//
// Template Length Metrics //
//=========================//

#[derive(Debug)]
pub struct AtomicTemplateLengthMetrics {
    // Histogram that represents the number of records that have a given
    // template length (up to the specified threshold).
    histogram: MutexBackedHistogram,
    records: AtomicRecordMetrics,
    summary: Option<SummaryMetrics>,
}

/// Primary struct used to compile stats regarding template length. Within this
/// struct, the histogram represents the distribution of records with a particular
/// template length up to a certain threshold. Any records that fall outside of
/// that range are ignored (as tallied in the `ignored` field). Similarly,
/// records that are processed are tallied in the `processed` field.
#[derive(Debug, Serialize, Deserialize)]
pub struct TemplateLengthMetrics {
    // Histogram that represents the number of records that have a given
    // template length (up to the specified threshold).
    histogram: Histogram,
    records: RecordMetrics,
    summary: Option<SummaryMetrics>,
}

impl AtomicTemplateLengthMetrics {
    pub fn clone_nonatomic(&self) -> TemplateLengthMetrics {
        TemplateLengthMetrics {
            histogram: self.histogram.clone_nonmutex(),
            records: self.records.clone_nonatomic(),
            summary: self.summary.clone(),
        }
    }
}

#[derive(Debug)]
pub struct TemplateLengthFacet {
    metrics: AtomicTemplateLengthMetrics,
}

impl TemplateLengthFacet {
    /// Creates a new `TemplateLengthHistogram` with default values.
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            metrics: AtomicTemplateLengthMetrics {
                histogram: MutexBackedHistogram::zero_based_with_capacity(capacity),
                records: AtomicRecordMetrics {
                    processed: AtomicUsize::new(0),
                    ignored: AtomicUsize::new(0),
                },
                summary: None,
            },
        }
    }
}

impl RecordBasedQualityCheckFacet for TemplateLengthFacet {
    fn name(&self) -> &'static str {
        "Template Length Metrics"
    }

    fn computational_load(&self) -> ComputationalLoad {
        ComputationalLoad::Light
    }

    fn process(&mut self, record: &Record) -> anyhow::Result<()> {
        let template_len = record.template_length() as usize;

        match self.metrics.histogram.increment(template_len) {
            Ok(()) => self
                .metrics
                .records
                .processed
                .fetch_add(1, Ordering::SeqCst),
            Err(_) => self.metrics.records.ignored.fetch_add(1, Ordering::SeqCst),
        };

        Ok(())
    }

    fn summarize(&mut self) -> anyhow::Result<()> {
        let processed = self.metrics.records.processed.load(Ordering::SeqCst);
        let ignored = self.metrics.records.ignored.load(Ordering::SeqCst);

        self.metrics.summary = Some(SummaryMetrics {
            template_length_unknown_pct: (self.metrics.histogram.get(0) as f64
                / (processed as f64 + ignored as f64))
                * 100.0,
            template_length_out_of_range_pct: (ignored as f64
                / (processed as f64 + ignored as f64))
                * 100.0,
        });

        Ok(())
    }

    fn aggregate(&self, results: &mut results::Results) {
        results.set_template_length(self.metrics.clone_nonatomic());
    }
}
