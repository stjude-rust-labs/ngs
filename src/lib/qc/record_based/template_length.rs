//! Functionality related to computing template lenght and related metrics.

use noodles_sam::alignment::Record;
use serde::{Deserialize, Serialize};

use crate::lib::{
    qc::{results, ComputationalLoad, RecordBasedQualityCheckFacet},
    utils::histogram::SimpleHistogram,
};

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SummaryMetrics {
    template_length_unknown_pct: f64,
    template_length_out_of_range_pct: f64,
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

/// Primary struct used to compile stats regarding template length. Within this
/// struct, the histogram represents the distribution of records with a particular
/// template length up to a certain threshold. Any records that fall outside of
/// that range are ignored (as tallied in the `ignored` field). Similarly,
/// records that are processed are tallied in the `processed` field.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct TemplateLengthFacet {
    // Histogram that represents the number of records that have a given
    // template length (up to the specified threshold).
    histogram: SimpleHistogram,
    records: RecordMetrics,
    summary: Option<SummaryMetrics>,
}

impl TemplateLengthFacet {
    /// Creates a new `TemplateLengthHistogram` with default values.
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            histogram: SimpleHistogram::zero_based_with_capacity(capacity),
            records: RecordMetrics {
                processed: 0,
                ignored: 0,
            },
            summary: None,
        }
    }

    /// Gets a value for the given bin within the histogram.
    #[allow(dead_code)]
    pub fn get(&self, bin: usize) -> usize {
        self.histogram.get(bin)
    }

    #[allow(dead_code)]
    /// Gets the number of processed records.
    pub fn get_processed_count(&self) -> usize {
        self.records.processed
    }

    /// Gets the number of ignored records.
    #[allow(dead_code)]
    pub fn get_ignored_count(&self) -> usize {
        self.records.ignored
    }
}

impl RecordBasedQualityCheckFacet for TemplateLengthFacet {
    fn name(&self) -> &'static str {
        "Template Length Metrics"
    }

    fn computational_load(&self) -> ComputationalLoad {
        ComputationalLoad::Light
    }

    fn process(&self, record: &Record) -> anyhow::Result<()> {
        let template_len = record.template_length() as usize;
        // match self.histogram.increment(template_len) {
        //     Ok(()) => self.records.processed += 1,
        //     Err(_) => self.records.ignored += 1,
        // }

        Ok(())
    }

    fn summarize(&mut self) -> anyhow::Result<()> {
        self.summary = Some(SummaryMetrics {
            template_length_unknown_pct: (self.histogram.get(0) as f64
                / (self.records.processed as f64 + self.records.ignored as f64))
                * 100.0,
            template_length_out_of_range_pct: (self.records.ignored as f64
                / (self.records.processed as f64 + self.records.ignored as f64))
                * 100.0,
        });

        Ok(())
    }

    fn aggregate(&self, results: &mut results::Results) {
        results.set_template_length(self.clone());
    }
}
