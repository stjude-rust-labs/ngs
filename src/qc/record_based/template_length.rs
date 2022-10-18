//! Functionality related to the Template Length quality control facet.

use noodles::sam::alignment::Record;
use serde::Deserialize;
use serde::Serialize;

use crate::qc::results;
use crate::qc::ComputationalLoad;
use crate::qc::RecordBasedQualityControlFacet;
use crate::utils::histogram::Histogram;

/// Summary statistics for the Template Length quality control facet.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SummaryMetrics {
    /// Percentage of records for which the template length was unknown (zero).
    pub template_length_unknown_pct: f64,

    /// Percentage of records for which the template length was greater than our
    /// histogram could support.
    pub template_length_out_of_range_pct: f64,
}

/// General metrics regarding records collected in the quality control facet.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct RecordMetrics {
    /// Number of records that were processed (and, as such, had template lengths
    /// that fell within our histogram's range).
    pub processed: usize,

    /// Number of records that were ignored (and, as such, had template lengths
    /// that fell outside of our histogram's range).
    pub ignored: usize,
}

/// Main struct for the Template Length quality control facet.
///
/// Within this struct, the histogram represents the distribution of records
/// with a particular template length up to a certain threshold. Any records
/// that fall outside of that range are ignored (as tallied in the `ignored`
/// field). Similarly, records that are processed are tallied in the `processed`
/// field.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct TemplateLengthFacet {
    /// Histogram that represents the number of records that have a given
    /// template length (up to the specified threshold).
    pub histogram: Histogram,

    /// General record metrics
    pub records: RecordMetrics,

    /// Summary statistics for the Template Length quality control facet.
    pub summary: Option<SummaryMetrics>,
}

impl TemplateLengthFacet {
    /// Creates a new [`TemplateLengthFacet`] with a specified capacity and
    /// otherwise default values.
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            histogram: Histogram::zero_based_with_capacity(capacity),
            records: RecordMetrics {
                processed: 0,
                ignored: 0,
            },
            summary: None,
        }
    }
}

impl RecordBasedQualityControlFacet for TemplateLengthFacet {
    fn name(&self) -> &'static str {
        "Template Length"
    }

    fn computational_load(&self) -> ComputationalLoad {
        ComputationalLoad::Light
    }

    fn process(&mut self, record: &Record) -> anyhow::Result<()> {
        let template_len = record.template_length() as usize;
        match self.histogram.increment(template_len) {
            Ok(()) => self.records.processed += 1,
            Err(_) => self.records.ignored += 1,
        }

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
        results.template_length = Some(self.clone());
    }
}
