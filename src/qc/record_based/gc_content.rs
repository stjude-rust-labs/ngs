//! Functionality related to the GC content quality control facet.

pub mod metrics;

use noodles::sam;
use rand::prelude::*;
use sam::alignment::Record;
use sam::record::sequence::Base;

use crate::qc::results;
use crate::qc::ComputationalLoad;
use crate::qc::RecordBasedQualityControlFacet;
use crate::utils::histogram::Histogram;

use self::metrics::GCContentMetrics;
use self::metrics::SummaryMetrics;

/// Truncates reads that are longer than this value by randomly selecting a
/// substring of this size.
pub const TRUNCATION_LENGTH: usize = 100;

/// Main struct for the GC content quality control facet.
#[derive(Default)]
pub struct GCContentFacet {
    /// The main metric counting struct.
    pub metrics: GCContentMetrics,
}

impl RecordBasedQualityControlFacet for GCContentFacet {
    fn name(&self) -> &'static str {
        "GC Content"
    }

    fn computational_load(&self) -> ComputationalLoad {
        ComputationalLoad::Light
    }

    fn process(&mut self, record: &Record) -> anyhow::Result<()> {
        // (1) Check the record's flags. If any of the flags aren't to our
        // liking, then we reject the record as an ignored flag record.
        let flags = record.flags();
        if flags.is_duplicate() || flags.is_secondary() {
            self.metrics.records.ignored_flags += 1;
            return Ok(());
        };

        // (2) Convert the BAM record to a SAM record so we can determine the
        // nucleobases and count up the A's, C's, G's, and T's. TODO: this could
        // be done strictly from the BAM without parsing into SAM.
        let sequence = record.sequence();
        let nucleobases = sequence.as_ref();
        let sequence_length = nucleobases.len();

        // (3) Checks whether the record is too short to check the GC bias for.
        // If the record is too short, it can only be a subset of the full range
        // of percentages from 0%-100%. We don't want that messing with our
        // results—all records should have a fair chance to generate between 0%
        // to 100% GC content.
        if sequence_length < TRUNCATION_LENGTH {
            self.metrics.records.ignored_too_short += 1;
            return Ok(());
        }

        // (4) Truncate the read TRUNCATION_LENGTH (generally 100 nucleobases)
        // so that we have a uniform chance of generating between 0%-100%. We
        // choose a random starting point within the record so as to not
        // introduce a bias towards the start of our record.
        let mut gc_this_read = 0usize;
        let offset = if TRUNCATION_LENGTH < sequence_length {
            let max_offset = sequence_length - TRUNCATION_LENGTH;
            ThreadRng::default().gen_range(0..max_offset)
        } else {
            0
        };

        // (5) Count up the A's, C's, G's, and T's.
        for i in 0..TRUNCATION_LENGTH {
            let nucleobase = nucleobases[offset + i];
            match nucleobase {
                Base::C | Base::G => {
                    gc_this_read += 1;
                    self.metrics.nucleobases.total_gc_count += 1;
                }
                Base::A | Base::T => self.metrics.nucleobases.total_at_count += 1,
                _ => self.metrics.nucleobases.total_other_count += 1,
            }
        }

        // (6) Calculate the GC content for this read and increment the
        // histogram accordingly.
        let gc_content_this_read_pct =
            ((gc_this_read as f64 / TRUNCATION_LENGTH as f64) * 100.0).round() as usize;
        self.metrics
            .histogram
            .increment(gc_content_this_read_pct)
            .unwrap();
        self.metrics.records.processed += 1;

        Ok(())
    }

    fn summarize(&mut self) -> anyhow::Result<()> {
        self.metrics.summary = Some(SummaryMetrics {
            gc_content_pct: (self.metrics.nucleobases.total_gc_count as f64
                / (self.metrics.nucleobases.total_gc_count
                    + self.metrics.nucleobases.total_at_count
                    + self.metrics.nucleobases.total_other_count) as f64)
                * 100.0,
            ignored_flags_pct: (self.metrics.records.ignored_flags as f64
                / (self.metrics.records.ignored_flags
                    + self.metrics.records.ignored_too_short
                    + self.metrics.records.processed) as f64)
                * 100.0,
            ignored_too_short_pct: (self.metrics.records.ignored_too_short as f64
                / (self.metrics.records.ignored_flags
                    + self.metrics.records.ignored_too_short
                    + self.metrics.records.processed) as f64)
                * 100.0,
        });

        Ok(())
    }

    fn aggregate(&self, results: &mut results::Results) {
        results.gc_content = Some(self.metrics.clone());
    }
}

impl Default for GCContentMetrics {
    fn default() -> Self {
        Self {
            histogram: Histogram::zero_based_with_capacity(100),
            nucleobases: Default::default(),
            records: Default::default(),
            summary: Default::default(),
        }
    }
}

#[cfg(test)]

mod tests {
    use super::*;

    #[test]
    pub fn it_defaults_with_zero_based_100_capacity_histogram() {
        let default = GCContentFacet::default();
        assert_eq!(default.metrics.histogram.range_start(), 0);
        assert_eq!(default.metrics.histogram.range_stop(), 100);
        assert_eq!(default.metrics.histogram.range_len(), 101);
    }
}
