//! Functionality related to computing template lenght and related metrics.

use noodles_bam::lazy::Record;
use serde::Serialize;

use crate::utils::histogram::{BinOutOfBounds, SimpleHistogram};

/// Primary struct used to compile stats regarding template length. Within this
/// struct, the histogram represents the distribution of records with a particular
/// template length up to a certain threshold. Any records that fall outside of
/// that range are ignored (as tallied in the `ignored` field). Similarly,
/// records that are processed are tallied in the `processed` field.
#[derive(Debug, Serialize)]
pub struct TemplateLengthHistogram {
    // Histogram that represents the number of records that have a given
    // template length (up to the specified threshold).
    histogram: SimpleHistogram,

    // Number of records that were processed (and, as such, had template lengths
    // that fell within our histogram's range).
    processed: usize,

    // Number of records that were ignored (and, as such, had template lengths
    // that fell outside of our histogram's range).
    ignored: usize,
}

impl TemplateLengthHistogram {
    /// Creates a new `TemplateLengthHistogram` with default values.
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            histogram: SimpleHistogram::zero_based_with_capacity(capacity),
            processed: 0,
            ignored: 0,
        }
    }

    /// Gets a value for the given bin within the histogram.
    pub fn get(&self, bin: usize) -> usize {
        self.histogram.get(bin)
    }

    #[allow(dead_code)]
    /// Gets the number of processed records.
    pub fn get_processed_count(&self) -> usize {
        self.processed
    }

    /// Gets the number of ignored records.
    pub fn get_ignored_count(&self) -> usize {
        self.ignored
    }

    /// Processes a record and updates the struct accordingly.
    pub fn process(&mut self, record: &Record) {
        let template_len = record.template_length() as usize;
        match self.histogram.increment(template_len) {
            Ok(()) => self.processed += 1,
            Err(BinOutOfBounds) => self.ignored += 1,
        }
    }
}
