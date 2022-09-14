//! Functionality related to computing template lenght and related metrics.

use std::{fs::File, io::Write, path::PathBuf};

use noodles_bam::lazy::Record;
use serde::Serialize;

use crate::utils::histogram::{BinOutOfBoundsError, SimpleHistogram};

use super::{ComputationalLoad, Error, QualityCheckFacet};

/// Primary struct used to compile stats regarding template length. Within this
/// struct, the histogram represents the distribution of records with a particular
/// template length up to a certain threshold. Any records that fall outside of
/// that range are ignored (as tallied in the `ignored` field). Similarly,
/// records that are processed are tallied in the `processed` field.
#[derive(Debug, Serialize)]
pub struct TemplateLengthFacet {
    // Histogram that represents the number of records that have a given
    // template length (up to the specified threshold).
    histogram: SimpleHistogram,

    // Number of records that were processed (and, as such, had template lengths
    // that fell within our histogram's range).
    processed: usize,

    // Number of records that were ignored (and, as such, had template lengths
    // that fell outside of our histogram's range).
    ignored: usize,

    template_length_unknown_pct: Option<f64>,
    template_length_out_of_range_pct: Option<f64>,
}

impl TemplateLengthFacet {
    /// Creates a new `TemplateLengthHistogram` with default values.
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            histogram: SimpleHistogram::zero_based_with_capacity(capacity),
            processed: 0,
            ignored: 0,
            template_length_unknown_pct: None,
            template_length_out_of_range_pct: None,
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
        self.processed
    }

    /// Gets the number of ignored records.
    #[allow(dead_code)]
    pub fn get_ignored_count(&self) -> usize {
        self.ignored
    }
}

impl QualityCheckFacet for TemplateLengthFacet {
    fn name(&self) -> &'static str {
        "Template Length Metrics"
    }

    fn computational_load(&self) -> ComputationalLoad {
        ComputationalLoad::Light
    }

    fn default(&self) -> bool {
        true
    }

    fn process(&mut self, record: &Record) -> Result<(), Error> {
        let template_len = record.template_length() as usize;
        match self.histogram.increment(template_len) {
            Ok(()) => self.processed += 1,
            Err(BinOutOfBoundsError) => self.ignored += 1,
        }

        Ok(())
    }

    fn summarize(&mut self) -> Result<(), Error> {
        self.template_length_unknown_pct = Some(
            (self.histogram.get(0) as f64 / (self.processed as f64 + self.ignored as f64)) * 100.0,
        );
        self.template_length_out_of_range_pct =
            Some((self.ignored as f64 / (self.processed as f64 + self.ignored as f64)) * 100.0);

        Ok(())
    }

    fn write(
        &self,
        output_prefix: String,
        directory: &std::path::Path,
    ) -> Result<(), std::io::Error> {
        let filename = output_prefix + ".template_lengths.json";
        let mut filepath = PathBuf::from(directory);
        filepath.push(filename);

        let mut file = File::create(filepath)?;
        let output = serde_json::to_string_pretty(&self).unwrap();
        file.write_all(output.as_bytes())?;
        Ok(())
    }
}
