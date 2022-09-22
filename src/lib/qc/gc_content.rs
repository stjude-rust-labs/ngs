//! Functionality related to computing GC content and related metrics.

use noodles_bam::lazy::Record;
use noodles_sam as sam;
use rand::prelude::*;
use sam::record::sequence::Base;
use serde::Serialize;

use crate::lib::utils::histogram::SimpleHistogram;

use super::{results, ComputationalLoad, Error, RecordBasedQualityCheckFacet};

const TRUNCATION_LENGTH: usize = 100;

#[derive(Clone, Debug, Default, Serialize)]
pub struct NucleobaseMetrics {
    // Total number of nucleobases read as a 'G' or a 'C'.
    total_gc_count: usize,

    // Total number of nucleobases read as an 'A' or a 'T'.
    total_at_count: usize,

    // Total number of nucleobases which were not read as an 'A', a 'C', a 'G',
    // or a 'T'.
    total_other_count: usize,
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct RecordMetrics {
    // Number of records that have been processed by this struct.
    processed: usize,

    // Number of records that were ignored because of their flags.
    ignored_flags: usize,

    // Number of records that were ignored because they were too short.
    ignored_too_short: usize,
}

#[derive(Clone, Debug, Serialize)]
pub struct SummaryMetrics {
    // Mean GC content for the given sample.
    gc_content_pct: f64,

    // Percentage of records that were ignored because of flags.
    ignored_flags_pct: f64,

    // Percentage of records that were ignored because they were too short.
    ignored_too_short_pct: f64,
}

/// Primary struct used to compile stats regarding GC content. Within this
/// struct, the histogram represents the number of reads which have 0% GC
/// content all the way up to 100% GC content. The other fields are for counting
/// the number of nucleobases which are G/C, the number of nucleobases that are
/// A/T, or the number of nucleobases that fall into other.
#[derive(Clone, Debug, Default, Serialize)]
pub struct GCContentMetrics {
    // Histogram that represents the number of reads which have 0% GC content
    // all the way up to 100% GC content.
    histogram: SimpleHistogram,

    // Struct holding all of the nucleobase metrics.
    pub nucleobases: NucleobaseMetrics,

    // Struct containing all of the status of processed/ignored records.
    pub records: RecordMetrics,

    summary: Option<SummaryMetrics>,
}

#[derive(Default)]
pub struct GCContentFacet {
    pub metrics: GCContentMetrics,
    rng: ThreadRng,
}

impl RecordBasedQualityCheckFacet for GCContentFacet {
    fn name(&self) -> &'static str {
        "GC Content Metrics"
    }

    fn computational_load(&self) -> ComputationalLoad {
        ComputationalLoad::Light
    }

    fn process(&mut self, record: &Record) -> Result<(), super::Error> {
        // (1) Check the record's flags. If any of the flags aren't to our
        // liking, then we reject the record as an ignored flag record.
        let flags = record.flags().unwrap();
        if flags.is_duplicate() || flags.is_secondary() {
            self.metrics.records.ignored_flags += 1;
            return Ok(());
        };

        // (2) Convert the BAM record to a SAM record so we can determine the
        // nucleobases and count up the A's, C's, G's, and T's. TODO: this could
        // be done strictly from the BAM without parsing into SAM.
        let sam_record: sam::record::Sequence = record.sequence().try_into().unwrap();
        let nucleobases = sam_record.as_ref();
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
            self.rng.gen_range(0..max_offset)
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

    fn summarize(&mut self) -> Result<(), Error> {
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

    fn aggregate_results(&self, results: &mut results::Results) {
        results.set_gc_content(self.metrics.clone());
    }
}
