//! Functionality related to computing GC content and related metrics.

use noodles_bam::lazy::Record;
use noodles_sam as sam;
use rand::prelude::*;
use sam::record::sequence::Base;
use serde::Serialize;

use crate::utils::histogram::SimpleHistogram;

/// Primary struct used to compile stats regarding GC content. Within this
/// struct, the histogram represents the number of reads which have 0% GC
/// content all the way up to 100% GC content. The other fields are for counting
/// the number of nucleotides which are G/C, the number of nucleotides that are
/// A/T, or the number of nucleotides that fall into other.
#[derive(Debug, Serialize)]
pub struct GCContentHistogram {
    // Histogram that represents the number of reads which have 0% GC content
    // all the way up to 100% GC content.
    histogram: SimpleHistogram,

    // Total number of nucleotides read as a 'G' or a 'C'.
    total_gc_count: usize,

    // Total number of nucleotides read as an 'A' or a 'T'.
    total_at_count: usize,

    // Total number of nucleotides which were not read as an 'A', a 'C', a 'G',
    // or a 'T'.
    total_other_count: usize,

    // Number of records that have been processed by this struct.
    processed: usize,

    // Number of records that were ignored because of their flags.
    ignored_flags: usize,

    // Number of records that were ignored because they were too short.
    ignored_too_short: usize,
}

impl GCContentHistogram {
    /// Creates a new GCContentHistogram with default values.
    pub fn new() -> Self {
        Self {
            histogram: SimpleHistogram::zero_based_with_capacity(100),
            total_gc_count: 0,
            total_at_count: 0,
            total_other_count: 0,
            processed: 0,
            ignored_flags: 0,
            ignored_too_short: 0,
        }
    }

    /// Gets the number of G/C nucleotides read so far by this struct.
    pub fn get_gc_count(&self) -> usize {
        self.total_gc_count
    }

    /// Gets the number of A/T nucleotides read so far by this struct.
    pub fn get_at_count(&self) -> usize {
        self.total_at_count
    }

    /// Gets the number of non-ACGT nucleotides read so far by this struct.
    pub fn get_other_count(&self) -> usize {
        self.total_other_count
    }

    /// Processes a single read and updates the struct accordingly.
    pub fn process(&mut self, record: &Record, rng: &mut ThreadRng) {
        const TRUNCATION_LENGTH: usize = 100;

        // (1) Check the record's flags. If any of the flags aren't to our
        // liking, then we reject the record as an ignored flag record.
        let flags = record.flags().unwrap();
        if flags.is_duplicate() || flags.is_secondary() || flags.is_unmapped() {
            self.ignored_flags += 1;
            return;
        };

        // (2) Convert the BAM record to a SAM record so we can determine the
        // nucleotides and count up the A's, C's, G's, and T's. TODO: this could
        // be done strictly from the BAM without parsing into SAM.
        let sam_record: sam::record::Sequence = record.sequence().try_into().unwrap();
        let nucleotides = sam_record.as_ref();
        let sequence_length = nucleotides.len();

        // (3) Checks whether the record is too short to check the GC bias for.
        // If the record is too short, it can only be a subset of the full range
        // of percentages from 0%-100%. We don't want that messing with our
        // results—all records should have a fair chance to generate between 0%
        // to 100% GC content.
        if sequence_length < TRUNCATION_LENGTH {
            self.ignored_too_short += 1;
            return;
        }

        // (4) Truncate the read TRUNCATION_LENGTH (generally 100 nucleotides)
        // so that we have a uniform chance of generating between 0%-100%. We
        // choose a random starting point within the record so as to not
        // introduce a bias towards the start of our record.
        let mut gc_this_read = 0usize;
        let offset = if TRUNCATION_LENGTH < sequence_length {
            let max_offset = sequence_length - TRUNCATION_LENGTH;
            rng.gen_range(0..max_offset)
        } else {
            0
        };

        // (5) Count up the A's, C's, G's, and T's.
        for i in 0..TRUNCATION_LENGTH {
            let nucleotide = nucleotides[offset + i];
            match nucleotide {
                Base::C | Base::G => {
                    gc_this_read += 1;
                    self.total_gc_count += 1;
                }
                Base::A | Base::T => self.total_at_count += 1,
                _ => self.total_other_count += 1,
            }
        }

        // (6) Calculate the GC content for this read and increment the
        // histogram accordingly.
        let gc_content_this_read_pct =
            ((gc_this_read as f64 / TRUNCATION_LENGTH as f64) * 100.0).round() as usize;
        self.histogram.increment(gc_content_this_read_pct).unwrap();
        self.processed += 1;
    }
}
