//! Functionality related to computing GC content and related metrics.

use noodles_bam::lazy::Record;
use noodles_sam as sam;
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
}

impl GCContentHistogram {
    /// Creates a new GCContentHistogram with default values.
    pub fn new() -> Self {
        Self {
            histogram: SimpleHistogram::zero_based_with_capacity(100),
            total_gc_count: 0,
            total_at_count: 0,
            total_other_count: 0,
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
    pub fn process(&mut self, record: &Record) {
        let flags = record.flags().unwrap();
        if flags.is_duplicate() || flags.is_secondary() || flags.is_unmapped() {
            return;
        };

        // (1) Convert the BAM record to a SAM record so we can determine the
        // nucleotides and count up the A's, C's, G's, and T's. TODO: this could
        // be done strictly from the BAM without parsing into SAM.
        let sam_record: sam::record::Sequence = record.sequence().try_into().unwrap();

        // (2) Initialize necessary counting variables
        let mut gc_this_read = 0usize;

        // (3) Count up the A's, C's, G's, and T's.
        for nucleotide in sam_record.as_ref() {
            match nucleotide {
                Base::C | Base::G => {
                    gc_this_read += 1;
                    self.total_gc_count += 1;
                }
                Base::A | Base::T => self.total_at_count += 1,
                _ => self.total_other_count += 1,
            }
        }

        // (4) Calculate the GC content for this read and increment the
        // histogram accordingly.
        // println!(
        //     "GC: {}, Total: {}",
        //     gc_this_read as f64,
        //     sam_record.as_ref().len() as f64
        // );
        let gc_content_this_read_pct =
            ((gc_this_read as f64 / sam_record.as_ref().len() as f64) * 100.0).floor() as usize;
        // println!("Pct: {}", gc_content_this_read_pct);
        self.histogram.increment(gc_content_this_read_pct).unwrap();
    }
}
