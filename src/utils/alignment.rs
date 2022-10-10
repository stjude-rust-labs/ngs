//! Utilities related to alignment of sequences.

use anyhow::bail;
use noodles::sam::record::{cigar::op::Kind, sequence::Base, Cigar};

use super::cigar::{consumes_reference, consumes_sequence};

/// Turns a condensed Cigar representation into a flattened representation. For
/// example, 10M will become a Vec of length 10 comprised completely of
/// Kind::MATCH. This utility is useful for generating a representation of a
/// Cigar string that we can step through alongside a sequence.
pub fn flatten(cigar: Cigar) -> Vec<Kind> {
    let mut result = Vec::new();

    for op in cigar.iter() {
        let len = op.len();
        let kind = op.kind();
        for _ in 1..=len {
            result.push(kind)
        }
    }

    result
}

/// Utility struct for stepping through a reference sequence, and record
/// sequence, and a Cigar string in unison.
pub struct ReferenceRecordStepThrough<'a> {
    reference_seq: &'a [u8],
    record_seq: &'a [Base],
    cigar: Vec<Kind>,
}

impl<'a> ReferenceRecordStepThrough<'a> {
    /// Creates a new [`ReferenceRecordStepThrough`].
    pub fn new(reference_seq: &'a [u8], record_seq: &'a [Base], cigar: Cigar) -> Self {
        Self {
            reference_seq,
            record_seq,
            cigar: flatten(cigar),
        }
    }

    /// Calculates the number of edits in the [`ReferenceRecordStepThrough`] by
    /// stepping through the genome and counting up all of the mismatched M's.
    /// Errors can occur if the reference or the sequence are not all the way
    /// consumed.
    pub fn edits(&self) -> anyhow::Result<usize> {
        let mut edits = 0;
        let mut record_ptr = 0;
        let mut reference_ptr = 0;

        for kind in self.cigar.iter().copied() {
            if kind == Kind::Match {
                let ref_base = self.reference_seq[reference_ptr] as char;
                let record_base: char = self.record_seq[record_ptr].into();
                if ref_base != record_base {
                    edits += 1;
                }
            }

            if consumes_reference(kind) {
                reference_ptr += 1;
            }

            if consumes_sequence(kind) {
                record_ptr += 1;
            }
        }

        if self.reference_seq.len() != reference_ptr {
            bail!("reference sequence was not fully consumed");
        } else if self.record_seq.len() != record_ptr {
            bail!("record sequence was not fully consumed");
        }

        Ok(edits)
    }
}
