use std::str::FromStr;

use noodles_sam::record::{cigar::op::Kind, sequence::Base, Cigar};

use super::cigar::{consumes_reference, consumes_sequence};

pub struct ReferenceRecordStrepthrough<'a> {
    reference_seq: &'a [u8],
    record_seq: &'a [Base],
    cigar: String,
}

pub fn flatten(cigar: Cigar) -> String {
    let mut result = String::new();

    for op in cigar.iter() {
        let len = op.len();
        let c: char = op.kind().into();
        for _ in 1..=len {
            result.push(c)
        }
    }

    result
}

impl<'a> ReferenceRecordStrepthrough<'a> {
    pub fn new(reference_seq: &'a [u8], record_seq: &'a [Base], cigar: Cigar) -> Self {
        Self {
            reference_seq,
            record_seq,
            cigar: flatten(cigar),
        }
    }

    pub fn edits(&self) -> usize {
        let mut edits = 0;
        let mut record_ptr = 0;
        let mut reference_ptr = 0;

        let chars = self.cigar.chars();

        for c in chars {
            let kind = Kind::from_str(c.to_string().as_str()).unwrap();

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

        edits
    }
}