//! Utilities related to CIGAR string processing.

use noodles::sam::record::cigar::op::Kind;

/// Reports whether a CIGAR operation consumes a reference base.
pub fn consumes_reference(kind: Kind) -> bool {
    matches!(
        kind,
        Kind::Match | Kind::Deletion | Kind::Skip | Kind::SequenceMatch | Kind::SequenceMismatch
    )
}

/// Reports whether a CIGAR operation consumes a sequence base.
pub fn consumes_sequence(kind: Kind) -> bool {
    matches!(
        kind,
        Kind::Match
            | Kind::Insertion
            | Kind::SoftClip
            | Kind::SequenceMatch
            | Kind::SequenceMismatch
    )
}
