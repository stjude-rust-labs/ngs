use noodles_sam::record::cigar::op::Kind;

pub fn consumes_reference(kind: Kind) -> bool {
    matches!(
        kind,
        Kind::Match | Kind::Deletion | Kind::Skip | Kind::SequenceMatch | Kind::SequenceMismatch
    )
}

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
