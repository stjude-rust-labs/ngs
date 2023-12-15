//! Module holding the logic for computing the consensus read length.

use anyhow::bail;
use serde::Serialize;
use std::collections::HashMap;

/// Struct holding the final results for an `ngs derive readlen` subcommand
/// call.
#[derive(Debug, Serialize)]
pub struct DerivedReadlenResult {
    /// Whether or not the `ngs derive readlen` subcommand succeeded.
    pub succeeded: bool,

    /// The concsensus read length, if available.
    pub consensus_read_length: Option<usize>,

    /// The majority vote percentage of the consensus read length, if available.
    pub majority_pct_detected: f64,

    /// Status of the evidence that supports (or does not support) this
    /// read length, if available.
    pub evidence: Vec<(usize, usize)>,
}

impl DerivedReadlenResult {
    /// Creates a new [`DerivedReadlenResult`].
    pub fn new(
        succeeded: bool,
        consensus_read_length: Option<usize>,
        majority_pct_detected: f64,
        evidence: Vec<(usize, usize)>,
    ) -> Self {
        DerivedReadlenResult {
            succeeded,
            consensus_read_length,
            majority_pct_detected,
            evidence,
        }
    }
}

/// Main method to evaluate the collected read lengths and
/// return a result for the consensus read length. This may fail, and the
/// resulting [`DerivedReadlenResult`] should be evaluated accordingly.
pub fn predict(
    read_lengths: HashMap<usize, usize>,
    num_samples: usize,
    majority_vote_cutoff: f64,
) -> Result<DerivedReadlenResult, anyhow::Error> {
    if num_samples <= 0 {
        bail!("No read lengths were detected in the file.");
    }

    // Sort the read lengths by their key for output.
    let mut read_lengths: Vec<(usize, usize)> = read_lengths.into_iter().collect();
    read_lengths.sort_by(|a, b| b.0.cmp(&a.0));

    let max_read_length = read_lengths[0].0;
    let max_count = read_lengths[0].1;

    let consensus_read_length = max_read_length;
    let majority_detected = max_count as f64 / num_samples as f64;

    let mut result =
        DerivedReadlenResult::new(false, None, majority_detected * 100.0, read_lengths);

    if majority_detected >= majority_vote_cutoff {
        result.succeeded = true;
        result.consensus_read_length = Some(consensus_read_length);
    }

    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_derive_readlen_from_empty_hashmap() {
        let read_lengths = HashMap::new();
        let result = predict(read_lengths, 0, 0.7);
        assert!(result.is_err());
    }

    #[test]
    fn test_derive_readlen_when_all_readlengths_equal() {
        let read_lengths = HashMap::from([(100, 10)]);
        let result = predict(read_lengths, 10, 1.0).unwrap();
        assert!(result.succeeded);
        assert_eq!(result.consensus_read_length, Some(100));
        assert_eq!(result.majority_pct_detected, 100.0);
        assert_eq!(result.evidence, Vec::from([(100, 10)]));
    }

    #[test]
    fn test_derive_readlen_success_when_not_all_readlengths_equal() {
        let read_lengths = HashMap::from([(101, 1000), (100, 5), (99, 5)]);
        let result = predict(read_lengths, 1010, 0.7).unwrap();
        assert!(result.succeeded);
        assert_eq!(result.consensus_read_length, Some(101));
        assert!(result.majority_pct_detected > 99.0);
        assert_eq!(result.evidence, Vec::from([(101, 1000), (100, 5), (99, 5)]));
    }

    #[test]
    fn test_derive_readlen_fail_when_not_all_readlengths_equal() {
        let read_lengths = HashMap::from([(101, 5), (100, 1000), (99, 5)]);
        let result = predict(read_lengths, 1010, 0.7).unwrap();
        assert!(!result.succeeded);
        assert_eq!(result.consensus_read_length, None);
        assert!(result.majority_pct_detected < 0.7);
        assert_eq!(result.evidence, Vec::from([(101, 5), (100, 1000), (99, 5)]));
    }
}
