//! Module holding the logic for computing the quality score encoding.

use anyhow::bail;
use serde::Serialize;
use std::collections::HashSet;

const MAX_VALID_PHRED_SCORE: u8 = 93;
const SANGER_MIN: u8 = 0;
const ILLUMINA_1_0_MIN: u8 = 26;
const ILLUMINA_1_3_MIN: u8 = 31;

/// Struct holding the final results for an `ngs derive encoding` subcommand
/// call.
#[derive(Debug, Serialize)]
pub struct DerivedEncodingResult {
    /// Whether or not the `ngs derive encoding` subcommand succeeded.
    pub succeeded: bool,

    /// The detected quality score encoding, if available.
    pub encoding: Option<String>,

    /// The minimum quality score observed.
    pub observed_min: u8,

    /// The maximum quality score observed.
    pub observed_max: u8,
}

impl DerivedEncodingResult {
    /// Creates a new [`DerivedEncodingResult`].
    pub fn new(
        succeeded: bool,
        encoding: Option<String>,
        observed_min: u8,
        observed_max: u8,
    ) -> Self {
        DerivedEncodingResult {
            succeeded,
            encoding,
            observed_min,
            observed_max,
        }
    }
}

/// Main method to evaluate the observed quality scores and
/// return a result for the derived encoding. This may fail, and the
/// resulting [`DerivedEncodingResult`] should be evaluated accordingly.
pub fn predict(score_set: HashSet<u8>) -> Result<DerivedEncodingResult, anyhow::Error> {
    if score_set.is_empty() {
        bail!("No quality scores were detected in the file.");
    }

    let observed_min = *score_set.iter().min().unwrap();
    let observed_max = *score_set.iter().max().unwrap();

    let mut result = DerivedEncodingResult::new(false, None, observed_min, observed_max);

    if observed_max > MAX_VALID_PHRED_SCORE {
        return anyhow::Ok(result);
    }
    match observed_min {
        ILLUMINA_1_3_MIN..=MAX_VALID_PHRED_SCORE => {
            result.succeeded = true;
            result.encoding = Some("Illumina 1.3".to_string());
        }
        ILLUMINA_1_0_MIN..=MAX_VALID_PHRED_SCORE => {
            result.succeeded = true;
            result.encoding = Some("Illumina 1.0".to_string());
        }
        SANGER_MIN..=MAX_VALID_PHRED_SCORE => {
            result.succeeded = true;
            result.encoding = Some("Sanger/Illumina 1.8".to_string());
        }
        _ => bail!("This shouldn't be possible!"),
    }

    anyhow::Ok(result)
}
