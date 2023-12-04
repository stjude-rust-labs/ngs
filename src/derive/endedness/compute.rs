//! Module holding the logic for computing the endedness of a BAM.

use anyhow::bail;
use serde::Serialize;
use std::collections::HashMap;

/// Struct holding the final results for an `ngs derive endedness` subcommand
/// call.
#[derive(Debug, Serialize)]
pub struct DerivedEndednessResult {
    /// Whether or not the `ngs derive endedness` subcommand succeeded.
    pub succeeded: bool,

    /// The endedness, if available.
    pub endedness: String,

    /// The f+l- read count.
    pub first: usize,

    /// The f-l+ read count.
    pub last: usize,

    /// The f+l+ read count.
    pub both: usize,

    /// The f-l- read count.
    pub neither: usize,
}

impl DerivedEndednessResult {
    /// Creates a new [`DerivedEndednessResult`].
    pub fn new(
        succeeded: bool,
        endedness: String,
        first: usize,
        last: usize,
        both: usize,
        neither: usize,
    ) -> Self {
        DerivedEndednessResult {
            succeeded,
            endedness,
            first,
            last,
            both,
            neither,
        }
    }
}

/// Main method to evaluate the collected ordering flags and
/// return a result for the endedness of the file. This may fail, and the
/// resulting [`DerivedEndednessResult`] should be evaluated accordingly.
pub fn predict(
    ordering_flags: HashMap<String, HashMap<String, usize>>,
    paired_deviance: f64,
) -> Result<DerivedEndednessResult, anyhow::Error> {
    let first = ordering_flags["overall"]["f+l-"];
    let last = ordering_flags["overall"]["f-l+"];
    let both = ordering_flags["overall"]["f+l+"];
    let neither = ordering_flags["overall"]["f-l-"];

    let mut result =
        DerivedEndednessResult::new(false, "Unknown".to_string(), first, last, both, neither);

    // all zeroes
    if first == 0 && last == 0 && both == 0 && neither == 0 {
        bail!("No reads were detected in the file.");
    }

    // only first present
    if first > 0 && last == 0 && both == 0 && neither == 0 {
        return Ok(result);
    }
    // only last present
    if first == 0 && last > 0 && both == 0 && neither == 0 {
        return Ok(result);
    }
    // only both present
    if first == 0 && last == 0 && both > 0 && neither == 0 {
        result.succeeded = true;
        result.endedness = "Single-End".to_string();
        return Ok(result);
    }
    // only neither present
    if first == 0 && last == 0 && both == 0 && neither > 0 {
        return Ok(result);
    }
    // first/last mixed with both/neither
    if (first > 0 || last > 0) && (both > 0 || neither > 0) {
        return Ok(result);
    }
    // any mix of both/neither, regardless of first/last
    if both > 0 && neither > 0 {
        return Ok(result);
    }

    // both and neither are now guarenteed to be 0
    // We only need to check first and last

    let first_frac = first as f64 / (first + last) as f64;
    let lower_limit = 0.5 - paired_deviance;
    let upper_limit = 0.5 + paired_deviance;
    if (first == last) || (lower_limit <= first_frac && first_frac <= upper_limit) {
        result.succeeded = true;
        result.endedness = "Paired-End".to_string();
        return Ok(result);
    }

    Ok(result)
}
