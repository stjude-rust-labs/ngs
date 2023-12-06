//! Module holding the logic for computing the endedness of a BAM.

use anyhow::bail;
use lazy_static::lazy_static;
use radix_trie::iter;
use radix_trie::Trie;
use serde::Serialize;
use std::collections::HashMap;
use std::rc::Rc;

lazy_static! {
    // Strings used to index into the HashMaps used to store the Read Group ordering flags.
    pub static ref OVERALL: String = String::from("overall");
    pub static ref UNKNOWN_READ_GROUP: String = String::from("unknown_read_group");

    // Strings used to index into the HashMaps used to store the ordering flag counts.
    pub static ref FIRST: String = String::from("f+l-");
    pub static ref LAST: String = String::from("f-l+");
    pub static ref BOTH: String = String::from("f+l+");
    pub static ref NEITHER: String = String::from("f-l-");
}

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

fn calculate_reads_per_template(
    read_names: Trie<String, Vec<Rc<String>>>,
) -> HashMap<Rc<String>, f64> {
    let mut total_reads: usize = 0;
    let mut templates_templates: usize = 0;
    let mut reads_per_template: HashMap<Rc<String>, f64> = HashMap::new();

    for read_groups in iter::Iter::new(read_names) {}

    reads_per_template
}

/// Main method to evaluate the collected ordering flags and
/// return a result for the endedness of the file. This may fail, and the
/// resulting [`DerivedEndednessResult`] should be evaluated accordingly.
pub fn predict(
    ordering_flags: HashMap<Rc<String>, HashMap<String, usize>>,
    read_names: Trie<String, Vec<Rc<String>>>,
    paired_deviance: f64,
) -> Result<DerivedEndednessResult, anyhow::Error> {
    let overall_ordering_flags = ordering_flags.get(&*OVERALL).unwrap();

    let overall_first = *overall_ordering_flags.get(&*FIRST).unwrap();
    let overall_last = *overall_ordering_flags.get(&*LAST).unwrap();
    let overall_both = *overall_ordering_flags.get(&*BOTH).unwrap();
    let overall_neither = *overall_ordering_flags.get(&*NEITHER).unwrap();

    let mut result = DerivedEndednessResult::new(
        false,
        "Unknown".to_string(),
        overall_first,
        overall_last,
        overall_both,
        overall_neither,
    );

    // all zeroes
    if overall_first == 0 && overall_last == 0 && overall_both == 0 && overall_neither == 0 {
        bail!("No reads were detected in the file.");
    }

    // only first present
    if overall_first > 0 && overall_last == 0 && overall_both == 0 && overall_neither == 0 {
        return Ok(result);
    }
    // only last present
    if overall_first == 0 && overall_last > 0 && overall_both == 0 && overall_neither == 0 {
        return Ok(result);
    }
    // only both present
    if overall_first == 0 && overall_last == 0 && overall_both > 0 && overall_neither == 0 {
        result.succeeded = true;
        result.endedness = "Single-End".to_string();
        return Ok(result);
    }
    // only neither present
    if overall_first == 0 && overall_last == 0 && overall_both == 0 && overall_neither > 0 {
        return Ok(result);
    }
    // first/last mixed with both/neither
    if (overall_first > 0 || overall_last > 0) && (overall_both > 0 || overall_neither > 0) {
        return Ok(result);
    }
    // any mix of both/neither, regardless of first/last
    if overall_both > 0 && overall_neither > 0 {
        return Ok(result);
    }

    // both and neither are now guarenteed to be 0
    // We only need to check first and last

    let first_frac = overall_first as f64 / (overall_first + overall_last) as f64;
    let lower_limit = 0.5 - paired_deviance;
    let upper_limit = 0.5 + paired_deviance;
    if (overall_first == overall_last) || (lower_limit <= first_frac && first_frac <= upper_limit) {
        result.succeeded = true;
        result.endedness = "Paired-End".to_string();
        return Ok(result);
    }

    Ok(result)
}
