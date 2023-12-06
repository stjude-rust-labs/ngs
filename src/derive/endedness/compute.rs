//! Module holding the logic for computing the endedness of a BAM.

use anyhow::bail;
use lazy_static::lazy_static;
use radix_trie::Trie;
use radix_trie::TrieCommon;
use serde::Serialize;
use std::collections::HashMap;
use std::collections::HashSet;
use std::rc::Rc;
use tracing::warn;

lazy_static! {
    // Strings used to index into the HashMaps used to store the Read Group ordering flags.
    pub static ref OVERALL: String = String::from("overall");
    pub static ref UNKNOWN_READ_GROUP: String = String::from("unknown_read_group");
}

pub struct OrderingFlagsCounts {
    pub first: usize,
    pub last: usize,
    pub both: usize,
    pub neither: usize,
}

impl OrderingFlagsCounts {
    pub fn new() -> Self {
        OrderingFlagsCounts {
            first: 0,
            last: 0,
            both: 0,
            neither: 0,
        }
    }
}

/// Struct holding the per read group results for an `ngs derive endedness`
/// subcommand call.
#[derive(Debug, Serialize)]
struct ReadGroupDerivedEndednessResult {
    /// Name of the read group.
    pub read_group: String,

    /// Whether or not an endedness was determined for this read group.
    pub succeeded: bool,

    /// The endedness of this read group or "Unknown".
    pub endedness: String,

    /// The f+l- read count.
    pub first: usize,

    /// The f-l+ read count.
    pub last: usize,

    /// The f+l+ read count.
    pub both: usize,

    /// The f-l- read count.
    pub neither: usize,

    /// The reads per template (RPT).
    /// Only available if `args.calc_rpt` is true.
    pub rpt: Option<f64>,
}

impl ReadGroupDerivedEndednessResult {
    /// Creates a new [`ReadGroupDerivedEndednessResult`].
    fn new(
        read_group: String,
        succeeded: bool,
        endedness: String,
        first: usize,
        last: usize,
        both: usize,
        neither: usize,
        rpt: Option<f64>,
    ) -> Self {
        ReadGroupDerivedEndednessResult {
            read_group,
            succeeded,
            endedness,
            first,
            last,
            both,
            neither,
            rpt,
        }
    }
}

/// Struct holding the final results for an `ngs derive endedness` subcommand
/// call.
#[derive(Debug, Serialize)]
pub struct DerivedEndednessResult {
    /// Whether or not the `ngs derive endedness` subcommand succeeded.
    pub succeeded: bool,

    /// The overall endedness of the file or "Unknown".
    pub endedness: String,

    /// The overall f+l- read count.
    pub first: usize,

    /// The overall f-l+ read count.
    pub last: usize,

    /// The overall f+l+ read count.
    pub both: usize,

    /// The overall f-l- read count.
    pub neither: usize,

    /// The overall reads per template (RPT).
    /// Only available if `args.calc_rpt` is true.
    pub rpt: Option<f64>,

    /// Vector of [`ReadGroupDerivedEndednessResult`]s.
    /// One for each read group in the BAM,
    /// and potentially one for any reads with an unknown read group.
    pub read_groups: Vec<ReadGroupDerivedEndednessResult>,
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
        rpt: Option<f64>,
        read_groups: Vec<ReadGroupDerivedEndednessResult>,
    ) -> Self {
        DerivedEndednessResult {
            succeeded,
            endedness,
            first,
            last,
            both,
            neither,
            rpt,
            read_groups,
        }
    }
}

fn calculate_reads_per_template(read_names: Trie<String, Vec<Rc<&str>>>) -> HashMap<Rc<&str>, f64> {
    let mut reads_per_template: HashMap<Rc<&str>, f64> = HashMap::new();
    let mut total_reads: usize = 0;
    let mut total_templates: usize = 0;
    let mut read_group_reads: HashMap<Rc<&str>, usize> = HashMap::new();
    let mut read_group_templates: HashMap<Rc<&str>, usize> = HashMap::new();

    for (read_name, read_groups) in read_names.iter() {
        let num_reads = read_groups.len();
        total_reads += num_reads;
        total_templates += 1;

        let read_group_set: HashSet<Rc<&str>> = read_groups.iter().cloned().collect();

        if read_group_set.len() == 1 {
            let read_group = read_group_set.iter().next().unwrap();

            read_group_reads
                .entry(read_group.clone())
                .and_modify(|e| *e += num_reads)
                .or_insert(num_reads);
            read_group_templates
                .entry(read_group.clone())
                .and_modify(|e| *e += 1)
                .or_insert(1);
        } else {
            warn!(
                "Read {} has multiple read groups: {:#?}",
                read_name, read_groups
            );
            for read_group in read_groups {
                let read_group = Rc::new(**read_group);
                read_group_reads
                    .entry(read_group)
                    .and_modify(|e| *e += 1)
                    .or_insert(1);
            }
            for read_group in read_group_set {
                read_group_templates
                    .entry(read_group)
                    .and_modify(|e| *e += 1)
                    .or_insert(1);
            }
        }
    }

    reads_per_template.insert(
        Rc::new(OVERALL.as_str()),
        total_reads as f64 / total_templates as f64,
    );

    for (read_group, num_reads) in read_group_reads.iter() {
        let num_templates = read_group_templates.get(read_group).unwrap();
        let rpt = *num_reads as f64 / *num_templates as f64;
        reads_per_template.insert(Rc::clone(read_group), rpt);
    }

    reads_per_template
}

fn predict_endedness(
    read_group_name: String,
    rg_ordering_flags: &OrderingFlagsCounts,
    paired_deviance: f64,
    reads_per_template: Option<&f64>,
    round_rpt: bool,
) -> Result<ReadGroupDerivedEndednessResult, anyhow::Error> {
    let first = rg_ordering_flags.first;
    let last = rg_ordering_flags.last;
    let both = rg_ordering_flags.both;
    let neither = rg_ordering_flags.neither;

    let mut result = ReadGroupDerivedEndednessResult::new(
        read_group_name,
        false,
        "Unknown".to_string(),
        first,
        last,
        both,
        neither,
        reads_per_template.copied(),
    );

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
        match reads_per_template {
            Some(rpt) => {
                if *rpt == 1.0 || (round_rpt && rpt.round() as usize == 1) {
                    result.succeeded = true;
                    result.endedness = String::from("Single-End");
                }
            }
            None => {
                result.succeeded = true;
                result.endedness = String::from("Single-End");
            }
        }
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
        match reads_per_template {
            Some(rpt) => {
                if *rpt == 2.0 || (round_rpt && rpt.round() as usize == 2) {
                    result.succeeded = true;
                    result.endedness = String::from("Paired-End");
                }
            }
            None => {
                result.succeeded = true;
                result.endedness = String::from("Paired-End");
            }
        }
    }

    return Ok(result);
}

/// Main method to evaluate the collected ordering flags and
/// return a result for the endedness of the file. This may fail, and the
/// resulting [`DerivedEndednessResult`] should be evaluated accordingly.
pub fn predict(
    ordering_flags: HashMap<Rc<&str>, OrderingFlagsCounts>,
    read_names: Trie<String, Vec<Rc<&str>>>,
    paired_deviance: f64,
    round_rpt: bool,
) -> Result<DerivedEndednessResult, anyhow::Error> {
    let mut rpts: HashMap<Rc<&str>, f64> = HashMap::new();
    if !read_names.is_empty() {
        rpts = calculate_reads_per_template(read_names);
    }

    let mut final_result =
        DerivedEndednessResult::new(false, "Unknown".to_string(), 0, 0, 0, 0, None, Vec::new());

    for (read_group, rg_ordering_flags) in ordering_flags.iter() {
        let result = predict_endedness(
            String::from(**read_group),
            rg_ordering_flags,
            paired_deviance,
            rpts.get(read_group),
            round_rpt,
        )?;
        if result.read_group == String::from("overall") {
            final_result.endedness = result.endedness;
            final_result.first = result.first;
            final_result.last = result.last;
            final_result.both = result.both;
            final_result.neither = result.neither;
            final_result.rpt = result.rpt;
            final_result.succeeded = result.succeeded;
        } else {
            final_result.read_groups.push(result);
        }
    }

    Ok(final_result)
}
