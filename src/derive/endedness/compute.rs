//! Module holding the logic for computing the endedness of a BAM.

use lazy_static::lazy_static;
use noodles::sam::header;
use serde::Serialize;
use std::collections::HashMap;
use std::collections::HashSet;
use std::sync::Arc;
use tracing::warn;

// Strings used to index into the HashMaps used to store the Read Group ordering flags.
// Lazy statics are used to save memory.
lazy_static! {
    /// String used to index into the HashMaps used to store the "overall" ordering flags.
    pub static ref OVERALL: Arc<String> = Arc::new(String::from("overall"));

    /// String used to index into th HashMaps used to store the "unknown_read_group" ordering flags.
    pub static ref UNKNOWN_READ_GROUP: Arc<String> = Arc::new(String::from("unknown_read_group"));
}

/// Struct holding the ordering flags for a single read group.
#[derive(Debug, Clone)]
pub struct OrderingFlagsCounts {
    /// The number of reads with the first in template flag set.
    pub first: usize,

    /// The number of reads with the last in template flag set.
    pub last: usize,

    /// The number of reads with both the first and last in template flags set.
    pub both: usize,

    /// The number of reads with neither the first nor last in template flags set.
    pub neither: usize,
}
impl OrderingFlagsCounts {
    /// Creates a new [`OrderingFlagsCounts`].
    pub fn new() -> Self {
        OrderingFlagsCounts {
            first: 0,
            last: 0,
            both: 0,
            neither: 0,
        }
    }
}

impl Default for OrderingFlagsCounts {
    fn default() -> Self {
        Self::new()
    }
}

/// Struct holding the per read group results for an `ngs derive endedness`
/// subcommand call.
#[derive(Debug, Serialize)]
pub struct ReadGroupDerivedEndednessResult {
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
        counts: OrderingFlagsCounts,
        rpt: Option<f64>,
    ) -> Self {
        ReadGroupDerivedEndednessResult {
            read_group,
            succeeded,
            endedness,
            first: counts.first,
            last: counts.last,
            both: counts.both,
            neither: counts.neither,
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
        counts: OrderingFlagsCounts,
        rpt: Option<f64>,
        read_groups: Vec<ReadGroupDerivedEndednessResult>,
    ) -> Self {
        DerivedEndednessResult {
            succeeded,
            endedness,
            first: counts.first,
            last: counts.last,
            both: counts.both,
            neither: counts.neither,
            rpt,
            read_groups,
        }
    }
}

/// Compares the read group tags found in the records
/// and the read groups found in the header.
/// Returns a vector of read group names that were found in the header
/// but not in the records.
pub fn validate_read_group_info(
    found_rgs: &HashSet<Arc<String>>,
    header: &header::Header,
) -> Vec<String> {
    let mut rgs_in_header_not_records = Vec::new();
    let mut rgs_in_records_not_header = Vec::new();

    for (rg_id, _) in header.read_groups() {
        if !found_rgs.contains(rg_id) {
            rgs_in_header_not_records.push(rg_id.to_string());
        }
    }
    if !rgs_in_header_not_records.is_empty() {
        warn!(
            "The following read groups were not found in the file: {:?}",
            rgs_in_header_not_records
        );
    }

    for rg_id in found_rgs {
        if !header.read_groups().contains_key(rg_id.as_str()) {
            rgs_in_records_not_header.push(rg_id.to_string());
        }
    }
    if !rgs_in_records_not_header.is_empty() {
        warn!(
            "The following read groups were not found in the header: {:?}",
            rgs_in_records_not_header
        );
    }

    rgs_in_header_not_records
}

fn calculate_reads_per_template(
    read_names: HashMap<String, Vec<Arc<String>>>,
) -> HashMap<Arc<String>, f64> {
    let mut reads_per_template: HashMap<Arc<String>, f64> = HashMap::new();
    let mut total_reads: usize = 0;
    let mut total_templates: usize = 0;
    let mut read_group_reads: HashMap<Arc<String>, usize> = HashMap::new();
    let mut read_group_templates: HashMap<Arc<String>, usize> = HashMap::new();

    let mut warning_count: usize = 0;

    for (read_name, read_groups) in read_names.iter() {
        let num_reads = read_groups.len();
        total_reads += num_reads;
        total_templates += 1;

        let read_group_set: HashSet<Arc<String>> = read_groups.iter().cloned().collect();

        if read_group_set.len() == 1 {
            let read_group = Arc::clone(read_group_set.iter().next().unwrap());

            read_group_reads
                .entry(Arc::clone(&read_group))
                .and_modify(|e| *e += num_reads)
                .or_insert(num_reads);
            read_group_templates
                .entry(read_group)
                .and_modify(|e| *e += 1)
                .or_insert(1);
        } else {
            warning_count += 1;
            match warning_count {
                1..=100 => {
                    warn!(
                        "QNAME: '{}' is in multiple read groups: {:?}",
                        read_name, read_group_set
                    );
                }
                101 => warn!(
                    "Too many warnings about QNAMEs in multiple read groups. Stopping warnings."
                ),
                _ => (),
            }

            for read_group in read_groups {
                read_group_reads
                    .entry(Arc::clone(read_group))
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

    if warning_count > 100 {
        warn!(
            "{} QNAMEs were found in multiple read groups.",
            warning_count
        );
    }

    reads_per_template.insert(
        Arc::clone(&OVERALL),
        total_reads as f64 / total_templates as f64,
    );

    for (read_group, num_reads) in read_group_reads.iter() {
        let num_templates = read_group_templates.get(read_group).unwrap();
        let rpt = *num_reads as f64 / *num_templates as f64;
        reads_per_template.insert(Arc::clone(read_group), rpt);
    }

    reads_per_template
}

fn predict_endedness(
    read_group_name: String,
    rg_ordering_flags: &OrderingFlagsCounts,
    paired_deviance: f64,
    reads_per_template: Option<&f64>,
    round_rpt: bool,
) -> ReadGroupDerivedEndednessResult {
    let first = rg_ordering_flags.first;
    let last = rg_ordering_flags.last;
    let both = rg_ordering_flags.both;
    let neither = rg_ordering_flags.neither;

    // all zeroes (Perform this check before creating the result struct
    // so that we don't have to clone the read group name)
    if first == 0 && last == 0 && both == 0 && neither == 0 {
        warn!(
            "No reads were detected in this read group: {}",
            read_group_name
        );
        return ReadGroupDerivedEndednessResult::new(
            read_group_name,
            false,
            "Unknown".to_string(),
            rg_ordering_flags.clone(),
            reads_per_template.copied(),
        );
    }

    let mut result = ReadGroupDerivedEndednessResult::new(
        read_group_name,
        false,
        "Unknown".to_string(),
        rg_ordering_flags.clone(),
        reads_per_template.copied(),
    );

    // only first present
    if first > 0 && last == 0 && both == 0 && neither == 0 {
        return result;
    }
    // only last present
    if first == 0 && last > 0 && both == 0 && neither == 0 {
        return result;
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
        return result;
    }
    // only neither present
    if first == 0 && last == 0 && both == 0 && neither > 0 {
        return result;
    }
    // first/last mixed with both/neither
    if (first > 0 || last > 0) && (both > 0 || neither > 0) {
        return result;
    }
    // any mix of both/neither, regardless of first/last
    if both > 0 && neither > 0 {
        return result;
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
    result
}

/// Main method to evaluate the collected ordering flags and
/// return a result for the endedness of the file. This may fail, and the
/// resulting [`DerivedEndednessResult`] should be evaluated accordingly.
pub fn predict(
    ordering_flags: HashMap<Arc<String>, OrderingFlagsCounts>,
    read_names: HashMap<String, Vec<Arc<String>>>,
    paired_deviance: f64,
    round_rpt: bool,
) -> DerivedEndednessResult {
    let mut rpts: HashMap<Arc<String>, f64> = HashMap::new();
    if !read_names.is_empty() {
        rpts = calculate_reads_per_template(read_names);
    }

    let mut final_result = DerivedEndednessResult::new(
        false,
        "Unknown".to_string(),
        OrderingFlagsCounts::new(),
        None,
        Vec::new(),
    );

    for (read_group, rg_ordering_flags) in ordering_flags.iter() {
        if (*read_group == *UNKNOWN_READ_GROUP)
            && (rg_ordering_flags.first == 0
                && rg_ordering_flags.last == 0
                && rg_ordering_flags.both == 0
                && rg_ordering_flags.neither == 0)
        {
            continue;
        }
        let result = predict_endedness(
            read_group.to_string(),
            rg_ordering_flags,
            paired_deviance,
            rpts.get(read_group),
            round_rpt,
        );
        if result.read_group == "overall" {
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

    final_result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_predict_endedness() {
        let mut ordering_flags: HashMap<Arc<String>, OrderingFlagsCounts> = HashMap::new();
        ordering_flags.insert(
            Arc::clone(&OVERALL),
            OrderingFlagsCounts {
                first: 1,
                last: 1,
                both: 0,
                neither: 0,
            },
        );
        let result = predict_endedness(
            "overall".to_string(),
            ordering_flags.get(&Arc::clone(&OVERALL)).unwrap(),
            0.0,
            None,
            false,
        );
        assert!(result.succeeded);
        assert_eq!(result.endedness, "Paired-End");
        assert_eq!(result.first, 1);
        assert_eq!(result.last, 1);
        assert_eq!(result.both, 0);
        assert_eq!(result.neither, 0);
        assert_eq!(result.rpt, None);
    }

    #[test]
    fn test_derive_endedness_from_all_zero_counts() {
        let mut ordering_flags: HashMap<Arc<String>, OrderingFlagsCounts> = HashMap::new();
        ordering_flags.insert(Arc::new(String::from("rg1")), OrderingFlagsCounts::new());
        let result = predict_endedness(
            String::from("rg1"),
            ordering_flags.get(&Arc::new(String::from("rg1"))).unwrap(),
            0.0,
            None,
            false,
        );
        assert!(!result.succeeded);
        assert_eq!(result.endedness, "Unknown");
        assert_eq!(result.first, 0);
        assert_eq!(result.last, 0);
        assert_eq!(result.both, 0);
        assert_eq!(result.neither, 0);
        assert_eq!(result.rpt, None);
    }

    #[test]
    fn test_derive_endedness_from_only_first() {
        let mut ordering_flags: HashMap<Arc<String>, OrderingFlagsCounts> = HashMap::new();
        ordering_flags.insert(
            Arc::clone(&OVERALL),
            OrderingFlagsCounts {
                first: 1,
                last: 0,
                both: 0,
                neither: 0,
            },
        );
        let result = predict(ordering_flags, HashMap::new(), 0.0, false);
        assert!(!result.succeeded);
        assert_eq!(result.endedness, "Unknown");
        assert_eq!(result.first, 1);
        assert_eq!(result.last, 0);
        assert_eq!(result.both, 0);
        assert_eq!(result.neither, 0);
        assert_eq!(result.rpt, None);
        assert_eq!(result.read_groups.len(), 0);
    }

    #[test]
    fn test_derive_endedness_from_only_last() {
        let mut ordering_flags: HashMap<Arc<String>, OrderingFlagsCounts> = HashMap::new();
        ordering_flags.insert(
            Arc::clone(&OVERALL),
            OrderingFlagsCounts {
                first: 0,
                last: 1,
                both: 0,
                neither: 0,
            },
        );
        let result = predict(ordering_flags, HashMap::new(), 0.0, false);
        assert!(!result.succeeded);
        assert_eq!(result.endedness, "Unknown");
        assert_eq!(result.first, 0);
        assert_eq!(result.last, 1);
        assert_eq!(result.both, 0);
        assert_eq!(result.neither, 0);
        assert_eq!(result.rpt, None);
        assert_eq!(result.read_groups.len(), 0);
    }

    #[test]
    fn test_derive_endedness_from_only_both() {
        let mut ordering_flags: HashMap<Arc<String>, OrderingFlagsCounts> = HashMap::new();
        ordering_flags.insert(
            Arc::clone(&OVERALL),
            OrderingFlagsCounts {
                first: 0,
                last: 0,
                both: 1,
                neither: 0,
            },
        );
        let result = predict(ordering_flags, HashMap::new(), 0.0, false);
        assert!(result.succeeded);
        assert_eq!(result.endedness, "Single-End");
        assert_eq!(result.first, 0);
        assert_eq!(result.last, 0);
        assert_eq!(result.both, 1);
        assert_eq!(result.neither, 0);
        assert_eq!(result.rpt, None);
        assert_eq!(result.read_groups.len(), 0);
    }

    #[test]
    fn test_derive_endedness_from_only_neither() {
        let mut ordering_flags: HashMap<Arc<String>, OrderingFlagsCounts> = HashMap::new();
        ordering_flags.insert(
            Arc::clone(&OVERALL),
            OrderingFlagsCounts {
                first: 0,
                last: 0,
                both: 0,
                neither: 1,
            },
        );
        let result = predict(ordering_flags, HashMap::new(), 0.0, false);
        assert!(!result.succeeded);
        assert_eq!(result.endedness, "Unknown");
        assert_eq!(result.first, 0);
        assert_eq!(result.last, 0);
        assert_eq!(result.both, 0);
        assert_eq!(result.neither, 1);
        assert_eq!(result.rpt, None);
        assert_eq!(result.read_groups.len(), 0);
    }

    #[test]
    fn test_derive_endedness_from_first_and_last() {
        let mut ordering_flags: HashMap<Arc<String>, OrderingFlagsCounts> = HashMap::new();
        ordering_flags.insert(
            Arc::clone(&OVERALL),
            OrderingFlagsCounts {
                first: 1,
                last: 1,
                both: 0,
                neither: 0,
            },
        );
        let result = predict(ordering_flags, HashMap::new(), 0.0, false);
        assert!(result.succeeded);
        assert_eq!(result.endedness, "Paired-End");
        assert_eq!(result.first, 1);
        assert_eq!(result.last, 1);
        assert_eq!(result.both, 0);
        assert_eq!(result.neither, 0);
        assert_eq!(result.rpt, None);
        assert_eq!(result.read_groups.len(), 0);
    }

    #[test]
    fn test_calculate_reads_per_template() {
        let mut read_names: HashMap<String, Vec<Arc<String>>> = HashMap::new();
        let rg_paired = Arc::new("rg_paired".to_string());
        let rg_single = Arc::new("rg_single".to_string());
        read_names.insert(
            "read1".to_string(),
            vec![Arc::clone(&rg_paired), Arc::clone(&rg_paired)],
        );
        read_names.insert(
            "read2".to_string(),
            vec![
                Arc::clone(&rg_paired),
                Arc::clone(&rg_paired),
                Arc::clone(&rg_single),
            ],
        );
        read_names.insert("read3".to_string(), vec![Arc::clone(&rg_single)]);
        read_names.insert(
            "read4".to_string(),
            vec![Arc::clone(&rg_paired), Arc::clone(&rg_paired)],
        );
        read_names.insert(
            "read5".to_string(),
            vec![
                Arc::clone(&rg_paired),
                Arc::clone(&rg_paired),
                Arc::clone(&rg_single),
            ],
        );
        let results = calculate_reads_per_template(read_names);
        assert_eq!(results.len(), 3);
        assert_eq!(results.get(&Arc::new("overall".to_string())).unwrap(), &2.2);
        assert_eq!(results.get(&Arc::clone(&rg_paired)).unwrap(), &2.0);
        assert_eq!(results.get(&Arc::clone(&rg_single)).unwrap(), &1.0);
    }

    #[test]
    fn test_derive_endedness_from_first_and_last_with_rpt() {
        let mut ordering_flags: HashMap<Arc<String>, OrderingFlagsCounts> = HashMap::new();
        let rg_paired = Arc::new("rg_paired".to_string());
        let rg_single = Arc::new("rg_single".to_string());
        ordering_flags.insert(
            Arc::clone(&OVERALL),
            OrderingFlagsCounts {
                first: 8,
                last: 8,
                both: 2,
                neither: 0,
            },
        );
        ordering_flags.insert(
            Arc::clone(&rg_paired),
            OrderingFlagsCounts {
                first: 8,
                last: 8,
                both: 0,
                neither: 0,
            },
        );
        ordering_flags.insert(
            Arc::clone(&rg_single),
            OrderingFlagsCounts {
                first: 0,
                last: 0,
                both: 2,
                neither: 0,
            },
        );
        let mut read_names: HashMap<String, Vec<Arc<String>>> = HashMap::new();
        read_names.insert(
            "read1".to_string(),
            vec![Arc::clone(&rg_paired), Arc::clone(&rg_paired)],
        );
        read_names.insert(
            "read2".to_string(),
            vec![
                Arc::clone(&rg_paired),
                Arc::clone(&rg_paired),
                Arc::clone(&rg_single),
            ],
        );
        read_names.insert("read3".to_string(), vec![Arc::clone(&rg_single)]);
        read_names.insert(
            "read4".to_string(),
            vec![Arc::clone(&rg_paired), Arc::clone(&rg_paired)],
        );
        read_names.insert(
            "read5".to_string(),
            vec![
                Arc::clone(&rg_paired),
                Arc::clone(&rg_paired),
                Arc::clone(&rg_single),
            ],
        );
        let result = predict(ordering_flags, read_names, 0.0, false);
        assert!(!result.succeeded);
        assert_eq!(result.endedness, "Unknown");
        assert_eq!(result.first, 8);
        assert_eq!(result.last, 8);
        assert_eq!(result.both, 2);
        assert_eq!(result.neither, 0);
        assert_eq!(result.rpt, Some(2.2));
        assert_eq!(result.read_groups.len(), 2);
        // We can't know which read group will be first in the vector.
        // But both should succeed.
        assert!(result.read_groups[0].succeeded && result.read_groups[1].succeeded);
    }
}
