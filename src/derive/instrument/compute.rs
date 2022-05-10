use std::collections::{HashMap, HashSet};

use regex::Regex;
use serde::Serialize;
use tracing::info;

#[derive(Debug, Serialize)]
pub struct InstrumentDetectionResults {
    pub possible_instruments: HashSet<String>,
    pub detected_at_least_one_machine: bool,
    pub initialized: bool,
}

impl InstrumentDetectionResults {
    pub fn new() -> Self {
        InstrumentDetectionResults {
            possible_instruments: HashSet::new(),
            detected_at_least_one_machine: false,
            initialized: false,
        }
    }

    pub fn update_instruments(&mut self, results: &HashSet<String>) {
        if self.initialized {
            // An initial base set has already been established, so take the
            // intersection of the existing possible instruments set and the
            // set being passed into the update function.
            let new_set = self
                .possible_instruments
                .intersection(&results)
                .map(|s| s.clone())
                .collect();
            self.possible_instruments = new_set;
        } else {
            // This is the first iteration, so we need to set our base set as
            // the first detected set of results.
            self.possible_instruments = results.clone();
            self.initialized = true;
        }

        // After we've updated the sets, we need to keep track of if we have
        // detected at least one machine to distinguish conflicting machines
        // from no machines detected at all.
        if !results.is_empty() {
            self.detected_at_least_one_machine = true;
        }
    }
}

#[derive(Debug, Serialize)]
pub struct DerivedInstrumentResult {
    pub succeeded: bool,
    pub instruments: Option<HashSet<String>>,
    pub confidence: String,
    pub evidence: Option<String>,
    pub comment: Option<String>,
}

impl DerivedInstrumentResult {
    pub fn from(
        succeeded: bool,
        instruments: Option<HashSet<String>>,
        confidence: String,
        evidence: Option<String>,
        comment: Option<String>,
    ) -> Self {
        DerivedInstrumentResult {
            succeeded,
            instruments,
            confidence,
            evidence,
            comment,
        }
    }
}

/// Computes the full set of possible instruments that could have generated the
/// value passed to the function. The `lookup_table` passed to this function is
/// constructed as a set of regex keys that map to machine that could have
/// generated a query that matches that regex. Effectively, this method iterates
/// through all of the keys, checks if the query matches the regex, and extends
/// the result HashSet with the values for that key if it does.
///
/// # Arguments
///
/// * `query` — A value to check against in the lookup HashMap.
/// * `instruments` — Lookup HashMap where each key represents a regex that
///   matches to the possible machines that generated the query contained within
///   the respective HashSet.
pub fn possible_instruments_for_query(
    query: String,
    lookup_table: &HashMap<&'static str, HashSet<&'static str>>,
) -> HashSet<String> {
    let mut result: HashSet<String> = HashSet::new();

    for (pattern, machines) in lookup_table {
        let re = Regex::new(pattern).unwrap();
        if re.is_match(query.as_str()) {
            let matching_machines: Vec<String> =
                machines.into_iter().map(|x| x.to_string()).collect();
            result.extend(matching_machines);
        }
    }

    info!(" [*] {}, Possible Instruments: {:?}", query, result);
    result
}

/// Given a HashSet of unique queries (usually a instrument ID or flowcell ID
/// parsed from a read name) that were detected from a SAM/BAM/CRAM file, return
/// a HashSet that contains all possible machines that could have generated that
/// list of queries.
///
/// This is done by iterating through the HashSet of machines that could have
/// produced each name and taking the intersection. It is possible, of course,
/// that there are multiple machines that generated the data contained within a
/// single file. In these cases, the result of this function would be an empty
/// HashSet, with no distinguishing factors between a failed lookup (a lookup
/// for which the query matched none of the regex keys) and conflicting
/// instrument names. Thus, we create a special return object here that also
/// allows for this special case to be flagged.
///
/// # Arguments
///
/// * `queries` — All of the queries detected in the SAM/BAM/CRAM files.
/// * `lookup_table` — Lookup HashMap where each key represents a regex that
///   matches to the possible machines that generated the name contained within
///   the respective HashSet.
pub fn predict_instrument(
    queries: HashSet<String>,
    lookup_table: &HashMap<&'static str, HashSet<&'static str>>,
) -> InstrumentDetectionResults {
    let mut result = InstrumentDetectionResults::new();

    for name in queries {
        let derived = possible_instruments_for_query(name, lookup_table);
        result.update_instruments(&derived);
    }

    result
}

pub fn resolve_instrument_prediction(
    iid_results: InstrumentDetectionResults,
    fcid_results: InstrumentDetectionResults,
) -> DerivedInstrumentResult {
    if iid_results.possible_instruments.is_empty() && iid_results.detected_at_least_one_machine {
        return DerivedInstrumentResult::from(
            false,
            None,
            "unknown".to_string(),
            Some("instrument id".to_string()),
            Some(
                "multiple instruments were detected in this file via the instrument id".to_string(),
            ),
        );
    }

    if fcid_results.possible_instruments.is_empty() && fcid_results.detected_at_least_one_machine {
        return DerivedInstrumentResult::from(
            false,
            None,
            "unknown".to_string(),
            Some("flowcell id".to_string()),
            Some("multiple instruments were detected in this file via the flowcell id".to_string()),
        );
    }

    // (1) if neither results turn up anything, then evaluating the results is
    // relatively easy: you just check to see if any machine was detected at
    // all, and that will let you know if there were multiple machines in the
    // run or not.
    if iid_results.possible_instruments.is_empty() && fcid_results.possible_instruments.is_empty() {
        return DerivedInstrumentResult::from(
            false,
            None,
            "unknown".to_string(),
            None,
            Some("no matching instruments were found".to_string()),
        );
    }

    // (2) So, if both aren't None, are iid_results _is_ none, then the fcid
    // results must not be empty. We can go ahead and issue a prediction based
    // on this medium to low confidence result.
    if iid_results.possible_instruments.is_empty() {
        let instruments = fcid_results.possible_instruments;
        let mut confidence = "medium";
        if instruments.len() > 1 {
            confidence = "low";
        }

        return DerivedInstrumentResult::from(
            true,
            Some(instruments),
            confidence.to_string(),
            Some("flowcell id".to_string()),
            None,
        );
    }

    // (3) Same as the block above, except now we are evaluating the opposite
    // (only the iid results contained some predicted machine).
    if fcid_results.possible_instruments.is_empty() {
        let instruments = iid_results.possible_instruments;
        let mut confidence = "medium";
        if instruments.len() > 1 {
            confidence = "low";
        }

        return DerivedInstrumentResult::from(
            true,
            Some(instruments),
            confidence.to_string(),
            Some("instrument id".to_string()),
            None,
        );
    }

    let overlapping_instruments: HashSet<String> = fcid_results
        .possible_instruments
        .intersection(&iid_results.possible_instruments)
        .map(|r| r.clone())
        .collect();

    if overlapping_instruments.is_empty() {
        return DerivedInstrumentResult::from(
            false,
            None,
            "high".to_string(),
            Some("instrument and flowcell id".to_string()),
            Some(
                "Case needs triaging, results from instrument id and \
                         flowcell id are mutually exclusive."
                    .to_string(),
            ),
        );
    }

    return DerivedInstrumentResult::from(
        true,
        Some(overlapping_instruments),
        "high".to_string(),
        Some("instrument and flowcell id".to_string()),
        None,
    );
}

/// Main method to evaluate the detected instrument names and flowcell names and
/// return a result for the derived instruments. This may fail, and the
/// resulting `DerivedInstrumentResult` should be evaluated accordingly.
pub fn predict(
    instrument_names: HashSet<String>,
    flowcell_names: HashSet<String>,
) -> DerivedInstrumentResult {
    let instruments = super::instruments::build_instrument_lookup_table();
    let flowcells = super::flowcells::build_flowcell_lookup_table();

    let iid_results = predict_instrument(instrument_names, &instruments);
    let fcid_results = predict_instrument(flowcell_names, &flowcells);

    resolve_instrument_prediction(iid_results, fcid_results)
}

#[cfg(test)]
mod tests {
    use crate::derive::instrument::{flowcells, instruments};

    use super::*;

    #[test]
    fn test_derive_instrument_from_invalid_instrument_name() {
        let instruments = instruments::build_instrument_lookup_table();
        let result = possible_instruments_for_query(String::from("NoMatchingName"), &instruments);
        assert!(result.is_empty());
    }

    #[test]
    fn test_derive_instrument_from_valid_instrument_name() {
        let instruments = instruments::build_instrument_lookup_table();
        let result = possible_instruments_for_query(String::from("A00000"), &instruments);
        assert_eq!(result.len(), 1);
        assert!(result.contains("NovaSeq"));
    }

    #[test]
    fn test_derive_instrument_from_invalid_flowcell_name() {
        let flowcells = flowcells::build_flowcell_lookup_table();
        let result = possible_instruments_for_query(String::from("NoMatchingName"), &flowcells);
        assert!(result.is_empty());
    }

    #[test]
    fn test_derive_instrument_from_valid_flowcell_name() {
        let flowcells = flowcells::build_flowcell_lookup_table();
        let result = possible_instruments_for_query(String::from("H00000RXX"), &flowcells);
        assert_eq!(result.len(), 1);
        assert!(result.contains("NovaSeq"));
    }

    #[test]
    fn test_derive_instrument_novaseq_succesfully() {
        let detected_iids = HashSet::from(["A00000".to_string()]);
        let detected_fcids = HashSet::from(["H00000RXX".to_string()]);
        let result = predict(detected_iids, detected_fcids);

        assert_eq!(result.succeeded, true);
        assert_eq!(
            result.instruments,
            Some(HashSet::from(["NovaSeq".to_string()]))
        );
        assert_eq!(result.confidence, "high".to_string());
        assert_eq!(
            result.evidence,
            Some("instrument and flowcell id".to_string())
        );
        assert_eq!(result.comment, None);
    }

    #[test]
    fn test_derive_instrument_conflicting_instrument_ids() {
        let detected_iids = HashSet::from(["A00000".to_string(), "D00000".to_string()]);
        let detected_fcids = HashSet::from(["H00000RXX".to_string()]);
        let result = predict(detected_iids, detected_fcids);

        assert_eq!(result.succeeded, false);
        assert_eq!(result.instruments, None);
        assert_eq!(result.confidence, "unknown".to_string());
        assert_eq!(result.evidence, Some("instrument id".to_string()));
        assert_eq!(
            result.comment,
            Some(
                "multiple instruments were detected in this file via the instrument id".to_string()
            )
        );
    }

    #[test]
    fn test_derive_instrument_conflicting_flowcell_ids() {
        let detected_iids = HashSet::from(["A00000".to_string()]);
        let detected_fcids = HashSet::from(["H00000RXX".to_string(), "B0000".to_string()]);
        let result = predict(detected_iids, detected_fcids);

        assert_eq!(result.succeeded, false);
        assert_eq!(result.instruments, None);
        assert_eq!(result.confidence, "unknown".to_string());
        assert_eq!(result.evidence, Some("flowcell id".to_string()));
        assert_eq!(
            result.comment,
            Some("multiple instruments were detected in this file via the flowcell id".to_string())
        );
    }

    #[test]
    fn test_derive_instrument_medium_instrument_evidence() {
        let detected_iids = HashSet::from(["A00000".to_string()]);
        let detected_fcids = HashSet::new();
        let result = predict(detected_iids, detected_fcids);

        assert_eq!(result.succeeded, true);
        assert_eq!(
            result.instruments,
            Some(HashSet::from(["NovaSeq".to_string()]))
        );
        assert_eq!(result.confidence, "medium".to_string());
        assert_eq!(result.evidence, Some("instrument id".to_string()));
        assert_eq!(result.comment, None);
    }

    #[test]
    fn test_derive_instrument_low_instrument_evidence() {
        let detected_iids = HashSet::from(["K00000".to_string()]);
        let detected_fcids = HashSet::new();
        let result = predict(detected_iids, detected_fcids);

        assert_eq!(result.succeeded, true);
        assert_eq!(
            result.instruments,
            Some(HashSet::from([
                "HiSeq 4000".to_string(),
                "HiSeq 3000".to_string()
            ]))
        );
        assert_eq!(result.confidence, "low".to_string());
        assert_eq!(result.evidence, Some("instrument id".to_string()));
        assert_eq!(result.comment, None);
    }

    #[test]
    fn test_derive_instrument_medium_flowcell_evidence() {
        let detected_iids = HashSet::new();
        let detected_fcids = HashSet::from(["H00000RXX".to_string()]);
        let result = predict(detected_iids, detected_fcids);

        assert_eq!(result.succeeded, true);
        assert_eq!(
            result.instruments,
            Some(HashSet::from(["NovaSeq".to_string()]))
        );
        assert_eq!(result.confidence, "medium".to_string());
        assert_eq!(result.evidence, Some("flowcell id".to_string()));
        assert_eq!(result.comment, None);
    }

    #[test]
    fn test_derive_instrument_low_flowcell_evidence() {
        let detected_iids = HashSet::new();
        let detected_fcids = HashSet::from(["H0000ADXX".to_string()]);
        let result = predict(detected_iids, detected_fcids);

        assert_eq!(result.succeeded, true);
        assert_eq!(
            result.instruments,
            Some(HashSet::from([
                "HiSeq 2000".to_string(),
                "HiSeq 1500".to_string(),
                "HiSeq 2500".to_string()
            ]))
        );
        assert_eq!(result.confidence, "low".to_string());
        assert_eq!(result.evidence, Some("flowcell id".to_string()));
        assert_eq!(result.comment, None);
    }

    #[test]
    fn test_derive_instrument_conflicting_flowcell_and_instrument_evidence() {
        let detected_iids = HashSet::from(["K00000".to_string()]);
        let detected_fcids = HashSet::from(["H00000RXX".to_string()]);
        let result = predict(detected_iids, detected_fcids);

        assert_eq!(result.succeeded, false);
        assert_eq!(result.instruments, None);
        assert_eq!(result.confidence, "high".to_string());
        assert_eq!(
            result.evidence,
            Some("instrument and flowcell id".to_string())
        );
        assert_eq!(result.comment, Some("Case needs triaging, results from instrument id and flowcell id are mutually exclusive.".to_string()));
    }

    #[test]
    fn test_derive_instrument_no_matches() {
        let detected_iids = HashSet::from(["QQQQQ".to_string()]);
        let detected_fcids = HashSet::from(["ZZZZZZ".to_string()]);
        let result = predict(detected_iids, detected_fcids);

        assert_eq!(result.succeeded, false);
        assert_eq!(result.instruments, None);
        assert_eq!(result.confidence, "unknown".to_string());
        assert_eq!(result.evidence, None);
        assert_eq!(
            result.comment,
            Some("no matching instruments were found".to_string())
        );
    }
}
