//! Combines the flowcell and instrument checks into a single workflow.

use regex::Regex;
use serde::Serialize;
use std::collections::{HashMap, HashSet};

use crate::derive::instrument::{flowcells, instruments};

/// Generalized struct for holding instrument detection results.
#[derive(Debug, Default, Serialize)]
pub struct InstrumentDetectionResults {
    /// The possible instruments contained within this result set.
    pub possible_instruments: Option<HashSet<String>>,
    /// Whether or not at least one machine has been detected.
    pub detected_at_least_one_machine: bool,
}

impl InstrumentDetectionResults {
    /// Updates the `InstrumentDetectionResults` with a set of instruments that
    /// have been detected.
    pub fn update_instruments(&mut self, instruments: &HashSet<String>) {
        self.possible_instruments = Some(match &self.possible_instruments {
            // An initial base set has already been established, so take the
            // intersection of the existing possible instruments set and the
            // set being passed into the update function.
            Some(r) => r.intersection(instruments).cloned().collect(),
            // This is the first iteration, so we need to set our base set as
            // the first detected set of results.
            None => instruments.clone(),
        });

        // After we've updated the sets, we need to keep track of if we have
        // detected at least one machine to distinguish conflicting machines
        // from no machines detected at all.
        if !instruments.is_empty() {
            self.detected_at_least_one_machine = true;
        }
    }
}

/// A query for a look-up table and the resulting hits from that table.
#[derive(Debug, Serialize)]
pub struct QueryResult {
    /// The query that was used to generate the result.
    pub query: String,

    /// The possible instruments that could have generated the query.
    pub result: HashSet<String>,
}

/// Utility struct for holding the results of look-ups for instrument and
/// flowcell names.
pub struct Queries {
    /// The results of the instrument name look-ups.
    pub instrument_name_queries: Vec<QueryResult>,

    /// The results of the flowcell name look-ups.
    pub flowcell_name_queries: Vec<QueryResult>,
}

/// Metrics related to how read records were processed.
#[derive(Debug, Default, Serialize)]
pub struct RecordMetrics {
    /// The total number of records that were processed.
    pub total_records: usize,

    /// The total number of records that couldn't be parsed
    /// due to a missing or invalid read name.
    pub bad_read_name: usize,

    /// The total number of records that contained a valid
    /// instrument name in their read name.
    pub found_instrument_name: usize,

    /// The total number of records that contained a valid
    /// flowcell name in their read name.
    pub found_flowcell_name: usize,
}

/// Struct holding the final results for an `ngs derive instrument` subcommand
/// call.
#[derive(Debug, Serialize)]
pub struct DerivedInstrumentResult {
    /// Whether or not the `ngs derive instrument` subcommand succeeded.
    pub succeeded: bool,

    /// The possible instruments detected by `ngs derive instrument`, if
    /// derivable.
    pub instruments: Option<HashSet<String>>,

    /// The level of confidence that the tool has concerning these results.
    pub confidence: String,

    /// Status of the evidence that supports (or lack thereof) these predicted
    /// instruments, if available.
    pub evidence: Option<String>,

    /// A general comment field, if available.
    pub comment: Option<String>,

    /// The results of the instrument name look-ups.
    pub instrument_name_queries: Vec<QueryResult>,

    /// The results of the flowcell name look-ups.
    pub flowcell_name_queries: Vec<QueryResult>,

    /// Metrics related to how read records were processed.
    pub records: RecordMetrics,
}

impl DerivedInstrumentResult {
    /// Creates a new [`DerivedInstrumentResult`].
    pub fn new(
        succeeded: bool,
        instruments: Option<HashSet<String>>,
        confidence: String,
        evidence: Option<String>,
        comment: Option<String>,
        queries: Queries,
        records: RecordMetrics,
    ) -> Self {
        DerivedInstrumentResult {
            succeeded,
            instruments,
            confidence,
            evidence,
            comment,
            instrument_name_queries: queries.instrument_name_queries,
            flowcell_name_queries: queries.flowcell_name_queries,
            records,
        }
    }
}

impl Default for DerivedInstrumentResult {
    fn default() -> Self {
        DerivedInstrumentResult {
            succeeded: false,
            instruments: None,
            confidence: "unknown".to_string(),
            evidence: None,
            comment: None,
            instrument_name_queries: Vec::new(),
            flowcell_name_queries: Vec::new(),
            records: RecordMetrics::default(),
        }
    }
}

/// Computes the full set of possible instruments that could have generated the
/// value passed to the function given the lookup table. This is intended to be
/// a general purpose method that will work with both flowcells and instrument
/// ids.
///
/// The `lookup_table` passed to this function is constructed as a set of regex
/// keys that map to machine that could have generated a query that matches that
/// regex. Effectively, this method iterates through all of the keys, checks if
/// the query matches the regex, and extends the result HashSet with the values
/// for that key if it does.
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
) -> QueryResult {
    let mut result_set: HashSet<String> = HashSet::new();

    for (pattern, machines) in lookup_table {
        let re = Regex::new(pattern).unwrap();
        if re.is_match(query.as_str()) {
            let matching_machines: Vec<String> = machines.iter().map(|x| x.to_string()).collect();
            result_set.extend(matching_machines);
        }
    }
    QueryResult {
        query,
        result: result_set,
    }
}

/// Given a HashSet of unique queries (usually a instrument ID or flowcell ID
/// parsed from a read name) that were detected from a SAM/BAM/CRAM file, return
/// a HashSet that contains all possible machines that could have generated that
/// list of queries and a vec recording the query look-ups that were made.
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
) -> (InstrumentDetectionResults, Vec<QueryResult>) {
    let mut result = InstrumentDetectionResults::default();
    let mut query_results = Vec::new();

    for name in queries {
        let derived = possible_instruments_for_query(name, lookup_table);
        result.update_instruments(&derived.result);
        query_results.push(derived);
    }

    (result, query_results)
}

/// Combines evidence from the instrument id detection and flowcell id detection
/// to produce a [`DerivedInstrumentResult`].
pub fn resolve_instrument_prediction(
    iid_results: InstrumentDetectionResults,
    fcid_results: InstrumentDetectionResults,
) -> DerivedInstrumentResult {
    let possible_instruments_by_iid = iid_results.possible_instruments.unwrap_or_default();
    let possible_instruments_by_fcid = fcid_results.possible_instruments.unwrap_or_default();

    let mut result = DerivedInstrumentResult::default();

    // (1) If the set of possible instruments as determined by the instrument id
    // is empty _and_ we have seen at least one machine, then the only possible
    // scenario is there are conflicting instrument ids.
    if possible_instruments_by_iid.is_empty() && iid_results.detected_at_least_one_machine {
        result.evidence = Some("instrument id".to_string());
        result.comment = Some(
            "multiple instruments were detected in this file via the instrument id".to_string(),
        );
        return result;
    }

    // (2) If the set of possible instruments as determined by the flowcell id
    // is empty _and_ we have seen at least one machine, then the only possible
    // scenario is there are conflicting flowcell ids.
    if possible_instruments_by_fcid.is_empty() && fcid_results.detected_at_least_one_machine {
        result.evidence = Some("flowcell id".to_string());
        result.comment =
            Some("multiple instruments were detected in this file via the flowcell id".to_string());
        return result;
    }

    // (3) if neither result turns up anything, then we can simply say that the
    // machine was not able to be detected.
    if possible_instruments_by_iid.is_empty() && possible_instruments_by_fcid.is_empty() {
        result.comment = Some("no matching instruments were found".to_string());
        return result;
    }

    // (4) If both aren't empty and iid_results _is_ empty, then the fcid
    // results must not be empty. We can go ahead and issue a prediction based
    // on this medium to low confidence result.
    if possible_instruments_by_iid.is_empty() {
        let instruments = possible_instruments_by_fcid;
        let confidence = match instruments.len() {
            1 => "medium",
            _ => "low",
        };

        result.succeeded = true;
        result.instruments = Some(instruments);
        result.confidence = confidence.to_string();
        result.evidence = Some("flowcell id".to_string());
        return result;
    }

    // (5) Same as the block above, except now we are evaluating the opposite
    // (only the iid results contained some predicted machine).
    if possible_instruments_by_fcid.is_empty() {
        let instruments = possible_instruments_by_iid;
        let confidence = match instruments.len() {
            1 => "medium",
            _ => "low",
        };

        result.succeeded = true;
        result.instruments = Some(instruments);
        result.confidence = confidence.to_string();
        result.evidence = Some("instrument id".to_string());
        return result;
    }

    let overlapping_instruments: HashSet<String> = possible_instruments_by_fcid
        .intersection(&possible_instruments_by_iid)
        .cloned()
        .collect();

    if overlapping_instruments.is_empty() {
        result.confidence = "high".to_string();
        result.evidence = Some("instrument and flowcell id".to_string());
        result.comment = Some(
            "Case needs triaging, results from instrument id and \
                         flowcell id are mutually exclusive."
                .to_string(),
        );
        return result;
    }

    result.succeeded = true;
    result.instruments = Some(overlapping_instruments);
    result.confidence = "high".to_string();
    result.evidence = Some("instrument and flowcell id".to_string());
    result
}

/// Main method to evaluate the detected instrument names and flowcell names and
/// return a result for the derived instruments. This may fail, and the
/// resulting [`DerivedInstrumentResult`] should be evaluated accordingly.
pub fn predict(
    instrument_names: HashSet<String>,
    flowcell_names: HashSet<String>,
) -> DerivedInstrumentResult {
    let instruments = instruments::build_instrument_lookup_table();
    let flowcells = flowcells::build_flowcell_lookup_table();

    let (iid_results, instrument_name_queries) = predict_instrument(instrument_names, &instruments);
    let (fcid_results, flowcell_name_queries) = predict_instrument(flowcell_names, &flowcells);

    let mut final_results = resolve_instrument_prediction(iid_results, fcid_results);
    final_results.instrument_name_queries = instrument_name_queries;
    final_results.flowcell_name_queries = flowcell_name_queries;
    final_results
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_derive_instrument_from_invalid_instrument_name() {
        let instruments = instruments::build_instrument_lookup_table();
        let result = possible_instruments_for_query(String::from("NoMatchingName"), &instruments);
        assert!(result.result.is_empty());
    }

    #[test]
    fn test_derive_instrument_from_valid_instrument_name() {
        let instruments = instruments::build_instrument_lookup_table();
        let result = possible_instruments_for_query(String::from("A00000"), &instruments);
        assert_eq!(result.result.len(), 1);
        assert!(result.result.contains("NovaSeq"));
    }

    #[test]
    fn test_derive_instrument_from_invalid_flowcell_name() {
        let flowcells = flowcells::build_flowcell_lookup_table();
        let result = possible_instruments_for_query(String::from("NoMatchingName"), &flowcells);
        assert!(result.result.is_empty());
    }

    #[test]
    fn test_derive_instrument_from_valid_flowcell_name() {
        let flowcells = flowcells::build_flowcell_lookup_table();
        let result = possible_instruments_for_query(String::from("H00000RXX"), &flowcells);
        assert_eq!(result.result.len(), 1);
        assert!(result.result.contains("NovaSeq"));
    }

    #[test]
    fn test_derive_instrument_novaseq_succesfully() {
        let detected_iids = HashSet::from(["A00000".to_string()]);
        let detected_fcids = HashSet::from(["H00000RXX".to_string()]);
        let result = predict(detected_iids, detected_fcids);

        assert!(result.succeeded);
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

        assert!(!result.succeeded);
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

        assert!(!result.succeeded);
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

        assert!(result.succeeded);
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

        assert!(result.succeeded);
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

        assert!(result.succeeded);
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

        assert!(result.succeeded);
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

        assert!(!result.succeeded);
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

        assert!(!result.succeeded);
        assert_eq!(result.instruments, None);
        assert_eq!(result.confidence, "unknown".to_string());
        assert_eq!(result.evidence, None);
        assert_eq!(
            result.comment,
            Some("no matching instruments were found".to_string())
        );
    }
}
