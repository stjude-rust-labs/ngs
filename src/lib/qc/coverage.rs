use std::collections::HashMap;

use noodles_sam::header::ReferenceSequence;
use serde::Serialize;

use crate::lib::utils::{genome::PRIMARY_CHROMOSOMES, histogram::SimpleHistogram};

use super::SequenceBasedQualityCheckFacet;

#[derive(Clone, Default, Serialize)]
pub struct CoverageMetrics {
    mean_coverage: HashMap<String, f64>,
    median_coverage: HashMap<String, f64>,
    median_over_mean_coverage: HashMap<String, f64>,
    ignored: HashMap<String, usize>,
    histograms: HashMap<String, SimpleHistogram>,
}

#[derive(Default)]
pub struct CoverageHistograms<'a> {
    storage: HashMap<&'a str, SimpleHistogram>,
}

#[derive(Default)]
pub struct CoverageFacet<'a> {
    by_position: CoverageHistograms<'a>,
    metrics: CoverageMetrics,
}

impl<'a> SequenceBasedQualityCheckFacet<'a> for CoverageFacet<'a> {
    fn name(&self) -> &'static str {
        "Coverage Metrics"
    }

    fn computational_load(&self) -> super::ComputationalLoad {
        super::ComputationalLoad::Heavy
    }

    fn supports_sequence_name(&self, name: &str) -> bool {
        PRIMARY_CHROMOSOMES.contains(&name)
    }

    fn process<'b>(
        &mut self,
        seq: &'b ReferenceSequence,
        record: &noodles_sam::alignment::Record,
    ) -> Result<(), super::Error>
    where
        'b: 'a,
    {
        let h = self
            .by_position
            .storage
            .entry(seq.name().as_str())
            .or_insert_with(|| SimpleHistogram::zero_based_with_capacity(usize::from(seq.len())));

        let record_start = usize::from(record.alignment_start().unwrap());
        let record_end = usize::from(record.alignment_end().unwrap());

        for i in record_start..=record_end {
            h.increment(i).unwrap();
        }

        Ok(())
    }

    fn summarize_seq(&mut self, seq: &ReferenceSequence) -> Result<(), super::Error> {
        let positions = self.by_position.storage.get(seq.name().as_str()).unwrap();
        let mut coverages = SimpleHistogram::zero_based_with_capacity(512);
        let mut ignored = 0;

        for i in positions.get_range_start()..=positions.get_range_stop() {
            if coverages.increment(positions.get(i)).is_err() {
                ignored += 1;
            }
        }

        let mean = coverages.mean();
        let median = coverages.median().unwrap();
        let median_over_mean = median / mean;

        // Saved for reporting.
        self.metrics
            .mean_coverage
            .insert(seq.name().to_string(), mean);
        self.metrics
            .median_coverage
            .insert(seq.name().to_string(), median);
        self.metrics
            .median_over_mean_coverage
            .insert(seq.name().to_string(), median_over_mean);
        self.metrics
            .histograms
            .insert(seq.name().to_string(), coverages);
        self.metrics.ignored.insert(seq.name().to_string(), ignored);

        // Removed to save memory.
        self.by_position.storage.remove(seq.name().as_str());
        Ok(())
    }

    fn aggregate_results(&self, results: &mut super::results::Results) {
        results.set_coverage(self.metrics.clone());
    }
}
