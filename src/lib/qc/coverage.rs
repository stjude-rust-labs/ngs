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
        super::ComputationalLoad::Moderate
    }

    fn supports_sequence_name(&self, name: &str) -> bool {
        PRIMARY_CHROMOSOMES.contains(&name)
    }

    fn setup_sequence(&mut self, _: &ReferenceSequence) -> anyhow::Result<()> {
        Ok(())
    }

    fn process_record<'b>(
        &mut self,
        seq: &'b ReferenceSequence,
        record: &noodles_sam::alignment::Record,
    ) -> anyhow::Result<()>
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

    fn teardown_sequence(&mut self, seq: &ReferenceSequence) -> anyhow::Result<()> {
        let positions = self.by_position.storage.get(seq.name().as_str()).unwrap();
        let mut coverages = SimpleHistogram::zero_based_with_capacity(1024);
        let mut ignored = 0;

        for i in positions.get_range_start()..=positions.get_range_stop() {
            if coverages.increment(positions.get(i)).is_err() {
                ignored += 1;
            }
        }

        let mean = coverages.mean();
        let median = coverages.median().unwrap();
        let median_over_mean = median / mean;

        // Removed to save memory.
        self.by_position.storage.remove(seq.name().as_str());

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

        Ok(())
    }

    fn aggregate_results(&mut self, results: &mut super::results::Results) {
        results.set_coverage(self.metrics.clone());
    }
}
