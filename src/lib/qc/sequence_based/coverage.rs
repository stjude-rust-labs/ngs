use std::{collections::HashMap, rc::Rc};

use noodles_sam::header::ReferenceSequence;
use serde::{Deserialize, Serialize};
use tracing::error;

use crate::lib::{
    qc::{results, ComputationalLoad, SequenceBasedQualityCheckFacet},
    utils::{
        genome::{get_primary_assembly, ReferenceGenome, Sequence},
        histogram::SimpleHistogram,
    },
};

#[derive(Clone, Default, Serialize, Deserialize)]
pub struct IgnoredMetrics {
    nonsensical_records: usize,
    pileup_too_large_positions: HashMap<String, usize>,
}

#[derive(Clone, Default, Serialize, Deserialize)]
pub struct CoverageMetrics {
    mean_coverage: HashMap<String, f64>,
    median_coverage: HashMap<String, f64>,
    median_over_mean_coverage: HashMap<String, f64>,
    ignored: IgnoredMetrics,
    histograms: HashMap<String, SimpleHistogram>,
}

#[derive(Default)]
pub struct CoverageHistograms<'a> {
    storage: HashMap<&'a str, SimpleHistogram>,
}

pub struct CoverageFacet<'a> {
    by_position: CoverageHistograms<'a>,
    metrics: CoverageMetrics,
    primary_assembly: Vec<Sequence>,
}

impl<'a> CoverageFacet<'a> {
    pub fn new(reference_genome: Rc<Box<dyn ReferenceGenome>>) -> Self {
        Self {
            by_position: CoverageHistograms::default(),
            metrics: CoverageMetrics::default(),
            primary_assembly: get_primary_assembly(reference_genome),
        }
    }
}

impl<'a> SequenceBasedQualityCheckFacet<'a> for CoverageFacet<'a> {
    fn name(&self) -> &'static str {
        "Coverage Metrics"
    }

    fn computational_load(&self) -> ComputationalLoad {
        ComputationalLoad::Moderate
    }

    fn supports_sequence_name(&self, name: &str) -> bool {
        self.primary_assembly
            .iter()
            .map(|s| s.name())
            .any(|x| x == name)
    }

    fn setup(&mut self, _: &ReferenceSequence) -> anyhow::Result<()> {
        Ok(())
    }

    fn process<'b>(
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
            if h.increment(i).is_err() {
                error!(
                    "Record crosses the sequence boundaries in an expected way. \
                    This usually means that the record is malformed. Please examine \
                    the record closely to ensure it fits within the sequence. \
                    Ignoring record. Read name: {}, Start Alignment: {}, End \
                    Alignment: {}, Cigar: {}",
                    record.read_name().unwrap(),
                    record.alignment_start().unwrap(),
                    record.alignment_end().unwrap(),
                    record.cigar()
                );
                self.metrics.ignored.nonsensical_records += 1;
            }
        }

        Ok(())
    }

    fn teardown(&mut self, sequence: &ReferenceSequence) -> anyhow::Result<()> {
        let positions = self
            .by_position
            .storage
            .get(sequence.name().as_str())
            .unwrap();
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
        self.by_position.storage.remove(sequence.name().as_str());

        // Saved for reporting.
        self.metrics
            .mean_coverage
            .insert(sequence.name().to_string(), mean);
        self.metrics
            .median_coverage
            .insert(sequence.name().to_string(), median);
        self.metrics
            .median_over_mean_coverage
            .insert(sequence.name().to_string(), median_over_mean);
        self.metrics
            .histograms
            .insert(sequence.name().to_string(), coverages);
        self.metrics
            .ignored
            .pileup_too_large_positions
            .insert(sequence.name().to_string(), ignored);

        Ok(())
    }

    fn aggregate(&mut self, results: &mut results::Results) {
        results.set_coverage(self.metrics.clone());
    }
}
