//! Functionality related to the Coverage quality control facet.

use std::{collections::HashMap, rc::Rc};

use noodles::sam::{
    alignment::Record,
    header::record::value::{map::ReferenceSequence, Map},
};
use serde::{Deserialize, Serialize};
use tracing::error;

use crate::{
    qc::{results, ComputationalLoad, SequenceBasedQualityCheckFacet},
    utils::{
        genome::{get_primary_assembly, ReferenceGenome, Sequence},
        histogram::SimpleHistogram,
    },
};

//=========//
// Metrics //
//=========//

/// Metrics related to ignored records or positions in the genome.
#[derive(Clone, Default, Serialize, Deserialize)]
pub struct IgnoredMetrics {
    /// The number of records that were considered non-sensical by this quality
    /// control facet.
    pub nonsensical_records: usize,

    /// The number of positions which had coverage that was too high such that
    /// it did not fit into our histogram.
    pub pileup_too_large_positions: HashMap<String, usize>,
}

/// Primary struct used to compile stats regarding coverage.
#[derive(Clone, Default, Serialize, Deserialize)]
pub struct CoverageMetrics {
    /// Hashmap containing the mean coverage for each sequence in the reference
    /// genome.
    pub mean_coverage: HashMap<String, f64>,

    /// Hashmap containing the median coverage for each sequence in the
    /// reference genome.
    pub median_coverage: HashMap<String, f64>,

    /// Hashmap containing the median over mean coverage for each sequence in
    /// the reference genome.
    pub median_over_mean_coverage: HashMap<String, f64>,

    /// Metrics recording various records or positions that were ignored during
    /// the analysis.
    pub ignored: IgnoredMetrics,

    /// Coverage distribution as a histogram per sequence.
    pub coverage_distribution_per_sequence: HashMap<String, SimpleHistogram>,
}

/// Main struct for the Coverage quality control facet.
pub struct CoverageFacet {
    /// Data structure for tallying up coverage across position for every
    /// sequence in the reference genome. This is eventually discarded in favor
    /// of a coverage distribution for the sequence.
    coverage_per_position: HashMap<String, SimpleHistogram>,

    /// Metrics related to the Coverage quality control facet.
    metrics: CoverageMetrics,

    /// Struct for caching all of the sequences considered part of the primary
    /// assembly.
    primary_assembly: Vec<Sequence>,
}

impl CoverageFacet {
    /// Creates a new [`CoverageFacet`].
    pub fn new(reference_genome: Rc<Box<dyn ReferenceGenome>>) -> Self {
        Self {
            coverage_per_position: HashMap::default(),
            metrics: CoverageMetrics::default(),
            primary_assembly: get_primary_assembly(reference_genome),
        }
    }
}

impl SequenceBasedQualityCheckFacet for CoverageFacet {
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

    fn setup(&mut self, _: &Map<ReferenceSequence>) -> anyhow::Result<()> {
        Ok(())
    }

    fn process<'b>(
        &mut self,
        seq: &'b Map<ReferenceSequence>,
        record: &Record,
    ) -> anyhow::Result<()> {
        let h = self
            .coverage_per_position
            .entry(seq.name().to_string())
            .or_insert_with(|| {
                SimpleHistogram::zero_based_with_capacity(usize::from(seq.length()))
            });

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

    fn teardown(&mut self, sequence: &Map<ReferenceSequence>) -> anyhow::Result<()> {
        let positions = match self.coverage_per_position.get(sequence.name().as_str()) {
            Some(s) => s,
            // In the None case, no records were inserted for this sequence.
            // This may be because the file is a mini-SAM/BAM/CRAM. If that's
            // the case, we just return Ok(()).
            None => return Ok(()),
        };

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
        self.coverage_per_position.remove(sequence.name().as_str());

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
            .coverage_distribution_per_sequence
            .insert(sequence.name().to_string(), coverages);
        self.metrics
            .ignored
            .pileup_too_large_positions
            .insert(sequence.name().to_string(), ignored);

        Ok(())
    }

    fn aggregate(&mut self, results: &mut results::Results) {
        results.coverage = Some(self.metrics.clone());
    }
}
