//! Functionality related to the Quality Scores quality control facet.

use std::collections::HashMap;

use noodles::sam::alignment::Record;
use serde::{Deserialize, Serialize};

use crate::{
    qc::{results, ComputationalLoad, RecordBasedQualityControlFacet},
    utils::histogram::Histogram,
};

/// Main struct for the Quality Scores quality control facet.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct QualityScoreFacet {
    /// Distribution of quality scores for each position in the records observed.
    pub scores: HashMap<usize, Histogram>,
}

/// Maximum quality score supported by the SAM specification.
///
/// Hopefully in the future, we can do something like this if `noodles` makes
/// this a public accessible const.
/// `const MAX_SCORE: usize = sam::record::quality_scores::score::MAX as usize;`
pub const MAX_SCORE: usize = 93;

impl RecordBasedQualityControlFacet for QualityScoreFacet {
    fn name(&self) -> &'static str {
        "Quality Score"
    }

    fn computational_load(&self) -> ComputationalLoad {
        ComputationalLoad::Moderate
    }

    fn process(&mut self, record: &Record) -> anyhow::Result<()> {
        for (i, val) in record.quality_scores().as_ref().iter().enumerate() {
            let histogram = self
                .scores
                .entry(i + 1) // indices are 0-based, we want this to be 1-based.
                .or_insert_with(|| Histogram::zero_based_with_capacity(self::MAX_SCORE));

            let score = u8::from(*val) as usize;
            histogram.increment(score).unwrap();
        }

        Ok(())
    }

    fn summarize(&mut self) -> anyhow::Result<()> {
        // Nothing to summarize here, as we simply report the histograms for
        // each position.

        Ok(())
    }

    fn aggregate(&self, results: &mut results::Results) {
        results.quality_scores = Some(self.clone());
    }
}
