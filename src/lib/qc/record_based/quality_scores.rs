use std::collections::HashMap;

use serde::Serialize;

use crate::lib::{
    qc::{results, ComputationalLoad, Error, RecordBasedQualityCheckFacet},
    utils::histogram::SimpleHistogram,
};

#[derive(Clone, Debug, Default, Serialize)]
pub struct QualityScoreFacet {
    scores: HashMap<usize, SimpleHistogram>,
}

const MAX_SCORE: usize = 93;
// Hopefully in the future, we can do something like this if `noodles` makes
// this a public accessible const.
// const MAX_SCORE: usize = sam::record::quality_scores::score::MAX as usize;

impl RecordBasedQualityCheckFacet for QualityScoreFacet {
    fn name(&self) -> &'static str {
        "Quality Score Metrics"
    }

    fn computational_load(&self) -> ComputationalLoad {
        ComputationalLoad::Moderate
    }

    fn process(&mut self, record: &noodles_sam::alignment::Record) -> Result<(), Error> {
        for (i, val) in record.quality_scores().as_ref().iter().enumerate() {
            let histogram = self
                .scores
                .entry(i + 1) // indices are 0-based, we want this to be 1-based.
                .or_insert_with(|| SimpleHistogram::zero_based_with_capacity(self::MAX_SCORE));

            let score = u8::from(*val) as usize;
            histogram.increment(score).unwrap();
        }

        Ok(())
    }

    fn summarize(&mut self) -> Result<(), Error> {
        // Nothing to summarize here, as we simply report the histograms for
        // each position.

        Ok(())
    }

    fn aggregate_results(&self, results: &mut results::Results) {
        results.set_quality_scores(self.clone());
    }
}
