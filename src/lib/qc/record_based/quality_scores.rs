use std::{collections::HashMap, sync::Mutex};

use serde::{Deserialize, Serialize};

use crate::lib::{
    qc::{results, ComputationalLoad, RecordBasedQualityCheckFacet},
    utils::histogram::{Histogram, MutexBackedHistogram},
};

//=======================//
// Quality Score Metrics //
//=======================//

#[derive(Debug, Default)]
pub struct AtomicQualityScoreMetrics {
    scores: Mutex<HashMap<usize, MutexBackedHistogram>>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct QualityScoreMetrics {
    pub scores: HashMap<usize, Histogram>,
}

impl AtomicQualityScoreMetrics {
    pub fn clone_nonatomic(&self) -> QualityScoreMetrics {
        let scores = self.scores.lock().unwrap();
        let mut new_scores = HashMap::new();

        for (key, value) in scores.iter() {
            new_scores.insert(*key, value.clone_nonmutex());
        }

        QualityScoreMetrics { scores: new_scores }
    }
}

#[derive(Debug, Default)]
pub struct QualityScoreFacet {
    metrics: AtomicQualityScoreMetrics,
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

    fn process(&mut self, record: &noodles_sam::alignment::Record) -> anyhow::Result<()> {
        for (i, val) in record.quality_scores().as_ref().iter().enumerate() {
            let histogram = self
                .metrics
                .scores
                .lock()
                .unwrap()
                .entry(i + 1) // indices are 0-based, we want this to be 1-based.
                .or_insert_with(|| MutexBackedHistogram::zero_based_with_capacity(self::MAX_SCORE));

            let score = u8::from(*val) as usize;
            // histogram.increment(score).unwrap();
        }

        Ok(())
    }

    fn summarize(&mut self) -> anyhow::Result<()> {
        // Nothing to summarize here, as we simply report the histograms for
        // each position.

        Ok(())
    }

    fn aggregate(&self, results: &mut results::Results) {
        results.set_quality_scores(self.metrics.clone_nonatomic());
    }
}
