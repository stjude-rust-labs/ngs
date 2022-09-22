use noodles_bam::lazy::Record;

pub use self::metrics::{GeneralMetricsFacet, SummaryMetrics};

use super::{results, ComputationalLoad, Error, QualityCheckFacet};

pub mod metrics;

impl QualityCheckFacet for GeneralMetricsFacet {
    fn name(&self) -> &'static str {
        "General Metrics"
    }

    fn computational_load(&self) -> ComputationalLoad {
        ComputationalLoad::Light
    }

    fn process(&mut self, record: &Record) -> Result<(), Error> {
        // (1) Count the number of reads in the file.
        self.records.total += 1;

        // (2) Compute metrics related to flags.
        if let Ok(flags) = record.flags() {
            if flags.is_unmapped() {
                self.records.unmapped += 1;
            }

            if flags.is_duplicate() {
                self.records.duplicate += 1;
            }

            if flags.is_secondary() {
                self.records.designation.secondary += 1;
            } else if flags.is_supplementary() {
                self.records.designation.supplementary += 1;
            } else {
                self.records.designation.primary += 1;

                if !flags.is_unmapped() {
                    self.records.primary_mapped += 1;
                }

                if flags.is_duplicate() {
                    self.records.primary_duplicate += 1;
                }

                if flags.is_segmented() {
                    self.records.paired += 1;

                    if flags.is_first_segment() {
                        self.records.read_1 += 1;
                    }

                    if flags.is_last_segment() {
                        self.records.read_2 += 1;
                    }

                    if !flags.is_unmapped() {
                        if flags.is_properly_aligned() {
                            self.records.proper_pair += 1;
                        }

                        if flags.is_mate_unmapped() {
                            self.records.singleton += 1;
                        } else {
                            self.records.mate_mapped += 1;
                        }
                    }
                }
            }
        }

        Ok(())
    }

    fn summarize(&mut self) -> Result<(), super::Error> {
        let summary = SummaryMetrics {
            duplication_pct: self.records.duplicate as f64 / self.records.total as f64 * 100.0,
            unmapped_pct: self.records.unmapped as f64 / self.records.total as f64 * 100.0,
        };

        self.summary = Some(summary);

        Ok(())
    }

    fn aggregate_results(&self, results: &mut results::Results) {
        results.set_general(self.clone())
    }
}
