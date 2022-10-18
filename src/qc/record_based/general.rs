//! Functionality related to the General quality control facet.

use noodles::sam;
use sam::alignment::Record;

use crate::qc::results;
use crate::qc::ComputationalLoad;
use crate::qc::RecordBasedQualityControlFacet;

use self::metrics::GeneralMetrics;
pub use self::metrics::SummaryMetrics;

pub mod metrics;

/// Main struct for the General quality control facet.
#[derive(Clone, Debug, Default)]
pub struct GeneralMetricsFacet {
    /// The main metric counting struct.
    pub metrics: GeneralMetrics,
}

impl RecordBasedQualityControlFacet for GeneralMetricsFacet {
    fn name(&self) -> &'static str {
        "General"
    }

    fn computational_load(&self) -> ComputationalLoad {
        ComputationalLoad::Light
    }

    fn process(&mut self, record: &Record) -> anyhow::Result<()> {
        // (1) Count the number of reads in the file.
        self.metrics.records.total += 1;

        // (2) Compute metrics related to flags.
        let flags = record.flags();
        if flags.is_unmapped() {
            self.metrics.records.unmapped += 1;
        }

        if flags.is_duplicate() {
            self.metrics.records.duplicate += 1;
        }

        if flags.is_secondary() {
            self.metrics.records.designation.secondary += 1;
        } else if flags.is_supplementary() {
            self.metrics.records.designation.supplementary += 1;
        } else {
            self.metrics.records.designation.primary += 1;

            if !flags.is_unmapped() {
                self.metrics.records.primary_mapped += 1;
            }

            if flags.is_duplicate() {
                self.metrics.records.primary_duplicate += 1;
            }

            if flags.is_segmented() {
                self.metrics.records.paired += 1;

                if flags.is_first_segment() {
                    self.metrics.records.read_1 += 1;
                }

                if flags.is_last_segment() {
                    self.metrics.records.read_2 += 1;
                }

                if !flags.is_unmapped() {
                    if flags.is_properly_aligned() {
                        self.metrics.records.proper_pair += 1;
                    }

                    if flags.is_mate_unmapped() {
                        self.metrics.records.singleton += 1;
                    } else {
                        self.metrics.records.mate_mapped += 1;

                        let reference_sequence_id = record.reference_sequence_id().unwrap();
                        let mate_reference_sequence_id =
                            record.mate_reference_sequence_id().unwrap();

                        if reference_sequence_id != mate_reference_sequence_id {
                            self.metrics.records.mate_reference_sequence_id_mismatch += 1;

                            let mapq = record
                                .mapping_quality()
                                .map(u8::from)
                                .unwrap_or(sam::record::mapping_quality::MISSING);

                            if mapq >= 5 {
                                self.metrics.records.mate_reference_sequence_id_mismatch_hq += 1;
                            }
                        }
                    }
                }
            }
        }

        // (3) Compute CIGAR accumulations
        let cigar = record.cigar();
        let read_one = record.flags().is_first_segment();
        for op in cigar.iter() {
            if read_one {
                *self
                    .metrics
                    .cigar
                    .read_one_cigar_ops
                    .entry(op.kind().to_string())
                    .or_insert(0) += 1
            } else {
                *self
                    .metrics
                    .cigar
                    .read_two_cigar_ops
                    .entry(op.kind().to_string())
                    .or_insert(0) += 1
            }
        }

        Ok(())
    }

    fn summarize(&mut self) -> anyhow::Result<()> {
        let summary = SummaryMetrics {
            duplication_pct: self.metrics.records.duplicate as f64
                / self.metrics.records.total as f64
                * 100.0,
            mapped_pct: (1.0
                - self.metrics.records.unmapped as f64 / self.metrics.records.total as f64)
                * 100.0,
            mate_reference_sequence_id_mismatch_pct: self
                .metrics
                .records
                .mate_reference_sequence_id_mismatch
                as f64
                / self.metrics.records.total as f64
                * 100.0,
            mate_reference_sequence_id_mismatch_hq_pct: self
                .metrics
                .records
                .mate_reference_sequence_id_mismatch_hq
                as f64
                / self.metrics.records.total as f64
                * 100.0,
        };

        self.metrics.summary = Some(summary);

        Ok(())
    }

    fn aggregate(&self, results: &mut results::Results) {
        results.general = Some(self.metrics.clone())
    }
}
