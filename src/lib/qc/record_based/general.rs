use std::sync::atomic::Ordering;

use noodles_sam::{self as sam};
use sam::alignment::Record;

use crate::lib::qc::{results, ComputationalLoad, RecordBasedQualityCheckFacet};

pub use self::metrics::{GeneralMetricsFacet, SummaryMetrics};

pub mod metrics;

impl RecordBasedQualityCheckFacet for GeneralMetricsFacet {
    fn name(&self) -> &'static str {
        "General Metrics"
    }

    fn computational_load(&self) -> ComputationalLoad {
        ComputationalLoad::Light
    }

    fn process(&mut self, record: &Record) -> anyhow::Result<()> {
        // (1) Count the number of reads in the file.
        self.metrics.records.total.fetch_add(1, Ordering::SeqCst);

        // (2) Compute metrics related to flags.
        let flags = record.flags();
        if flags.is_unmapped() {
            self.metrics.records.unmapped.fetch_add(1, Ordering::SeqCst);
        }

        if flags.is_duplicate() {
            self.metrics
                .records
                .duplicate
                .fetch_add(1, Ordering::SeqCst);
        }

        if flags.is_secondary() {
            self.metrics
                .records
                .designation
                .secondary
                .fetch_add(1, Ordering::SeqCst);
        } else if flags.is_supplementary() {
            self.metrics
                .records
                .designation
                .supplementary
                .fetch_add(1, Ordering::SeqCst);
        } else {
            self.metrics
                .records
                .designation
                .primary
                .fetch_add(1, Ordering::SeqCst);

            if !flags.is_unmapped() {
                self.metrics
                    .records
                    .primary_mapped
                    .fetch_add(1, Ordering::SeqCst);
            }

            if flags.is_duplicate() {
                self.metrics
                    .records
                    .primary_duplicate
                    .fetch_add(1, Ordering::SeqCst);
            }

            if flags.is_segmented() {
                self.metrics.records.paired.fetch_add(1, Ordering::SeqCst);

                if flags.is_first_segment() {
                    self.metrics.records.read_1.fetch_add(1, Ordering::SeqCst);
                }

                if flags.is_last_segment() {
                    self.metrics.records.read_2.fetch_add(1, Ordering::SeqCst);
                }

                if !flags.is_unmapped() {
                    if flags.is_properly_aligned() {
                        self.metrics
                            .records
                            .proper_pair
                            .fetch_add(1, Ordering::SeqCst);
                    }

                    if flags.is_mate_unmapped() {
                        self.metrics
                            .records
                            .singleton
                            .fetch_add(1, Ordering::SeqCst);
                    } else {
                        self.metrics
                            .records
                            .mate_mapped
                            .fetch_add(1, Ordering::SeqCst);

                        let reference_sequence_id = record.reference_sequence_id().unwrap();
                        let mate_reference_sequence_id =
                            record.mate_reference_sequence_id().unwrap();

                        if reference_sequence_id != mate_reference_sequence_id {
                            self.metrics
                                .records
                                .mate_reference_sequence_id_mismatch
                                .fetch_add(1, Ordering::SeqCst);

                            let mapq = record
                                .mapping_quality()
                                .map(u8::from)
                                .unwrap_or(sam::record::mapping_quality::MISSING);

                            if mapq >= 5 {
                                self.metrics
                                    .records
                                    .mate_reference_sequence_id_mismatch_hq
                                    .fetch_add(1, Ordering::SeqCst);
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
                    .lock()
                    .unwrap()
                    .entry(op.kind().to_string())
                    .or_insert(0) += 1
            } else {
                *self
                    .metrics
                    .cigar
                    .read_two_cigar_ops
                    .lock()
                    .unwrap()
                    .entry(op.kind().to_string())
                    .or_insert(0) += 1
            }
        }

        Ok(())
    }

    fn summarize(&mut self) -> anyhow::Result<()> {
        let total = self.metrics.records.total.load(Ordering::SeqCst);
        let duplicate = self.metrics.records.duplicate.load(Ordering::SeqCst);
        let unmapped = self.metrics.records.unmapped.load(Ordering::SeqCst);
        let mate_reference_sequence_id_mismatch = self
            .metrics
            .records
            .mate_reference_sequence_id_mismatch
            .load(Ordering::SeqCst);

        let mate_reference_sequence_id_mismatch_hq = self
            .metrics
            .records
            .mate_reference_sequence_id_mismatch_hq
            .load(Ordering::SeqCst);

        let summary = SummaryMetrics {
            duplication_pct: duplicate as f64 / total as f64 * 100.0,
            unmapped_pct: unmapped as f64 / total as f64 * 100.0,
            mate_reference_sequence_id_mismatch_pct: mate_reference_sequence_id_mismatch as f64
                / total as f64
                * 100.0,
            mate_reference_sequence_id_mismatch_hq_pct: mate_reference_sequence_id_mismatch_hq
                as f64
                / total as f64
                * 100.0,
        };

        self.metrics.summary = Some(summary);

        Ok(())
    }

    fn aggregate(&self, results: &mut results::Results) {
        results.set_general(self.metrics.clone_nonatomic())
    }
}
