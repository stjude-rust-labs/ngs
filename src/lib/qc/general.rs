use noodles_bam::lazy::Record;
use noodles_sam as sam;

pub use self::metrics::{GeneralMetricsFacet, SummaryMetrics};

use super::{results, ComputationalLoad, Error, RecordBasedQualityCheckFacet};

pub mod metrics;

impl RecordBasedQualityCheckFacet for GeneralMetricsFacet {
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

                            let reference_sequence_id =
                                record.reference_sequence_id().unwrap().unwrap();
                            let mate_reference_sequence_id =
                                record.mate_reference_sequence_id().unwrap().unwrap();

                            if reference_sequence_id != mate_reference_sequence_id {
                                self.records.mate_reference_sequence_id_mismatch += 1;

                                let mapq = record
                                    .mapping_quality()
                                    .map(|x| x.unwrap())
                                    .map(u8::from)
                                    .unwrap_or(sam::record::mapping_quality::MISSING);

                                if mapq >= 5 {
                                    self.records.mate_reference_sequence_id_mismatch_hq += 1;
                                }
                            }
                        }
                    }
                }
            }
        }

        // (3) Compute CIGAR accumulations
        let cigar: sam::record::Cigar = record.cigar().try_into().unwrap();
        let read_one = record.flags().unwrap().is_first_segment();
        for op in cigar.iter() {
            if read_one {
                *self
                    .records
                    .read_one_cigar_ops
                    .entry(op.kind().to_string())
                    .or_insert(0) += 1
            } else {
                *self
                    .records
                    .read_two_cigar_ops
                    .entry(op.kind().to_string())
                    .or_insert(0) += 1
            }
        }
        Ok(())
    }

    fn summarize(&mut self) -> Result<(), super::Error> {
        let summary = SummaryMetrics {
            duplication_pct: self.records.duplicate as f64 / self.records.total as f64 * 100.0,
            unmapped_pct: self.records.unmapped as f64 / self.records.total as f64 * 100.0,
            mate_reference_sequence_id_mismatch_pct: self
                .records
                .mate_reference_sequence_id_mismatch
                as f64
                / self.records.total as f64
                * 100.0,
            mate_reference_sequence_id_mismatch_hq_pct: self
                .records
                .mate_reference_sequence_id_mismatch_hq
                as f64
                / self.records.total as f64
                * 100.0,
        };

        self.summary = Some(summary);

        Ok(())
    }

    fn aggregate_results(&self, results: &mut results::Results) {
        results.set_general(self.clone())
    }
}
