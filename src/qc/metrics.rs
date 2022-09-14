use std::{fs::File, io::Write, path::PathBuf};

use noodles_bam::lazy::Record;
use serde::Serialize;

use super::{ComputationalLoad, Error, QualityCheckFacet};

#[derive(Debug, Serialize)]
pub struct ReadDesignation {
    primary: usize,
    secondary: usize,
    supplementary: usize,
}

#[derive(Debug, Serialize)]
pub struct ReadCountMetrics {
    total: usize,
    unmapped: usize,
    duplicates: usize,
    designations: ReadDesignation,
}

#[derive(Debug, Serialize)]
pub struct SummaryMetrics {
    duplication_pct: f64,
    unmapped_pct: f64,
}

#[derive(Debug, Serialize)]
pub struct SummaryMetricsFacet {
    read_count_metrics: ReadCountMetrics,
    summary: Option<SummaryMetrics>,
}

impl SummaryMetricsFacet {
    pub fn default() -> Self {
        SummaryMetricsFacet {
            read_count_metrics: ReadCountMetrics {
                total: 0,
                unmapped: 0,
                duplicates: 0,
                designations: ReadDesignation {
                    primary: 0,
                    secondary: 0,
                    supplementary: 0,
                },
            },
            summary: None,
        }
    }
}

impl QualityCheckFacet for SummaryMetricsFacet {
    fn name(&self) -> &'static str {
        "Summary Metrics"
    }

    fn computational_load(&self) -> ComputationalLoad {
        ComputationalLoad::Light
    }

    fn default(&self) -> bool {
        true
    }

    fn process(&mut self, record: &Record) -> Result<(), Error> {
        // (1) Count the number of reads in the file.
        self.read_count_metrics.total += 1;

        // (2) Compute metrics related to flags.
        if let Ok(s) = record.flags() {
            if s.is_duplicate() {
                self.read_count_metrics.duplicates += 1;
            }

            if s.is_unmapped() {
                self.read_count_metrics.unmapped += 1;
            }

            if s.is_secondary() {
                self.read_count_metrics.designations.secondary += 1;
            } else if s.is_supplementary() {
                self.read_count_metrics.designations.supplementary += 1;
            } else {
                self.read_count_metrics.designations.primary += 1;
            }
        }

        Ok(())
    }

    fn summarize(&mut self) -> Result<(), super::Error> {
        let summary = SummaryMetrics {
            duplication_pct: self.read_count_metrics.duplicates as f64
                / self.read_count_metrics.total as f64
                * 100.0,
            unmapped_pct: self.read_count_metrics.unmapped as f64
                / self.read_count_metrics.total as f64
                * 100.0,
        };

        self.summary = Some(summary);

        Ok(())
    }

    fn write(
        &self,
        output_prefix: String,
        directory: &std::path::Path,
    ) -> Result<(), std::io::Error> {
        let metrics_filename = output_prefix + ".metrics.json";
        let mut metrics_filepath = PathBuf::from(directory);
        metrics_filepath.push(metrics_filename);

        let mut file = File::create(metrics_filepath).unwrap();
        let output = serde_json::to_string_pretty(&self).unwrap();
        file.write_all(output.as_bytes())?;

        Ok(())
    }
}
impl SummaryMetricsFacet {}
