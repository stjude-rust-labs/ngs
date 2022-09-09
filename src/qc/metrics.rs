use noodles_bam::lazy::Record;
use serde::Serialize;

use crate::qc::template_length::TemplateLengthHistogram;

use super::gc_content::GCContentHistogram;

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
    gc_content_pct: f64,
    template_length_unknown_pct: f64,
    template_length_out_of_range_pct: f64,
}

#[derive(Debug, Serialize)]
pub struct QualityCheckMetrics {
    read_count_metrics: ReadCountMetrics,
    summary_metrics: Option<SummaryMetrics>,
}

impl QualityCheckMetrics {
    pub fn new() -> Self {
        QualityCheckMetrics {
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
            summary_metrics: None,
        }
    }

    pub fn process(&mut self, record: &Record) {
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
    }

    pub fn summarize(
        &mut self,
        template_length_histogram: &TemplateLengthHistogram,
        gc_content: &GCContentHistogram,
    ) {
        let summary = SummaryMetrics {
            duplication_pct: self.read_count_metrics.duplicates as f64
                / self.read_count_metrics.total as f64
                * 100.0,
            unmapped_pct: self.read_count_metrics.unmapped as f64
                / self.read_count_metrics.total as f64
                * 100.0,
            gc_content_pct: (gc_content.get_gc_count() as f64
                / (gc_content.get_gc_count()
                    + gc_content.get_at_count()
                    + gc_content.get_other_count()) as f64)
                * 100.0,
            template_length_unknown_pct: template_length_histogram.get(0) as f64
                / self.read_count_metrics.total as f64
                * 100.0,
            template_length_out_of_range_pct: template_length_histogram.get_ignored_count() as f64
                / self.read_count_metrics.total as f64
                * 100.0,
        };

        self.summary_metrics = Some(summary)
    }
}
