use std::{
    fs::{self, File},
    io::{self, Write},
    path::{Path, PathBuf},
};

use serde::{Deserialize, Serialize};

use super::{
    record_based::{features, gc_content, general, quality_scores, template_length},
    sequence_based::{coverage, edits},
};

#[derive(Default, Serialize, Deserialize)]
pub struct Results {
    general: Option<general::metrics::GeneralMetrics>,
    features: Option<features::Metrics>,
    gc_content: Option<gc_content::GCContentMetrics>,
    template_length: Option<template_length::TemplateLengthFacet>,
    coverage: Option<coverage::CoverageMetrics>,
    edits: Option<edits::EditMetrics>,
    quality_scores: Option<quality_scores::QualityScoreFacet>,
}

impl Results {
    pub fn write(&self, output_prefix: String, directory: &Path) -> Result<(), io::Error> {
        let features_filename = output_prefix + ".results.json";
        let mut features_filepath = PathBuf::from(directory);
        features_filepath.push(features_filename);

        let mut file = File::create(features_filepath)?;
        let output = serde_json::to_string_pretty(&self).unwrap();
        file.write_all(output.as_bytes())?;

        Ok(())
    }

    pub fn read(filepath: impl AsRef<Path>) -> anyhow::Result<Results> {
        let path = filepath.as_ref();
        let contents = fs::read_to_string(path)?;
        Ok(serde_json::from_str(&contents)?)
    }

    //=========//
    // General //
    //=========//

    #[allow(dead_code)]
    pub fn general(&self) -> Option<&general::metrics::GeneralMetrics> {
        self.general.as_ref()
    }

    pub fn set_general(&mut self, general: general::GeneralMetricsFacet) {
        self.general = Some(general.metrics);
    }

    //==========//
    // Features //
    //==========//

    #[allow(dead_code)]
    pub fn features(&self) -> Option<&features::Metrics> {
        self.features.as_ref()
    }

    pub fn set_features(&mut self, features: features::Metrics) {
        self.features = Some(features);
    }

    //============//
    // GC Content //
    //============//

    #[allow(dead_code)]
    pub fn gc_content(&self) -> Option<&gc_content::GCContentMetrics> {
        self.gc_content.as_ref()
    }

    pub fn set_gc_content(&mut self, gc_content: gc_content::GCContentMetrics) {
        self.gc_content = Some(gc_content);
    }

    //=================//
    // Template Length //
    //=================//

    #[allow(dead_code)]
    pub fn template_length(&self) -> Option<&template_length::TemplateLengthFacet> {
        self.template_length.as_ref()
    }

    pub fn set_template_length(&mut self, template_length: template_length::TemplateLengthFacet) {
        self.template_length = Some(template_length);
    }

    //==========//
    // Coverage //
    //==========//

    #[allow(dead_code)]
    pub fn coverage(&self) -> Option<&coverage::CoverageMetrics> {
        self.coverage.as_ref()
    }

    pub fn set_coverage(&mut self, coverage: coverage::CoverageMetrics) {
        self.coverage = Some(coverage);
    }

    //=======//
    // Edits //
    //=======//

    #[allow(dead_code)]
    pub fn edits(&self) -> Option<&edits::EditMetrics> {
        self.edits.as_ref()
    }

    pub fn set_edits(&mut self, edits: edits::EditMetrics) {
        self.edits = Some(edits);
    }

    //================//
    // Quality Scores //
    //================//

    pub fn quality_scores(&self) -> Option<&quality_scores::QualityScoreFacet> {
        self.quality_scores.as_ref()
    }

    pub fn set_quality_scores(&mut self, quality_scores: quality_scores::QualityScoreFacet) {
        self.quality_scores = Some(quality_scores);
    }
}
