use std::{
    fs::File,
    io::{self, Write},
    path::{Path, PathBuf},
};

use serde::Serialize;

use super::{coverage, features, gc_content, general, template_length};

#[derive(Default, Serialize)]
pub struct Results {
    general: Option<general::GeneralMetricsFacet>,
    features: Option<features::Metrics>,
    gc_content: Option<gc_content::GCContentMetrics>,
    template_length: Option<template_length::TemplateLengthFacet>,
    coverage: Option<coverage::CoverageMetrics>,
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

    pub fn set_template_length(&mut self, template_length: template_length::TemplateLengthFacet) {
        self.template_length = Some(template_length);
    }

    pub fn set_gc_content(&mut self, gc_content: gc_content::GCContentMetrics) {
        self.gc_content = Some(gc_content);
    }

    pub fn set_general(&mut self, general: general::GeneralMetricsFacet) {
        self.general = Some(general);
    }

    pub fn set_features(&mut self, features: features::Metrics) {
        self.features = Some(features);
    }

    pub fn set_coverage(&mut self, coverage: coverage::CoverageMetrics) {
        self.coverage = Some(coverage);
    }
}
