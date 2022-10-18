//! Functionality related to the aggregation of results across all quality
//! control facets.

use std::fs;
use std::fs::File;
use std::io;
use std::io::Write;
use std::path::Path;
use std::path::PathBuf;

use serde::Deserialize;
use serde::Serialize;

use super::record_based::features;
use super::record_based::gc_content;
use super::record_based::general;
use super::record_based::quality_scores;
use super::record_based::template_length;
use super::sequence_based::coverage;
use super::sequence_based::edits;

/// Main struct for collecting _all_ quality control facet results.
#[derive(Default, Serialize, Deserialize)]
pub struct Results {
    /// The quality control results from the General facet.
    pub general: Option<general::metrics::GeneralMetrics>,

    /// The quality control results from the Features facet.
    pub features: Option<features::Metrics>,

    /// The quality control results from the GC Content facet.
    pub gc_content: Option<gc_content::metrics::GCContentMetrics>,

    /// The quality control results from the Template Length facet.
    pub template_length: Option<template_length::TemplateLengthFacet>,

    /// The quality control results from the Quality Scores facet.
    pub quality_scores: Option<quality_scores::QualityScoreFacet>,

    /// The quality control results from the Coverage facet.
    pub coverage: Option<coverage::CoverageMetrics>,

    /// The quality control results from the Edits facet.
    pub edits: Option<edits::EditMetrics>,
}

impl Results {
    /// Attempts to write the [`Results`] struct to a file within the specified
    /// directory.
    pub fn write(&self, output_prefix: String, directory: &Path) -> Result<(), io::Error> {
        let features_filename = output_prefix + ".results.json";
        let mut features_filepath = PathBuf::from(directory);
        features_filepath.push(features_filename);

        let mut file = File::create(features_filepath)?;
        let output = serde_json::to_string_pretty(&self).unwrap();
        file.write_all(output.as_bytes())?;

        Ok(())
    }

    /// Attempts to read a [`Results`] struct from a file.
    pub fn read(filepath: impl AsRef<Path>) -> anyhow::Result<Results> {
        let path = filepath.as_ref();
        let contents = fs::read_to_string(path)?;
        Ok(serde_json::from_str(&contents)?)
    }
}
