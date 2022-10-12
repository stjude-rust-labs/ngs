//! Functionality related to the `ngs qc` subcommand.

use std::{path::PathBuf, rc::Rc};

use noodles::sam::{
    self,
    header::record::value::{map::ReferenceSequence, Map},
    Header,
};
use sam::alignment::Record;

use crate::utils::genome::ReferenceGenome;

use self::{
    record_based::{
        features::{FeatureNames, GenomicFeaturesFacet},
        gc_content::GCContentFacet,
        general::GeneralMetricsFacet,
        quality_scores::QualityScoreFacet,
        template_length::TemplateLengthFacet,
    },
    sequence_based::{coverage::CoverageFacet, edits::EditsFacet},
};

pub mod command;
pub mod record_based;
pub mod results;
pub mod sequence_based;

//==============================================//
// Dynamic allocation of quality control facets //
//==============================================//

/// Dynamically compiles the record-based quality control facets that should be
/// run for this invocation of the command line tool.
pub fn get_record_based_qc_facets<'a>(
    features_gff: Option<PathBuf>,
    feature_names: &'a FeatureNames,
    header: &'a Header,
    reference_genome: Rc<Box<dyn ReferenceGenome>>,
) -> anyhow::Result<Vec<Box<dyn RecordBasedQualityControlFacet + 'a>>> {
    // Default facets that are loaded within the qc subcommand.
    let mut facets: Vec<Box<dyn RecordBasedQualityControlFacet>> = vec![
        Box::new(GeneralMetricsFacet::default()),
        Box::new(TemplateLengthFacet::with_capacity(1024)),
        Box::new(GCContentFacet::default()),
        Box::new(QualityScoreFacet::default()),
    ];

    // Optionally load the Genomic Features facet if the GFF file is provided.
    if let Some(s) = features_gff {
        facets.push(Box::new(GenomicFeaturesFacet::try_from(
            s,
            feature_names,
            header,
            reference_genome,
        )?));
    }

    Ok(facets)
}

/// Dynamically compiles the sequence-based quality control facets that should
/// be run for this invocation of the command line tool.
pub fn get_sequence_based_qc_facets<'a>(
    reference_fasta: Option<PathBuf>,
    reference_genome: Rc<Box<dyn ReferenceGenome>>,
) -> anyhow::Result<Vec<Box<dyn SequenceBasedQualityControlFacet + 'a>>> {
    // Default facets that are loaded within the qc subcommand.
    let mut facets: Vec<Box<dyn SequenceBasedQualityControlFacet>> =
        vec![Box::new(CoverageFacet::new(reference_genome))];

    // Optionally load the Edits facet if a reference FASTA is provided.
    if let Some(fasta) = reference_fasta {
        facets.push(Box::new(EditsFacet::try_from(fasta)?))
    }

    Ok(facets)
}

//====================//
// Computational Load //
//====================//

/// An indicator of how computationally heavy an analysis is.
#[derive(Debug)]
pub enum ComputationalLoad {
    /// Less than 1 second per million records processed.
    Light,

    /// Between 1 and 5 seconds per million records processed.
    Moderate,

    /// More than 5 seconds per million records processed.
    Heavy,
}

//====================================//
// Record-based Quality Control Facet //
//====================================//

/// A struct is [`RecordBasedQualityControlFacet`] if it is a record-based
/// quality control facet.
pub trait RecordBasedQualityControlFacet {
    //================//
    // Static methods //
    //================//

    /// Name of the record-based quality control facet.
    fn name(&self) -> &'static str;

    /// Computational load of the record-based quality control facet.
    fn computational_load(&self) -> ComputationalLoad;

    //===================//
    // Lifecycle methods //
    //===================//

    /// Processes a record within this quality control facet.
    fn process(&mut self, record: &Record) -> anyhow::Result<()>;

    /// Summarizes the results of the quality control facet once all records
    /// have been processed.
    fn summarize(&mut self) -> anyhow::Result<()>;

    /// Adds the results of this quality control facet to the global
    /// [`results::Results`] object for writing to a file.
    fn aggregate(&self, results: &mut results::Results);
}

//======================================//
// Sequence-based Quality Control Facet //
//======================================//

/// A struct is [`SequenceBasedQualityControlFacet`] if it is a sequence-based
/// quality control facet.
pub trait SequenceBasedQualityControlFacet {
    //================//
    // Static methods //
    //================//

    /// Name of the sequence-based quality control facet.
    fn name(&self) -> &'static str;

    /// Computational load of the sequence-based quality control facet.
    fn computational_load(&self) -> ComputationalLoad;

    //===================//
    // Lifecycle methods //
    //===================//

    /// Informs on whether a sequence is supported by this quality-control
    /// facet.
    fn supports_sequence_name(&self, name: &str) -> bool;

    /// Sets up a quality control facet for a given sequence.
    fn setup(&mut self, sequence: &Map<ReferenceSequence>) -> anyhow::Result<()>;

    /// Processes a sequence for a quality control facet.
    fn process<'b>(
        &mut self,
        seq: &'b Map<ReferenceSequence>,
        record: &Record,
    ) -> anyhow::Result<()>;

    /// Tears down any machinery that was built up for this sequence within the
    /// quality control facet.
    fn teardown(&mut self, sequence: &Map<ReferenceSequence>) -> anyhow::Result<()>;

    /// Adds the results of this quality control facet to the global
    /// [`results::Results`] object for writing to a file.
    fn aggregate(&mut self, results: &mut results::Results);
}
