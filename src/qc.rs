//! Functionality related to the `ngs qc` subcommand.

use noodles::sam::{
    self,
    header::record::value::{map::ReferenceSequence, Map},
};
use sam::alignment::Record;

pub mod command;
pub mod record_based;
pub mod results;
pub mod sequence_based;

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

/// A struct is [`RecordBasedQualityCheckFacet`] if it is a record-based
/// quality control facet.
pub trait RecordBasedQualityCheckFacet {
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

/// A struct is [`SequenceBasedQualityCheckFacet`] if it is a sequence-based
/// quality control facet.
pub trait SequenceBasedQualityCheckFacet {
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
