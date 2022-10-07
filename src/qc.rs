use noodles_sam as sam;
use sam::{alignment::Record, header::ReferenceSequence};

pub mod command;
pub mod record_based;
pub mod results;
pub mod sequence_based;

#[derive(Debug)]
#[allow(dead_code)]
pub enum ComputationalLoad {
    Light,
    Moderate,
    Heavy,
}

pub trait RecordBasedQualityCheckFacet {
    // Static methods
    fn name(&self) -> &'static str;
    fn computational_load(&self) -> ComputationalLoad;

    // Lifecycle methods
    fn process(&mut self, record: &Record) -> anyhow::Result<()>;
    fn summarize(&mut self) -> anyhow::Result<()>;
    fn aggregate(&self, results: &mut results::Results);
}

pub trait SequenceBasedQualityCheckFacet<'a> {
    // Static methods
    fn name(&self) -> &'static str;
    fn computational_load(&self) -> ComputationalLoad;

    // Lifecycle methods
    fn supports_sequence_name(&self, name: &str) -> bool;
    fn setup(&mut self, sequence: &ReferenceSequence) -> anyhow::Result<()>;
    fn process<'b>(
        &mut self,
        seq: &'b ReferenceSequence,
        record: &sam::alignment::Record,
    ) -> anyhow::Result<()>
    where
        'b: 'a;
    fn teardown(&mut self, sequence: &ReferenceSequence) -> anyhow::Result<()>;
    fn aggregate(&mut self, results: &mut results::Results);
}
