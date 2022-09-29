use noodles_sam as sam;
use sam::{alignment::Record, header::ReferenceSequence};

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
    fn name(&self) -> &'static str;
    fn computational_load(&self) -> ComputationalLoad;
    fn process(&mut self, record: &Record) -> anyhow::Result<()>;
    fn summarize(&mut self) -> anyhow::Result<()>;
    fn aggregate_results(&self, results: &mut results::Results);
}

pub trait SequenceBasedQualityCheckFacet<'a> {
    fn name(&self) -> &'static str;
    fn computational_load(&self) -> ComputationalLoad;
    fn supports_sequence_name(&self, name: &str) -> bool;
    fn setup_sequence(&mut self, seq: &ReferenceSequence) -> anyhow::Result<()>;
    fn process_record<'b>(
        &mut self,
        seq: &'b ReferenceSequence,
        record: &sam::alignment::Record,
    ) -> anyhow::Result<()>
    where
        'b: 'a;
    fn teardown_sequence(&mut self, seq: &ReferenceSequence) -> anyhow::Result<()>;
    fn aggregate_results(&mut self, results: &mut results::Results);
}
