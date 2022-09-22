use core::fmt;

use noodles_bam::lazy::Record;

pub mod features;
pub mod gc_content;
pub mod general;
pub mod results;
pub mod template_length;

#[derive(Debug)]
pub struct Error {
    pub message: String,
}

impl Error {
    pub fn new<I>(message: I) -> Self
    where
        I: Into<String>,
    {
        Error {
            message: message.into(),
        }
    }
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.message)
    }
}

#[derive(Debug)]
#[allow(dead_code)]
pub enum ComputationalLoad {
    Light,
    Moderate,
    Heavy,
}

pub trait QualityCheckFacet {
    fn name(&self) -> &'static str;
    fn computational_load(&self) -> ComputationalLoad;
    fn process(&mut self, record: &Record) -> Result<(), Error>;
    fn summarize(&mut self) -> Result<(), Error>;
    fn aggregate_results(&self, results: &mut results::Results);
}
