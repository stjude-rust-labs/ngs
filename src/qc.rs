use core::fmt;
use std::{io, path::Path};

use noodles_bam::lazy::Record;

pub mod features;
pub mod gc_content;
pub mod metrics;
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
    fn default(&self) -> bool;
    fn process(&mut self, record: &Record) -> Result<(), Error>;
    fn summarize(&mut self) -> Result<(), Error>;
    fn write(&self, output_prefix: String, directory: &Path) -> Result<(), io::Error>;
}
