//! Utilities related to the parsing of arguments.

use std::fmt::Display;

use noodles::{bgzf::writer::CompressionLevel, sam::record::MappingQuality};
use tracing::debug;

//===================//
// Number of Records //
//===================//

/// Utility enum to designate whether we are reviewing all records in the file
/// or just some of them.
pub enum NumberOfRecords {
    /// Designates that we should review _all_ of the records in the file.
    All,

    /// Designates that we should review _some_ of the records in the file. The
    /// exact count of records is stored in the `usize`.
    Some(usize),
}

impl From<Option<usize>> for NumberOfRecords {
    fn from(num_records: Option<usize>) -> Self {
        match num_records {
            Some(n) => {
                debug!("Reading a maximum of {} records.", n);
                NumberOfRecords::Some(n)
            }
            None => {
                debug!("Reading all available records.");
                NumberOfRecords::All
            }
        }
    }
}

//======================//
// Compression Strategy //
//======================//

/// An enum representing the compression strategy to follow.
#[derive(clap::ValueEnum, Clone, Debug)]
pub enum CompressionStrategy {
    /// Compress the file as much as possible (maximum gzip compression).
    Best,

    /// Balance the compression level and the speed of the compression process.
    Balanced,

    /// Compress the file as quickly as possible (minimum gzip compression without
    /// turning off compression altogether).
    Fastest,
}

impl Display for CompressionStrategy {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Balanced => write!(f, "balanced"),
            Self::Best => write!(f, "best"),
            Self::Fastest => write!(f, "fastest"),
        }
    }
}

impl From<CompressionStrategy> for CompressionLevel {
    fn from(strategy: CompressionStrategy) -> Self {
        match strategy {
            CompressionStrategy::Balanced => CompressionLevel::try_from(6)
                .expect("can create a compression level with value of 6"),
            CompressionStrategy::Best => CompressionLevel::best(),
            CompressionStrategy::Fastest => CompressionLevel::fast(),
        }
    }
}

//=============//
// Arg Parsers //
//=============//

/// Utility method to parse command line floats and ensure they are
/// within the range [MIN, MAX].
pub fn arg_in_range(arg: f32, range: std::ops::RangeInclusive<f32>) -> anyhow::Result<f32> {
    match range.contains(&arg) {
        true => Ok(arg),
        false => anyhow::bail!(
            "Value must be between {} and {}",
            range.start(),
            range.end()
        ),
    }
}

// TODO dead code, not used. Doesn't work as written.
/// Utility method to parse command line integers and ensure they are
/// within the range [0, 255) and return them as MappingQualities.
pub fn parse_min_mapq(s: &str) -> Result<Option<MappingQuality>, std::num::ParseIntError> {
    let value = s.parse()?;
    match value {
        0..=254 => Ok(Some(MappingQuality::new(value).unwrap())),
        255 => Ok(None),
    }
}
