//! Utilities for the Features quality control facet.

use std::str::FromStr;

/// Strand orientation (forward or reverse). Useful when parsing a GFF or other
/// features files.
#[derive(Clone, PartialEq, Eq)]
pub enum Strand {
    /// Features that fall on the _forward_ strand.
    Forward,

    /// Features that fall of the _reverse_ strand.
    Reverse,
}

/// Error that occurs when a strand cannot be parsed from a `&str`.
#[derive(Debug)]
pub struct StrandParseError(String);

impl std::fmt::Display for StrandParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "attempted to parse strand from value: {}", self.0)
    }
}

impl std::error::Error for StrandParseError {}

impl FromStr for Strand {
    type Err = StrandParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "+" => Ok(Strand::Forward),
            "-" => Ok(Strand::Reverse),
            _ => Err(StrandParseError(String::from(s))),
        }
    }
}

/// Utility struct that contains both a feature name and the strand that the
/// feature is contained on. This is used when building the
/// [`rust_lapper::Lapper`] interval lookup.
#[derive(Clone, PartialEq, Eq)]
pub struct FeatureNameStrand {
    name: String,
    strand: Strand,
}

impl FeatureNameStrand {
    /// Creates a new [`FeatureNameStrand`].
    pub fn new(name: String, strand: Strand) -> Self {
        FeatureNameStrand { name, strand }
    }

    /// Get a reference to the feature name strand's name.
    pub fn name(&self) -> &str {
        self.name.as_ref()
    }

    /// Get a reference to the feature name strand's strand.
    pub fn strand(&self) -> &Strand {
        &self.strand
    }
}
