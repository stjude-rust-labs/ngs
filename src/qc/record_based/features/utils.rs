//! Utilities for the Features quality control facet.

/// Utility struct that contains both a feature name and the strand that the
/// feature is contained on. This is used when building the
/// [`rust_lapper::Lapper`] interval lookup.
#[derive(Clone, PartialEq, Eq)]
pub struct FeatureNameStrand {
    name: String,
    strand: String,
}

impl FeatureNameStrand {
    /// Creates a new [`FeatureNameStrand`].
    pub fn new(name: String, strand: String) -> Self {
        FeatureNameStrand { name, strand }
    }

    /// Get a reference to the feature name strand's name.
    pub fn name(&self) -> &str {
        self.name.as_ref()
    }

    /// Get a reference to the feature name strand's strand.
    pub fn strand(&self) -> &str {
        self.strand.as_ref()
    }
}
