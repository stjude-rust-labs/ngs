#[derive(Clone, PartialEq, Eq)]
pub struct FeatureNameStrand {
    name: String,
    strand: String,
}

impl FeatureNameStrand {
    pub fn new(name: String, strand: String) -> Self {
        FeatureNameStrand { name, strand }
    }

    /// Get a reference to the feature name strand's name.
    pub fn name(&self) -> &str {
        self.name.as_ref()
    }

    /// Get a reference to the feature name strand's strand.
    #[allow(dead_code)]
    pub fn strand(&self) -> &str {
        self.strand.as_ref()
    }
}
