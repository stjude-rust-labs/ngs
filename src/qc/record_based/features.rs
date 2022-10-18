//! Functionality related to the Features quality control facet.

use std::collections::HashMap;
use std::path::PathBuf;
use std::rc::Rc;

use anyhow::bail;
use anyhow::Context;
use noodles::sam;
use rust_lapper::Interval;
use rust_lapper::Lapper;
use sam::alignment::Record;
use sam::Header;
use tracing::debug;

pub mod metrics;
pub mod utils;

use crate::qc::record_based::features::utils::Strand;
use crate::qc::results;
use crate::qc::ComputationalLoad;
use crate::qc::RecordBasedQualityControlFacet;
use crate::utils::formats;
use crate::utils::genome::get_primary_assembly;
use crate::utils::genome::ReferenceGenome;

pub use self::metrics::Metrics;
pub use self::metrics::SummaryMetrics;
pub use self::utils::FeatureNameStrand;

//=================//
// Utility structs //
//=================//

/// A utility struct for passing feature name arguments from the command line
/// around more easily.
pub struct FeatureNames {
    /// Name of the 5' UTR region feature for the gene model used.
    pub five_prime_utr_feature_name: String,

    /// Name of the 3' UTR region feature for the gene model used.
    pub three_prime_utr_feature_name: String,

    /// Name of the coding sequence region feature for the gene model used.
    pub coding_sequence_feature_name: String,

    /// Name of the exon region feature for the gene model used.
    pub exon_feature_name: String,

    /// Name of the gene region feature for the gene model used.
    pub gene_feature_name: String,
}

impl FeatureNames {
    /// Creates a new [`FeatureNames`].
    pub fn new<I>(
        five_prime_utr_feature_name: I,
        three_prime_utr_feature_name: I,
        coding_sequence_feature_name: I,
        exon_feature_name: I,
        gene_feature_name: I,
    ) -> Self
    where
        I: Into<String>,
    {
        FeatureNames {
            five_prime_utr_feature_name: five_prime_utr_feature_name.into(),
            three_prime_utr_feature_name: three_prime_utr_feature_name.into(),
            coding_sequence_feature_name: coding_sequence_feature_name.into(),
            exon_feature_name: exon_feature_name.into(),
            gene_feature_name: gene_feature_name.into(),
        }
    }
}

//========================//
// Genomic Features Facet //
//========================//

/// Main struct for the Features quality control facet.
pub struct GenomicFeaturesFacet<'a> {
    /// Store of the cached exonic translation regions.
    pub exonic_translation_regions: HashMap<String, Lapper<usize, FeatureNameStrand>>,

    /// Store of the cached gene regions.
    pub gene_regions: HashMap<String, Lapper<usize, FeatureNameStrand>>,

    /// Feature names that correspond to the respective features in the gene
    /// model. These are passed in on the command line.
    pub feature_names: &'a FeatureNames,

    /// The SAM header. Useful for looking up sequence ids for records during
    /// processing.
    pub header: &'a Header,

    /// The main metric counting struct.
    pub metrics: Metrics,

    /// Cached set of primary chromosomes (for ignoring records aligned to
    /// non-primary sequences).
    pub primary_chromosome_names: Vec<String>,
}

impl<'a> RecordBasedQualityControlFacet for GenomicFeaturesFacet<'a> {
    fn name(&self) -> &'static str {
        "Genomic Features"
    }

    fn computational_load(&self) -> ComputationalLoad {
        ComputationalLoad::Moderate
    }

    fn process(&mut self, record: &Record) -> anyhow::Result<()> {
        // (1) Parse the read name.
        let read_name = match record.read_name() {
            Some(name) => name,
            _ => bail!("Could not parse read name"),
        };

        // (2) Parse the flags so we can see if the read is mapped.
        let flags = record.flags();

        // (3) If the read is unmapped, just return—no need to throw an error.
        if flags.is_unmapped() {
            self.metrics.records.ignored_flags += 1;
            return Ok(());
        }

        // (4) Parse the reference sequence id from the record.
        let id = match record.reference_sequence_id() {
            Some(id) => id,
            _ => {
                bail!(
                    "Could not parse reference sequence id for read: {}",
                    read_name
                )
            }
        };

        // (5) Map the parsed reference sequence id to a reference sequence name.
        let seq_name = match self
            .header
            .reference_sequences()
            .get_index(id)
            .map(|(_, rs)| Some(rs.name().as_str()))
        {
            Some(Some(name)) => name,
            _ => {
                bail!(
                    "Could not map reference sequence id to header for read: {}",
                    read_name
                )
            }
        };

        if !self
            .primary_chromosome_names
            .iter()
            .any(|s: &String| s == seq_name)
        {
            self.metrics.records.ignored_nonprimary_chromosome += 1;
            return Ok(());
        }

        // (6) Calculate the start and end position of this read. This will
        // later be used for lookup within our feature map.
        let start = match record.alignment_start() {
            Some(s) => usize::from(s),
            _ => bail!("Could not parse record's start position."),
        };

        let cigar = record.cigar();

        let end = start + cigar.alignment_span();

        // (7) Tally up which features this record contributes to.
        let mut counted_as_five_prime_utr = false;
        let mut counted_as_three_prime_utr = false;
        let mut counted_as_coding_sequence = false;

        // (7a) Tally up exonic translations.
        if let Some(utrs) = self.exonic_translation_regions.get(seq_name) {
            for utr in utrs.find(start, end + 1) {
                let f = &utr.val;

                if !counted_as_five_prime_utr
                    && f.name() == self.feature_names.five_prime_utr_feature_name
                {
                    counted_as_five_prime_utr = true;
                    self.metrics.exonic_translation_regions.utr_five_prime_count += 1;
                } else if !counted_as_three_prime_utr
                    && f.name() == self.feature_names.three_prime_utr_feature_name
                {
                    counted_as_three_prime_utr = true;
                    self.metrics
                        .exonic_translation_regions
                        .utr_three_prime_count += 1;
                } else if !counted_as_coding_sequence
                    && f.name() == self.feature_names.coding_sequence_feature_name
                {
                    counted_as_coding_sequence = true;
                    self.metrics
                        .exonic_translation_regions
                        .coding_sequence_count += 1;
                }
            }
        }

        // (7b) Tally up gene regions.
        if let Some(genics) = self.gene_regions.get(seq_name) {
            let mut has_gene = false;
            let mut has_exon = false;
            for gene in genics.find(start, end + 1) {
                let f = &gene.val;

                if f.name() == self.feature_names.gene_feature_name {
                    has_gene = true;
                } else if f.name() == self.feature_names.exon_feature_name {
                    has_exon = true;
                }

                if has_gene && has_exon {
                    break;
                }
            }

            if has_gene {
                if has_exon {
                    self.metrics.gene_regions.exonic_count += 1;
                } else {
                    self.metrics.gene_regions.intronic_count += 1;
                }
            } else {
                self.metrics.gene_regions.intergenic_count += 1;
            }
        }

        self.metrics.records.processed += 1;
        Ok(())
    }

    fn summarize(&mut self) -> anyhow::Result<()> {
        self.metrics.summary = Some(SummaryMetrics {
            ignored_flags_pct: (self.metrics.records.ignored_flags as f64
                / (self.metrics.records.ignored_flags
                    + self.metrics.records.ignored_nonprimary_chromosome
                    + self.metrics.records.processed) as f64)
                * 100.0,
            ignored_nonprimary_chromosome_pct: (self.metrics.records.ignored_nonprimary_chromosome
                as f64
                / (self.metrics.records.ignored_flags
                    + self.metrics.records.ignored_nonprimary_chromosome
                    + self.metrics.records.processed) as f64)
                * 100.0,
        });

        Ok(())
    }

    fn aggregate(&self, results: &mut results::Results) {
        results.features = Some(self.metrics.clone())
    }
}

impl<'a> GenomicFeaturesFacet<'a> {
    /// Tries to create a [`GenomicFeaturesFacet`] from a set of provided
    /// arguments. May fail if there are issues opening the GFF file.
    pub fn try_from(
        src: PathBuf,
        feature_names: &'a FeatureNames,
        header: &'a Header,
        reference_genome: Rc<Box<dyn ReferenceGenome>>,
    ) -> anyhow::Result<Self> {
        let mut gff = formats::gff::open(&src)
            .with_context(|| format!("Could not open GFF: {}", src.display()))?;

        let mut exonic_translations = HashMap::new();
        let mut gene_regions = HashMap::new();

        debug!("Reading all records in GFF.");
        let mut gff_records = Vec::new();
        for result in gff.records() {
            let record = result.unwrap();
            gff_records.push(record);
        }

        debug!("Tabulating GFF features.");
        let primary_assembly_sequence_names: Vec<String> = get_primary_assembly(reference_genome)
            .iter()
            .map(|s| String::from(s.name()))
            .collect();

        for parent_seq_name in primary_assembly_sequence_names.iter() {
            let mut utr_features: Vec<Interval<usize, FeatureNameStrand>> = Vec::new();
            let mut gene_region_features: Vec<Interval<usize, FeatureNameStrand>> = Vec::new();

            for record in &gff_records {
                let seq_name = record.reference_sequence_name();
                if parent_seq_name == seq_name {
                    let ty = record.ty();
                    let strand = record.strand();
                    let start: usize = record.start().into();
                    let stop: usize = record.end().into();

                    let feature_name_strand = FeatureNameStrand::new(
                        ty.to_string(),
                        strand.to_string().parse::<Strand>()?,
                    );

                    let interval = Interval {
                        start,
                        stop,
                        val: feature_name_strand,
                    };

                    if ty == feature_names.five_prime_utr_feature_name
                        || ty == feature_names.three_prime_utr_feature_name
                        || ty == feature_names.coding_sequence_feature_name
                    {
                        utr_features.push(interval);
                    } else if ty == feature_names.exon_feature_name
                        || ty == feature_names.gene_feature_name
                    {
                        gene_region_features.push(interval);
                    }
                }
            }

            debug!(
                "{} has {} gene region features and {} exonic translation features.",
                parent_seq_name,
                gene_region_features.len(),
                utr_features.len()
            );

            exonic_translations.insert(parent_seq_name.to_string(), Lapper::new(utr_features));
            gene_regions.insert(
                parent_seq_name.to_string(),
                Lapper::new(gene_region_features),
            );
        }

        debug!("Finalizing GFF features lookup.");

        Ok(Self {
            exonic_translation_regions: exonic_translations,
            gene_regions,
            feature_names,
            header,
            metrics: Metrics::default(),
            primary_chromosome_names: primary_assembly_sequence_names,
        })
    }
}
