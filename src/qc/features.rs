use std::{
    collections::HashMap,
    fs::File,
    io::{self, Write},
    path::{Path, PathBuf},
};

use noodles_bam::lazy::Record;
use noodles_sam as sam;
use rust_lapper::{Interval, Lapper};
use sam::Header;
use tracing::debug;

use crate::{commands::qc::FeatureNames, formats, utils::genome::PRIMARY_CHROMOSOMES};

use self::{
    metrics::{Metrics, SummaryMetrics},
    name_strand::FeatureNameStrand,
};

use super::{ComputationalLoad, Error, QualityCheckFacet};

pub mod metrics;
pub mod name_strand;

pub struct GenomicFeaturesFacet<'a> {
    exonic_translation_regions: HashMap<String, Lapper<usize, FeatureNameStrand>>,
    gene_regions: HashMap<String, Lapper<usize, FeatureNameStrand>>,
    feature_names: &'a FeatureNames,
    header: &'a Header,
    metrics: Metrics,
}

impl<'a> QualityCheckFacet for GenomicFeaturesFacet<'a> {
    fn name(&self) -> &'static str {
        "Genomic Features Metrics"
    }

    fn computational_load(&self) -> ComputationalLoad {
        ComputationalLoad::Moderate
    }

    fn process(&mut self, record: &Record) -> Result<(), super::Error> {
        // (1) Parse the read name.
        let read_name = match record.read_name() {
            Ok(Some(name)) => name,
            _ => return Err(Error::new("Could not parse read name")),
        };

        // (2) Parse the flags so we can see if the read is mapped.
        let flags = match record.flags() {
            Ok(flags) => flags,
            Err(_) => {
                return Err(Error::new(format!(
                    "Could not parse flags for read: {}",
                    read_name
                )))
            }
        };

        // (3) If the read is unmapped, just return—no need to throw an error.
        if flags.is_unmapped() {
            self.metrics.records.ignored_flags += 1;
            return Ok(());
        }

        // (4) Parse the reference sequence id from the record.
        let id = match record.reference_sequence_id() {
            Ok(Some(id)) => id,
            _ => {
                return Err(Error::new(format!(
                    "Could not parse reference sequence id for read: {}",
                    read_name
                )))
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
                return Err(Error::new(format!(
                    "Could map reference sequence id to header for read: {}",
                    read_name
                )))
            }
        };

        if !PRIMARY_CHROMOSOMES.contains(&seq_name) {
            self.metrics.records.ignored_nonprimary_chromosome += 1;
            return Ok(());
        }

        // (6) Calculate the start and end position of this read. This will
        // later be used for lookup within our feature map.
        let start = match record.alignment_start() {
            Ok(Some(s)) => usize::from(s),
            _ => return Err(Error::new("Could not parse record's start position.")),
        };

        let cigar = sam::record::Cigar::try_from(record.cigar())
            .expect("Could not parse BAM's CIGAR string.");

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

    fn summarize(&mut self) -> Result<(), Error> {
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

    fn write(&self, output_prefix: String, directory: &Path) -> Result<(), io::Error> {
        let features_filename = output_prefix + ".features.json";
        let mut features_filepath = PathBuf::from(directory);
        features_filepath.push(features_filename);

        let mut file = File::create(features_filepath)?;
        let output = serde_json::to_string_pretty(&self.get_metrics()).unwrap();
        file.write_all(output.as_bytes())?;

        Ok(())
    }
}

impl<'a> GenomicFeaturesFacet<'a> {
    pub fn get_metrics(&self) -> &Metrics {
        &self.metrics
    }

    pub fn try_from(
        src: &str,
        feature_names: &'a FeatureNames,
        header: &'a Header,
    ) -> io::Result<Self> {
        let mut gff = formats::gff::open(src)?;

        let mut exonic_translations = HashMap::new();
        let mut gene_regions = HashMap::new();

        debug!("Reading all records in GFF.");
        let mut gff_records = Vec::new();
        for result in gff.records() {
            let record = result.unwrap();
            gff_records.push(record);
        }

        debug!("Tabulating GFF features.");
        for parent_seq_name in PRIMARY_CHROMOSOMES {
            let mut utr_features: Vec<Interval<usize, FeatureNameStrand>> = Vec::new();
            let mut gene_region_features: Vec<Interval<usize, FeatureNameStrand>> = Vec::new();

            for record in &gff_records {
                let seq_name = record.reference_sequence_name();
                if parent_seq_name == seq_name {
                    let ty = record.ty();
                    let strand = record.strand();
                    let start: usize = record.start().into();
                    let stop: usize = record.end().into();

                    let feature_name_strand =
                        FeatureNameStrand::new(ty.to_string(), strand.to_string());

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
            metrics: Metrics::new(),
        })
    }
}
