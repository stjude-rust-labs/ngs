//! Functionality related to the `ngs qc` subcommand.

use std::{num::NonZeroUsize, path::PathBuf, rc::Rc};

use anyhow::bail;
use itertools::Itertools;
use noodles::sam;
use noodles::sam::header::record::value::map::ReferenceSequence;
use noodles::sam::header::record::value::Map;
use noodles::sam::Header;
use sam::alignment::Record;

use crate::utils::genome::ReferenceGenome;

use self::record_based::features::FeatureNames;
use self::record_based::features::GenomicFeaturesFacet;
use self::record_based::gc_content::GCContentFacet;
use self::record_based::general::GeneralMetricsFacet;
use self::record_based::quality_scores::QualityScoreFacet;
use self::record_based::template_length::TemplateLengthFacet;
use self::sequence_based::coverage::CoverageFacet;
use self::sequence_based::edits::EditsFacet;

pub mod command;
pub mod record_based;
pub mod results;
pub mod sequence_based;

//==============================================//
// Dynamic allocation of quality control facets //
//==============================================//

type RecordBasedQualityControlFacetBoxedVec<'a> = Vec<Box<dyn RecordBasedQualityControlFacet + 'a>>;
type SequenceBasedQualityControlFacetBoxedVec<'a> =
    Vec<Box<dyn SequenceBasedQualityControlFacet + 'a>>;

/// Dynamically compiles the all of the quality control facets that should be
/// run for this invocation of the command line tool.
///
/// This method starts by defining the base set of QC facets that will be run
/// based on the arguments provided on the command line. Next, filtering is done
/// based on the arguments provided on the command line.
pub fn get_qc_facets<'a>(
    features_gff: Option<PathBuf>,
    feature_names: Option<&'a FeatureNames>,
    header: Option<&'a Header>,
    reference_fasta: Option<PathBuf>,
    reference_genome: Rc<Box<dyn ReferenceGenome>>,
    only_facet: Option<String>,
    vafs_file_path: Option<PathBuf>,
) -> anyhow::Result<(
    RecordBasedQualityControlFacetBoxedVec<'_>,
    SequenceBasedQualityControlFacetBoxedVec<'_>,
)> {
    // (1) Define the full list of facets that are supported for the
    // record-based quality control facets.

    // Default facets that are loaded within the qc subcommand.
    let mut record_based_facets: Vec<Box<dyn RecordBasedQualityControlFacet>> = vec![
        Box::new(GeneralMetricsFacet::default()),
        Box::new(TemplateLengthFacet::with_capacity(1024)),
        Box::new(GCContentFacet::default()),
        Box::new(QualityScoreFacet::default()),
    ];

    // Optionally load the Genomic Features facet if the GFF file is provided.
    if let Some(features_src) = features_gff {
        if let Some(feature_names) = feature_names {
            if let Some(header) = header {
                record_based_facets.push(Box::new(GenomicFeaturesFacet::try_from(
                    features_src,
                    feature_names,
                    header,
                    Rc::clone(&reference_genome),
                )?));
            }
        }
    }

    // (2) Define the full list of facets that are supported for the
    // sequence-based quality control facets.

    // Default facets that are loaded within the qc subcommand.
    let mut sequence_based_facets: Vec<Box<dyn SequenceBasedQualityControlFacet>> =
        vec![Box::new(CoverageFacet::new(
            Rc::clone(&reference_genome),
            NonZeroUsize::new(50_000).unwrap(),
        ))];

    // Optionally load the Edits facet if a reference FASTA is provided.
    if let Some(fasta) = reference_fasta {
        sequence_based_facets.push(Box::new(EditsFacet::try_from(&fasta, vafs_file_path)?));
    }

    // (3) If `only_facet` is provided, we need to (a) filter out all of the
    // quality control facets except the one that is provided, (b) error out if
    // no quality control facets match the provided argument, and (c) return
    // with the limited list otherwise.

    if let Some(only) = only_facet {
        let record_based_filtered = record_based_facets
            .into_iter()
            .filter(|x| x.name().eq_ignore_ascii_case(&only))
            .collect_vec();

        let sequence_based_filtered = sequence_based_facets
            .into_iter()
            .filter(|x| x.name().eq_ignore_ascii_case(&only))
            .collect_vec();

        let selected_facets_count = record_based_filtered.len() + sequence_based_filtered.len();

        match selected_facets_count {
            0 => bail!("No facets matched the specified `--only` flag: {}", only),
            1 => return Ok((record_based_filtered, sequence_based_filtered)),
            _ => bail!(
                "Too many facets matched the specified `--only` flag: {}. This is a \
                very strange error, and it should be reported on the Github issues page.",
                only
            ),
        }
    }

    Ok((record_based_facets, sequence_based_facets))
}

//====================//
// Computational Load //
//====================//

/// An indicator of how computationally heavy an analysis is.
#[derive(Debug)]
pub enum ComputationalLoad {
    /// Less than 1 second per million records processed.
    Light,

    /// Between 1 and 5 seconds per million records processed.
    Moderate,

    /// More than 5 seconds per million records processed.
    Heavy,
}

//====================================//
// Record-based Quality Control Facet //
//====================================//

/// A struct is [`RecordBasedQualityControlFacet`] if it is a record-based
/// quality control facet.
pub trait RecordBasedQualityControlFacet {
    //================//
    // Static methods //
    //================//

    /// Name of the record-based quality control facet.
    fn name(&self) -> &'static str;

    /// Computational load of the record-based quality control facet.
    fn computational_load(&self) -> ComputationalLoad;

    //===================//
    // Lifecycle methods //
    //===================//

    /// Processes a record within this quality control facet.
    fn process(&mut self, record: &Record) -> anyhow::Result<()>;

    /// Summarizes the results of the quality control facet once all records
    /// have been processed.
    fn summarize(&mut self) -> anyhow::Result<()>;

    /// Adds the results of this quality control facet to the global
    /// [`results::Results`] object for writing to a file.
    fn aggregate(&self, results: &mut results::Results);
}

//======================================//
// Sequence-based Quality Control Facet //
//======================================//

/// A struct is [`SequenceBasedQualityControlFacet`] if it is a sequence-based
/// quality control facet.
pub trait SequenceBasedQualityControlFacet {
    //================//
    // Static methods //
    //================//

    /// Name of the sequence-based quality control facet.
    fn name(&self) -> &'static str;

    /// Computational load of the sequence-based quality control facet.
    fn computational_load(&self) -> ComputationalLoad;

    //===================//
    // Lifecycle methods //
    //===================//

    /// Informs on whether a sequence is supported by this quality-control
    /// facet.
    fn supports_sequence_name(&self, name: &str) -> bool;

    /// Sets up a quality control facet for a given sequence.
    fn setup(&mut self, sequence: &Map<ReferenceSequence>) -> anyhow::Result<()>;

    /// Processes a sequence for a quality control facet.
    fn process<'b>(
        &mut self,
        seq: &'b Map<ReferenceSequence>,
        record: &Record,
    ) -> anyhow::Result<()>;

    /// Tears down any machinery that was built up for this sequence within the
    /// quality control facet.
    fn teardown(&mut self, sequence: &Map<ReferenceSequence>) -> anyhow::Result<()>;

    /// Adds the results of this quality control facet to the global
    /// [`results::Results`] object for writing to a file.
    fn aggregate(&mut self, results: &mut results::Results);
}

#[cfg(test)]
mod tests {

    use crate::utils::genome::get_reference_genome;

    use super::*;

    #[test]
    pub fn it_returns_the_correct_number_of_facets_by_default() {
        let (record_based, sequence_based) = get_qc_facets(
            None,
            None,
            None,
            None,
            Rc::new(get_reference_genome("GRCh38_no_alt_AnalysisSet").unwrap()),
            None,
            None,
        )
        .unwrap();

        assert_eq!(record_based.len(), 4);
        assert_eq!(sequence_based.len(), 1);
    }

    #[test]
    pub fn it_returns_the_correct_number_of_facets_when_only_is_specified() {
        let (record_based, sequence_based) = get_qc_facets(
            None,
            None,
            None,
            None,
            Rc::new(get_reference_genome("GRCh38_no_alt_AnalysisSet").unwrap()),
            Some(String::from("GC Content")),
            None,
        )
        .unwrap();

        assert_eq!(record_based.len(), 1);
        assert_eq!(sequence_based.len(), 0);
        assert!(record_based.get(0).unwrap().name() == "GC Content");
    }
}
