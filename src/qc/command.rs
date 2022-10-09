//! Functionality related to the `ngs qc` command itself.

use std::{fs::File, path::PathBuf, rc::Rc};

use anyhow::{bail, Context};
use clap::Args;
use noodles::bam::{self as bam, bai};
use noodles::core::{Position, Region};
use noodles::sam::Header;
use num_format::{Locale, ToFormattedString};
use tracing::{debug, info};

use crate::{
    qc::results::Results,
    utils::{
        formats::sam::parse_header,
        genome::{get_all_sequences, get_reference_genome, ReferenceGenome},
    },
};

use super::{
    record_based::{
        features::{FeatureNames, GenomicFeaturesFacet},
        gc_content::GCContentFacet,
        general::GeneralMetricsFacet,
        quality_scores::QualityScoreFacet,
        template_length::TemplateLengthFacet,
    },
    sequence_based::{coverage::CoverageFacet, edits::EditsFacet},
    RecordBasedQualityControlFacet, SequenceBasedQualityControlFacet,
};

//====================================//
// Command line parsing utility types //
//====================================//

/// Utility enum to designate whether we are reviewing all records in the file
/// or just some of them. TODO: this may be generally useful, so we may want to
/// pull this out into its own utility module.
pub enum NumberOfRecords {
    /// Designates that we should review _all_ of the records in the file.
    All,

    /// Designates that we should review _some_ of the records in the file. The
    /// exact count of records is stored in the `usize`.
    Some(usize),
}

//========================//
// Command line arguments //
//========================//

/// Clap arguments for the `ngs qc` subcommand.
#[derive(Args)]
pub struct QcArgs {
    /// Source BAM file.
    #[arg(value_name = "BAM")]
    src: PathBuf,

    /// Supported reference genome used as the basis for analysis.
    reference_genome: String,

    /// Features GFF file (some metrics only supported if present).
    #[arg(short = 'f', long, value_name = "PATH")]
    features_gff: Option<PathBuf>,

    /// Number of records to process in the first pass.
    #[arg(short = 'n', long, value_name = "USIZE")]
    num_records: Option<usize>,

    /// Directory to output files to. Defaults to current working directory.
    #[arg(short = 'o', long, value_name = "PATH")]
    output_directory: Option<PathBuf>,

    /// Output prefix for the files that will be created. Defaults to the name
    /// of the file.
    #[arg(short = 'p', long, value_name = "STRING")]
    output_prefix: Option<String>,

    /// Reference FASTA file (some metrics only supported if present).
    #[arg(short = 'r', long, value_name = "PATH")]
    reference_fasta: Option<PathBuf>,

    /// Name of the feature that represents a five prime UTR region in the GFF
    /// file. Defaults to the respective GENCODE feature name.
    #[arg(long, value_name = "STRING", default_value = "five_prime_UTR")]
    five_prime_utr_feature_name: String,

    /// Name of the feature that represents a three prime UTR region in the GFF
    /// file. Defaults to the respective GENCODE feature name.
    #[arg(long, value_name = "STRING", default_value = "three_prime_UTR")]
    three_prime_utr_feature_name: String,

    /// Name of the feature that represents a coding sequence region in the GFF
    /// file. Defaults to the respective GENCODE feature name.
    #[arg(long, value_name = "STRING", default_value = "CDS")]
    coding_sequence_feature_name: String,

    /// Name of the feature that represents an exonic region in the GFF file.
    /// Defaults to the repective GENCODE feature name.
    #[arg(long, value_name = "STRING", default_value = "exon")]
    exon_feature_name: String,

    /// Name of the feature that represents a gene region in the GFF file.
    /// Defaults to the repective GENCODE feature name.
    #[arg(long, value_name = "STRING", default_value = "gene")]
    gene_feature_name: String,
}

//==============================================//
// Dynamic allocation of quality control facets //
//==============================================//

/// Dynamically compiles the record-based quality control facets that should be
/// run for this invocation of the command line tool.
pub fn get_record_based_qc_facets<'a>(
    features_gff: Option<PathBuf>,
    feature_names: &'a FeatureNames,
    header: &'a Header,
    reference_genome: Rc<Box<dyn ReferenceGenome>>,
) -> anyhow::Result<Vec<Box<dyn RecordBasedQualityControlFacet + 'a>>> {
    // Default facets that are loaded within the qc subcommand.
    let mut facets: Vec<Box<dyn RecordBasedQualityControlFacet>> = vec![
        Box::new(GeneralMetricsFacet::default()),
        Box::new(TemplateLengthFacet::with_capacity(1024)),
        Box::new(GCContentFacet::default()),
        Box::new(QualityScoreFacet::default()),
    ];

    // Optionally load the Genomic Features facet if the GFF file is provided.
    if let Some(s) = features_gff {
        facets.push(Box::new(GenomicFeaturesFacet::try_from(
            s,
            feature_names,
            header,
            reference_genome,
        )?));
    }

    Ok(facets)
}

/// Dynamically compiles the sequence-based quality control facets that should
/// be run for this invocation of the command line tool.
pub fn get_sequence_based_qc_facets<'a>(
    reference_fasta: Option<PathBuf>,
    reference_genome: Rc<Box<dyn ReferenceGenome>>,
) -> anyhow::Result<Vec<Box<dyn SequenceBasedQualityControlFacet + 'a>>> {
    // Default facets that are loaded within the qc subcommand.
    let mut facets: Vec<Box<dyn SequenceBasedQualityControlFacet>> =
        vec![Box::new(CoverageFacet::new(reference_genome))];

    // Optionally load the Edits facet if a reference FASTA is provided.
    if let Some(fasta) = reference_fasta {
        facets.push(Box::new(EditsFacet::try_from(fasta)?))
    }

    Ok(facets)
}

//==============================//
// Prepares the `qc` subcommand //
//==============================//

/// Prepares the arguments for running the main `qc` subcommand.
pub fn qc(args: QcArgs) -> anyhow::Result<()> {
    info!("Starting qc command...");
    debug!("Arguments:");

    //=============//
    // Source Path //
    //=============//

    let src: PathBuf = args.src;
    debug!("  [*] Source: {}", src.display());

    //==================//
    // Reference Genome //
    //==================//

    let provided_reference_genome = args.reference_genome;

    let reference_genome = match get_reference_genome(&provided_reference_genome) {
        Some(s) => Rc::new(s),
        None => bail!(
            "reference genome is not supported: {}. \
            Did you set the correct reference genome?. \
            Use the `list reference-genomes` subcommand to see supported reference genomes.",
            provided_reference_genome,
        ),
    };
    debug!("  [*] Reference genome: {}", provided_reference_genome);

    //=================//
    // Reference FASTA //
    //=================//

    let reference_fasta = args.reference_fasta;
    debug!("  [*] Reference FASTA: {:?}", reference_fasta);

    //==============//
    // Features GFF //
    //==============//

    let features_gff = args.features_gff;
    debug!("  [*] Features GFF : {:?}", features_gff);

    //===============//
    // Output Prefix //
    //===============//

    // Default is the name of the file.
    let output_prefix = args.output_prefix.unwrap_or_else(|| {
        src.file_name()
            .unwrap()
            .to_os_string()
            .into_string()
            .unwrap()
    });
    debug!("  [*] Output prefix: {}", output_prefix);

    //==========================//
    // Feature GFF Column Names //
    //==========================//

    let feature_names = FeatureNames::new(
        args.five_prime_utr_feature_name,
        args.three_prime_utr_feature_name,
        args.coding_sequence_feature_name,
        args.exon_feature_name,
        args.gene_feature_name,
    );

    //==================//
    // Output Directory //
    //==================//

    let output_directory = match args.output_directory {
        Some(p) => p,
        None => std::env::current_dir()?,
    };
    debug!("  [*] Output directory: {}", output_directory.display());

    //===================//
    // Number of Records //
    //===================//

    let num_records = match args.num_records {
        Some(n) => {
            debug!("Reading a maximum of {} records in the first pass.", n);
            NumberOfRecords::Some(n)
        }
        None => {
            debug!("Reading all available records in the first pass.");
            NumberOfRecords::All
        }
    };

    app(
        src,
        reference_fasta,
        features_gff,
        reference_genome,
        output_prefix,
        output_directory,
        num_records,
        feature_names,
    )
}

//==============//
// Main program //
//==============//

/// Runs the main program for the `qc` subcommand.
#[allow(clippy::too_many_arguments)]
fn app(
    src: PathBuf,
    reference_fasta: Option<PathBuf>,
    features_gff: Option<PathBuf>,
    reference_genome: Rc<Box<dyn ReferenceGenome>>,
    output_prefix: String,
    output_directory: PathBuf,
    num_records: NumberOfRecords,
    feature_names: FeatureNames,
) -> anyhow::Result<()> {
    //=====================================================//
    // Preprocessing: set up file handles and prepare file //
    //=====================================================//

    let mut reader = File::open(&src).map(bam::Reader::new)?;
    // This check is here simply so that, if the BAM index does not exist, we
    // don't complete the first pass before erroring out. It's not strictly
    // needed for this first pass as we aren't doing random access throughout
    // the file.
    let _ = bai::read(&src.with_extension("bam.bai")).with_context(|| "bam index")?;

    let ht = reader.read_header()?;
    let header = parse_header(ht);

    let reference_sequences = reader.read_reference_sequences()?;

    if !output_directory.exists() {
        std::fs::create_dir_all(output_directory.clone())
            .expect("Could not create output directory.");
    }

    //=====================================================//
    // Preprocessing: reference sequence concordance check //
    //=====================================================//

    let supported_sequences = get_all_sequences(Rc::clone(&reference_genome));

    for (sequence, _) in reference_sequences {
        if !supported_sequences
            .iter()
            .map(|s| s.name())
            .any(|x| x == sequence)
        {
            bail!(
                "Sequence \"{}\" not found in specified reference genome. \
                Did you set the correct reference genome?",
                sequence
            );
        }
    }

    //===========================================================//
    // First pass: print out which facets we're going to analyze //
    //===========================================================//

    let mut record_facets = get_record_based_qc_facets(
        features_gff,
        &feature_names,
        &header,
        Rc::clone(&reference_genome),
    )?;
    info!("First pass with the following facets enabled:");
    for facet in &record_facets {
        info!("  [*] {}, {:?}", facet.name(), facet.computational_load());
    }

    //====================================================================//
    // First pass: processes every record, accumulating QC stats as we go //
    //====================================================================//

    debug!("Starting first pass for QC stats.");
    let mut record_count = 0;

    for result in reader.records() {
        let record = result?;

        for facet in &mut record_facets {
            facet.process(&record)?;
        }

        record_count += 1;
        if record_count % 1_000_000 == 0 {
            info!(
                "  [*] Processed {} records.",
                record_count.to_formatted_string(&Locale::en),
            );
        }

        if let NumberOfRecords::Some(n) = num_records {
            if record_count > n {
                break;
            }
        }
    }

    info!(
        "Processed {} records in the first pass.",
        record_count.to_formatted_string(&Locale::en)
    );

    //================================//
    // First pass: summarize qc stats //
    //================================//

    info!("Summarizing quality control facets for the first pass.");
    for facet in &mut record_facets {
        facet.summarize()?;
    }

    //============================================================//
    // Second pass: print out which facets we're going to analyze //
    //============================================================//

    let mut sequence_facets =
        get_sequence_based_qc_facets(reference_fasta, Rc::clone(&reference_genome))?;
    info!("Second pass with the following facets enabled:");
    for facet in &sequence_facets {
        info!("  [*] {}, {:?}", facet.name(), facet.computational_load());
    }

    //===================================================//
    // Second pass: set up file handles and prepare file //
    //===================================================//

    let mut reader = File::open(&src).map(bam::Reader::new)?;
    let index = bai::read(&src.with_extension("bam.bai")).with_context(|| "bam index")?;

    for (name, seq) in header.reference_sequences() {
        let start = Position::MIN;
        let end = Position::try_from(usize::from(seq.length()))?;

        info!("Starting sequence {} ", name);
        let mut processed = 0;

        debug!("  [*] Setting up sequence.");
        for facet in &mut sequence_facets {
            if facet.supports_sequence_name(name) {
                facet.setup(seq)?;
            }
        }

        let query = reader.query(
            header.reference_sequences(),
            &index,
            &Region::new(name, start..=end),
        )?;

        debug!("  [*] Processing records from sequence.");
        for result in query {
            let record = result?;
            for facet in &mut sequence_facets {
                if facet.supports_sequence_name(name) {
                    facet.process(seq, &record)?;
                }
            }

            processed += 1;

            if processed % 1_000_000 == 0 {
                info!(
                    "  [*] Processed {} records for this sequence.",
                    record_count.to_formatted_string(&Locale::en),
                );
            }
        }

        debug!("  [*] Tearing down sequence.");
        for facet in &mut sequence_facets {
            if facet.supports_sequence_name(name) {
                facet.teardown(seq)?;
            }
        }
    }

    //=====================================//
    // Finalize: write all results to file //
    //=====================================//

    let mut results = Results::default();

    for facet in &record_facets {
        facet.aggregate(&mut results);
    }

    for facet in &mut sequence_facets {
        facet.aggregate(&mut results);
    }

    results.write(output_prefix, &output_directory)?;

    Ok(())
}
