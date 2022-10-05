use std::{
    fs::File,
    path::PathBuf,
    rc::Rc,
    sync::{Arc, Mutex},
};

use anyhow::{bail, Context};
use clap::{value_parser, Arg, ArgMatches, Command};
use itertools::Itertools;
use noodles_bam::{self as bam, bai};
use noodles_core::{Position, Region};
use noodles_sam::{alignment::Record, Header};
use num_format::{Locale, ToFormattedString};
use rayon::prelude::*;
use tracing::{debug, info};

use crate::lib::{
    qc::{
        record_based::features::{FeatureNames, GenomicFeaturesFacet},
        record_based::gc_content::GCContentFacet,
        record_based::general::metrics::GeneralMetricsFacet,
        record_based::quality_scores::QualityScoreFacet,
        record_based::template_length::TemplateLengthFacet,
        results::Results,
        // sequence_based::coverage::CoverageFacet,
        // sequence_based::edits::EditsFacet,
        RecordBasedQualityCheckFacet,
        SequenceBasedQualityCheckFacet,
    },
    utils::{
        formats::sam::parse_header,
        genome::{get_all_sequences, get_reference_genome, ReferenceGenome},
    },
};

//============================================//
// Dynamic allocation of quality check facets //
//============================================//

/// Dynamically compiles the record-based quality check facets that should be run for this
/// invocation of the command line tool.
pub fn get_record_based_qc_facets<'a>(
    features_gff: Option<&str>,
    feature_names: &'a FeatureNames,
    header: &'a Header,
    reference_genome: Rc<Box<dyn ReferenceGenome>>,
) -> anyhow::Result<Arc<Vec<Mutex<Box<dyn RecordBasedQualityCheckFacet + Send + Sync + 'a>>>>> {
    // Default facets that are loaded within the qc subcommand.
    let mut facets: Vec<Mutex<Box<dyn RecordBasedQualityCheckFacet + Send + Sync>>> = vec![
        Mutex::new(Box::new(GeneralMetricsFacet::default())),
        Mutex::new(Box::new(TemplateLengthFacet::with_capacity(1024))),
        Mutex::new(Box::new(GCContentFacet::default())),
        Mutex::new(Box::new(QualityScoreFacet::default())),
    ];

    // Optionally load the Genomic Features facet if the GFF file is provided.
    if let Some(s) = features_gff {
        facets.push(Mutex::new(Box::new(GenomicFeaturesFacet::try_from(
            s,
            feature_names,
            header,
            reference_genome,
        )?)));
    }

    Ok(Arc::new(facets))
}

/// Dynamically compiles the sequence-based quality check facets that should be run for this
/// invocation of the command line tool.
// pub fn get_sequence_based_qc_facets<'a>(
//     reference_fasta: Option<&PathBuf>,
//     header: &'a Header,
//     reference_genome: Rc<Box<dyn ReferenceGenome>>,
// ) -> anyhow::Result<Vec<Box<dyn SequenceBasedQualityCheckFacet<'a> + 'a>>> {
//     // Default facets that are loaded within the qc subcommand.
//     let mut facets: Vec<Box<dyn SequenceBasedQualityCheckFacet<'_>>> =
//         vec![Box::new(CoverageFacet::new(reference_genome))];

//     // Optionally load the Edits facet if a reference FASTA is provided.
//     if let Some(fasta) = reference_fasta {
//         facets.push(Box::new(EditsFacet::try_from(fasta, header)?))
//     }

//     Ok(facets)
// }

//========================//
// Command line arguments //
//========================//

/// Gets the command line arguments for the `qc` subcommand.
pub fn get_command<'a>() -> Command<'a> {
    Command::new("qc")
        .about("Generates quality control metrics for BAM files.")
        .arg(
            Arg::new("src")
                .help("Source BAM file to perform QC on.")
                .value_parser(value_parser!(PathBuf))
                .required(true),
        )
        .arg(
            Arg::new("reference-fasta")
                .long("--reference-fasta")
                .short('r')
                .help("Reference FASTA for edit lookups.")
                .value_parser(value_parser!(PathBuf))
                .takes_value(true),
        )
        .arg(
            Arg::new("reference-genome")
                .long("--reference-genome")
                .help("Reference genome used as the basis for the file.")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::new("features-gff")
                .long("--features-gff")
                .short('f')
                .help("Features GFF file.")
                .takes_value(true),
        )
        .arg(
            Arg::new("output-prefix")
                .long("--output-prefix")
                .short('p')
                .help("Output prefix for the files that will be created.")
                .takes_value(true),
        )
        .arg(
            Arg::new("output-directory")
                .long("--output-directory")
                .short('o')
                .help("The directory to output files to.")
                .required(false)
                .takes_value(true),
        )
        .arg(
            Arg::new("num-records")
                .long("--num-records")
                .short('n')
                .help("Number of records to process in the first pass.")
                .takes_value(true)
                .required(false),
        )
        .arg(
            Arg::new("five-prime-utr-feature-name")
                .long("--five-prime-utr-feature-name")
                .help(concat!(
                    "Name of the feature that represents a five prime ",
                    "UTR region in the GFF file. The default is to use the ",
                    "GENCODE feature name."
                ))
                .takes_value(true)
                .default_value("five_prime_UTR"),
        )
        .arg(
            Arg::new("three-prime-utr-feature-name")
                .long("--three-prime-utr-feature-name")
                .help(concat!(
                    "Name of the feature that represents a three prime ",
                    "UTR region in the GFF file. The default is to use the ",
                    "GENCODE feature name."
                ))
                .takes_value(true)
                .default_value("three_prime_UTR"),
        )
        .arg(
            Arg::new("coding-sequence-feature-name")
                .long("--coding-sequence-feature-name")
                .help(concat!(
                    "Name of the feature that represents a coding sequence ",
                    "region in the GFF file. The default is to use the ",
                    "GENCODE feature name."
                ))
                .takes_value(true)
                .default_value("CDS"),
        )
        .arg(
            Arg::new("exon-feature-name")
                .long("--exon-feature-name")
                .help(concat!(
                    "Name of the feature that represents an exonic ",
                    "region in the GFF file. The default is to use the ",
                    "GENCODE feature name."
                ))
                .takes_value(true)
                .default_value("exon"),
        )
        .arg(
            Arg::new("gene-feature-name")
                .long("--gene-feature-name")
                .help(concat!(
                    "Name of the feature that represents an gene ",
                    "region in the GFF file. The default is to use the ",
                    "GENCODE feature name."
                ))
                .takes_value(true)
                .default_value("gene"),
        )
}

//==============================//
// Prepares the `qc` subcommand //
//==============================//

/// Prepares the arguments for running the main `qc` subcommand.
pub fn qc(matches: &ArgMatches) -> anyhow::Result<()> {
    info!("Starting qc command...");
    debug!("Arguments:");

    //=============//
    // Source Path //
    //=============//

    let src: &PathBuf = matches
        .get_one("src")
        .expect("Could not parse the arguments that were passed in for src.");
    debug!("  [*] Source: {}", src.display());

    //==================//
    // Reference Genome //
    //==================//

    let provided_reference_genome = matches
        .get_one::<String>("reference-genome")
        .expect("Did not receive a reference genome.");

    let reference_genome = match get_reference_genome(provided_reference_genome) {
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

    let reference_fasta = matches.get_one("reference-fasta");
    debug!("  [*] Reference FASTA: {:?}", reference_fasta);

    //==============//
    // Features GFF //
    //==============//

    let features_gff = matches.value_of("features-gff");
    debug!("  [*] Features GFF : {:?}", features_gff);

    //===============//
    // Output Prefix //
    //===============//

    // Default is the name of the file.
    let output_prefix = matches
        .value_of("output-prefix")
        .unwrap_or_else(|| src.file_name().unwrap().to_str().unwrap());
    debug!("  [*] Output prefix: {}", output_prefix);

    //==========================//
    // Feature GFF Column Names //
    //==========================//

    let five_prime_utr_feature_name = matches
        .value_of("five-prime-utr-feature-name")
        .expect("Could not parse the five prime UTR feature name.");

    let three_prime_utr_feature_name = matches
        .value_of("three-prime-utr-feature-name")
        .expect("Could not parse the three prime UTR feature name.");

    let coding_sequence_feature_name = matches
        .value_of("coding-sequence-feature-name")
        .expect("Could not parse the coding sequence feature name.");

    let exon_feature_name = matches
        .value_of("exon-feature-name")
        .expect("Could not parse the exon feature name.");

    let gene_feature_name = matches
        .value_of("gene-feature-name")
        .expect("Could not parse the gene feature name.");

    let feature_names = FeatureNames::new(
        five_prime_utr_feature_name,
        three_prime_utr_feature_name,
        coding_sequence_feature_name,
        exon_feature_name,
        gene_feature_name,
    );

    //==================//
    // Output Directory //
    //==================//

    let output_directory = if let Some(m) = matches.value_of("output-directory") {
        PathBuf::from(m)
    } else {
        std::env::current_dir().expect("Could not retrieve the current working directory.")
    };
    debug!("  [*] Output directory: {}", output_directory.display());

    //===================//
    // Number of Records //
    //===================//

    let num_records = if let Some(m) = matches.value_of("num-records") {
        let res = m.parse::<i64>()?;
        debug!("Reading a maximum of {} records in the first pass.", res);
        res
    } else {
        debug!("Reading all available records in the first pass.");
        -1
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
    src: &PathBuf,
    reference_fasta: Option<&PathBuf>,
    features_gff: Option<&str>,
    reference_genome: Rc<Box<dyn ReferenceGenome>>,
    output_prefix: &str,
    output_directory: PathBuf,
    num_records: i64,
    feature_names: FeatureNames,
) -> anyhow::Result<()> {
    //=====================================================//
    // Preprocessing: set up file handles and prepare file //
    //=====================================================//

    let mut reader = File::open(src).map(bam::Reader::new)?;

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

    let record_facets = get_record_based_qc_facets(
        features_gff,
        &feature_names,
        &header,
        Rc::clone(&reference_genome),
    )?;
    info!("First pass with the following facets enabled:");
    for result in record_facets.iter() {
        let facet = result.lock().unwrap();
        info!("  [*] {}, {:?}", facet.name(), facet.computational_load());
    }

    //====================================================================//
    // First pass: processes every record, accumulating QC stats as we go //
    //====================================================================//

    debug!("Starting first pass for QC stats.");
    let mut record_count = 0;

    for chunk in &reader.records().into_iter().chunks(1_000_000) {
        let records: Vec<Record> = chunk.map(|x| x.unwrap()).collect();
        let results: Vec<_> = records
            .par_iter()
            .map(|record| -> anyhow::Result<()> {
                for facet in record_facets.iter() {
                    facet.lock().unwrap().process(record)?;
                }

                Ok(())
            })
            .collect();

        record_count += 1_000_000;
        info!(
            "  [*] Processed {} records.",
            record_count.to_formatted_string(&Locale::en),
        );

        if num_records > -1 && record_count >= num_records {
            break;
        }
    }

    // for result in reader.records() {
    //     let record = result?;

    //     let mut facets = record_facets.lock().unwrap();
    //     for facet in &mut facets.iter_mut() {
    //         facet.process(&record)?;
    //     }

    //     record_count += 1;
    //     if record_count % 1_000_000 == 0 && record_count > 0 {
    //         info!(
    //             "  [*] Processed {} records.",
    //             record_count.to_formatted_string(&Locale::en),
    //         );
    //     }

    //     if num_records > -1 && record_count >= num_records {
    //         break;
    //     }
    // }

    info!(
        "Processed {} records in the first pass.",
        record_count.to_formatted_string(&Locale::en)
    );

    //================================//
    // First pass: summarize qc stats //
    //================================//

    info!("Summarizing quality control facets for the first pass.");
    // for facet in record_facets.iter() {
    //     facet.summarize()?;
    // }

    //============================================================//
    // Second pass: print out which facets we're going to analyze //
    //============================================================//

    // let mut sequence_facets =
    //     get_sequence_based_qc_facets(reference_fasta, &header, Rc::clone(&reference_genome))?;
    // info!("Second pass with the following facets enabled:");
    // for facet in &sequence_facets {
    //     info!("  [*] {}, {:?}", facet.name(), facet.computational_load());
    // }

    //===================================================//
    // Second pass: set up file handles and prepare file //
    //===================================================//

    // let mut reader = File::open(src).map(bam::Reader::new)?;
    // let index = bai::read(src.with_extension("bam.bai")).with_context(|| "bam index")?;

    // for (name, seq) in header.reference_sequences() {
    //     let start = Position::MIN;
    //     let end = Position::try_from(usize::from(seq.len()))?;

    //     info!("Starting sequence {} ", name);
    //     let mut processed = 0;

    //     debug!("  [*] Setting up sequence.");
    //     for facet in &mut sequence_facets {
    //         if facet.supports_sequence_name(name) {
    //             facet.setup(seq)?;
    //         }
    //     }

    //     let query = reader.query(
    //         header.reference_sequences(),
    //         &index,
    //         &Region::new(name, start..=end),
    //     )?;

    //     debug!("  [*] Processing records from sequence.");
    //     for result in query {
    //         let record = result?;
    //         for facet in &mut sequence_facets {
    //             if facet.supports_sequence_name(name) {
    //                 facet.process(seq, &record)?;
    //             }
    //         }

    //         processed += 1;

    //         if processed % 1_000 == 0 {
    //             info!("  [*] Processed {} records for this sequence.", processed);
    //         }
    //     }

    //     debug!("  [*] Tearing down sequence.");
    //     for facet in &mut sequence_facets {
    //         if facet.supports_sequence_name(name) {
    //             facet.teardown(seq)?;
    //         }
    //     }
    // }

    //=====================================//
    // Finalize: write all results to file //
    //=====================================//

    let mut results = Results::default();

    for result in record_facets.iter() {
        let facet = result.lock().unwrap();
        facet.aggregate(&mut results);
    }

    // for facet in &mut sequence_facets {
    //     facet.aggregate(&mut results);
    // }

    results.write(String::from(output_prefix), &output_directory)?;

    Ok(())
}
