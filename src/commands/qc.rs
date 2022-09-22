use futures::TryStreamExt;
use noodles_sam::Header;
use tokio::fs::File;

use std::path::PathBuf;

use clap::{Arg, ArgMatches, Command};
use noodles_bam as bam;
use num_format::{Locale, ToFormattedString};
use tracing::{debug, error, info};

use crate::lib::{
    qc::{
        features::GenomicFeaturesFacet, gc_content::GCContentFacet,
        general::metrics::GeneralMetricsFacet, results::Results,
        template_length::TemplateLengthFacet, RecordBasedQualityCheckFacet,
    },
    utils::formats::sam::parse_header,
};

/// A utility struct for passing feature name arguments from the command line
/// around more easily.
pub struct FeatureNames {
    pub five_prime_utr_feature_name: String,
    pub three_prime_utr_feature_name: String,
    pub coding_sequence_feature_name: String,
    pub exon_feature_name: String,
    pub gene_feature_name: String,
}

impl FeatureNames {
    pub fn new<I>(
        five_prime_utr_feature_name: I,
        three_prime_utf_feature_name: I,
        coding_sequence_feature_name: I,
        exon_feature_name: I,
        gene_feature_name: I,
    ) -> Self
    where
        I: Into<String>,
    {
        FeatureNames {
            five_prime_utr_feature_name: five_prime_utr_feature_name.into(),
            three_prime_utr_feature_name: three_prime_utf_feature_name.into(),
            coding_sequence_feature_name: coding_sequence_feature_name.into(),
            exon_feature_name: exon_feature_name.into(),
            gene_feature_name: gene_feature_name.into(),
        }
    }
}

/// Dynamically compiles the quality check facets that should be run for this
/// invocation of the command line tool.
///
/// # Arguments
///
/// * `features_gff` — Optionally, the path to a GFF file (if using the Genomic
///   Features quality check facet).
/// * `feature_names` — Feature names passed in from the command line. Necessary
///   for looking up the correct features from the GFF file.
/// * `header` — The SAM header; useful for looking up the sequence names from
///   the respective sequence ids.
pub fn get_facets<'a>(
    features_gff: Option<&str>,
    feature_names: &'a FeatureNames,
    header: &'a Header,
) -> anyhow::Result<Vec<Box<dyn RecordBasedQualityCheckFacet + 'a>>> {
    // Default facets that are loaded within the qc subcommand.
    let mut facets: Vec<Box<dyn RecordBasedQualityCheckFacet>> = vec![
        Box::new(GeneralMetricsFacet::default()),
        Box::new(TemplateLengthFacet::with_capacity(1024)),
        Box::new(GCContentFacet::default()),
    ];

    // Optionally load the Genomic Features facet if the GFF file is provided.
    if let Some(s) = features_gff {
        facets.push(Box::new(GenomicFeaturesFacet::try_from(
            s,
            feature_names,
            header,
        )?));
    }

    Ok(facets)
}

/// Gets the command line arguments for the `qc` subcommand.
pub fn get_command<'a>() -> Command<'a> {
    Command::new("qc")
        .about("Generates quality control metrics for BAM files.")
        .arg(
            Arg::new("src")
                .help("Source BAM file to perform QC on.")
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
                .required(true)
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
            Arg::new("max-records")
                .long("--max-records")
                .short('m')
                .help("Maximum number of records to process.")
                .takes_value(true)
                .required(false),
        )
        .arg(
            Arg::new("five-prime-utr-feature-name")
                .long("--five-prime-utr-feature-name")
                .help(concat!(
                    "Name of the feature that represents a five prime",
                    "UTR region in the GFF file. The default is to use the",
                    "GENCODE feature name."
                ))
                .takes_value(true)
                .default_value("five_prime_UTR"),
        )
        .arg(
            Arg::new("three-prime-utr-feature-name")
                .long("--three-prime-utr-feature-name")
                .help(concat!(
                    "Name of the feature that represents a three prime",
                    "UTR region in the GFF file. The default is to use the",
                    "GENCODE feature name."
                ))
                .takes_value(true)
                .default_value("three_prime_UTR"),
        )
        .arg(
            Arg::new("coding-sequence-feature-name")
                .long("--coding-sequence-feature-name")
                .help(concat!(
                    "Name of the feature that represents a coding sequence",
                    "region in the GFF file. The default is to use the",
                    "GENCODE feature name."
                ))
                .takes_value(true)
                .default_value("CDS"),
        )
        .arg(
            Arg::new("exon-feature-name")
                .long("--exon-feature-name")
                .help(concat!(
                    "Name of the feature that represents an exonic",
                    "region in the GFF file. The default is to use the",
                    "GENCODE feature name."
                ))
                .takes_value(true)
                .default_value("exon"),
        )
        .arg(
            Arg::new("gene-feature-name")
                .long("--gene-feature-name")
                .help(concat!(
                    "Name of the feature that represents an gene",
                    "region in the GFF file. The default is to use the",
                    "GENCODE feature name."
                ))
                .takes_value(true)
                .default_value("gene"),
        )
}

/// Prepares the arguments for running the main `qc` subcommand.
pub fn qc(matches: &ArgMatches) -> anyhow::Result<()> {
    info!("Starting qc command...");
    let src = matches
        .value_of("src")
        .expect("Could not parse the arguments that were passed in for src.");

    let features_gff = matches.value_of("features-gff");

    let output_prefix = matches
        .value_of("output-prefix")
        .expect("Did not receive any output prefix from args.");

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

    let output_directory = if let Some(m) = matches.value_of("output-directory") {
        PathBuf::from(m)
    } else {
        std::env::current_dir().expect("Could not retrieve the current working directory.")
    };

    let max_records = if let Some(m) = matches.value_of("max-records") {
        let res = m.parse::<i64>().unwrap();
        debug!("Reading a maximum of {} records.", res);
        res
    } else {
        debug!("Reading all available records.");
        -1
    };

    if !output_directory.exists() {
        std::fs::create_dir_all(output_directory.clone())
            .expect("Could not create output directory.");
    }

    let rt = tokio::runtime::Builder::new_multi_thread()
        .enable_all()
        .build()
        .unwrap();

    let app = app(
        src,
        features_gff,
        output_prefix,
        output_directory,
        max_records,
        feature_names,
    );

    rt.block_on(app)
}

/// Runs the `qc` subcommand.
///
/// # Arguments
///
/// * `src` — The filepath to the NGS file to run QC on.
/// * `features_gff` — Optionally, the path to a GFF gene model file. Useful
///   when you want to run the Genomic Features facet.
/// * `output_prefix` — Output prefix for all files generated by this
///   subcommand.
/// * `output_directory` — Root folder to output all of the results files to.
/// * `max_records` — Maximum number of records to process. Anything less than 0
///   is considered infinite.
/// * `feature_names` — Feature names for lookup within the GFF file.
async fn app(
    src: &str,
    features_gff: Option<&str>,
    output_prefix: &str,
    output_directory: PathBuf,
    max_records: i64,
    feature_names: FeatureNames,
) -> anyhow::Result<()> {
    let mut reader = File::open(src).await.map(bam::AsyncReader::new)?;

    let ht = reader.read_header().await?;
    let header = parse_header(ht);

    reader.read_reference_sequences().await?;

    let mut facets = get_facets(features_gff, &feature_names, &header)?;
    info!("");
    info!("Running with the following facets enabled:");
    info!("");
    for facet in &facets {
        info!(" [*] {}, {:?}", facet.name(), facet.computational_load());
    }
    info!("");

    //===============================================================//
    // Processes each of the records, accumulating QC stats as we go //
    //===============================================================//

    debug!("Accumulating QC stats.");
    let mut record_count = 0;
    let mut records = reader.lazy_records();

    while let Some(record) = records.try_next().await? {
        for facet in &mut facets {
            match facet.process(&record) {
                Ok(_) => {}
                Err(e) => {
                    error!("[{}] {}", facet.name(), e.message);
                    std::process::exit(1);
                }
            }
        }

        record_count += 1;
        if record_count % 1_000_000 == 0 && record_count > 0 {
            info!(
                "Processed {} records.",
                record_count.to_formatted_string(&Locale::en),
            );
        }

        if max_records > -1 && record_count >= max_records {
            break;
        }
    }

    info!(
        "Processed {} records.",
        record_count.to_formatted_string(&Locale::en)
    );

    info!("Summarizing quality control facets.");
    for facet in &mut facets {
        facet.summarize().unwrap();
    }

    let mut results = Results::default();

    for facet in &facets {
        facet.aggregate_results(&mut results);
    }

    results.write(String::from(output_prefix), &output_directory)?;

    Ok(())
}
