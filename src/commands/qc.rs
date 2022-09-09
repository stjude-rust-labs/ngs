use futures::TryStreamExt;

use std::{io, path::PathBuf};

use clap::{Arg, ArgMatches, Command};
use noodles_bam as bam;
use num_format::{Locale, ToFormattedString};
use tokio::{fs::File, io::AsyncWriteExt};
use tracing::{debug, info};

use crate::qc::{
    gc_content::GCContentHistogram, metrics::QualityCheckMetrics,
    template_length::TemplateLengthHistogram,
};

pub fn get_command<'a>() -> Command<'a> {
    Command::new("qc")
        .about("Generates quality control metrics for BAM files.")
        .arg(Arg::new("src").help("Source file.").index(1).required(true))
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
}

async fn app(
    src: &str,
    output_prefix: &str,
    output_directory: PathBuf,
    max_records: i64,
) -> io::Result<()> {
    let mut reader = File::open(src).await.map(bam::AsyncReader::new)?;

    let mut metrics = QualityCheckMetrics::new();
    let mut template_length_histogram = TemplateLengthHistogram::with_capacity(1024);
    let mut gc_content_histogram = GCContentHistogram::new();

    reader.read_header().await?;
    reader.read_reference_sequences().await?;

    //===============================================================//
    // Processes each of the records, accumulating QC stats as we go //
    //===============================================================//

    let mut record_count = 0;
    let mut records = reader.lazy_records();
    while let Some(record) = records.try_next().await? {
        metrics.process(&record);
        template_length_histogram.process(&record);
        gc_content_histogram.process(&record);

        record_count += 1;
        if record_count % 1_000_000 == 0 && record_count > 0 {
            debug!(
                "Processed {} records.",
                record_count.to_formatted_string(&Locale::en)
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

    metrics.summarize(&template_length_histogram, &gc_content_histogram);

    //==================================//
    // Output all of the relevant files //
    //==================================//

    // (1) Write main metrics file.
    let metrics_filename = String::from(output_prefix) + ".metrics.json";
    let mut metrics_filepath = output_directory.clone();
    metrics_filepath.push(metrics_filename);

    let mut file = File::create(metrics_filepath).await?;
    let output = serde_json::to_string_pretty(&metrics).unwrap();
    file.write_all(output.as_bytes()).await?;

    // (2) Write the template length file.
    let template_length_filename = String::from(output_prefix) + ".template_lengths.json";
    let mut template_length_filepath = output_directory.clone();
    template_length_filepath.push(template_length_filename);

    let mut file = File::create(template_length_filepath).await?;
    let output = serde_json::to_string_pretty(&template_length_histogram).unwrap();
    file.write_all(output.as_bytes()).await?;

    // (3) Write the GC content file.
    let gc_content_filename = String::from(output_prefix) + ".gc_content.json";
    let mut gc_content_filepath = output_directory.clone();
    gc_content_filepath.push(gc_content_filename);

    let mut file = File::create(gc_content_filepath).await?;
    let output = serde_json::to_string_pretty(&gc_content_histogram).unwrap();
    file.write_all(output.as_bytes()).await?;

    Ok(())
}

pub fn qc(matches: &ArgMatches) -> io::Result<()> {
    info!("Starting qc command...");
    let src = matches
        .value_of("src")
        .expect("Could not parse the arguments that were passed in for src.");

    let output_prefix = matches
        .value_of("output-prefix")
        .expect("Did not receive any output prefix from args.");

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

    let app = app(src, output_prefix, output_directory, max_records);
    rt.block_on(app)
}
