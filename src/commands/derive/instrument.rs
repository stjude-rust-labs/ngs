use futures::TryStreamExt;
use std::{collections::HashSet, io};

use crate::{
    derive::instrument::{compute, reads::IlluminaReadName},
    errors::{exit, ExitCode},
};

use clap::{Arg, ArgMatches, Command};
use noodles_bam as bam;
use tokio::fs::File;
use tracing::{error, info};

pub fn get_command<'a>() -> Command<'a> {
    Command::new("instrument")
        .about("Derives the instrument used to produce the file (only Illumina is supported)")
        .arg(Arg::new("src").help("Source file.").index(1).required(true))
        .arg(
            Arg::new("first-n-reads")
                .short('n')
                .long("first-n-reads")
                .takes_value(true)
                .help(
                    "Only consider the first n reads in the file. If less \
                      than or equal to zero, the whole file will be read.",
                ),
        )
}

/// Runs the bulk of the derive instrument subcommand in an `async` context.
async fn app(src: &str, first_n_reads: Option<usize>) -> io::Result<()> {
    let mut instrument_names = HashSet::new();
    let mut flowcell_names = HashSet::new();

    let mut reader = File::open(src).await.map(bam::AsyncReader::new)?;
    reader.read_header().await?;
    reader.read_reference_sequences().await?;

    // (1) Collect instrument names and flowcell names from reads within the
    // file. Support for sampling only a portion of the reads is provided.
    let mut samples = 0;
    let mut sample_max = 0;

    if let Some(s) = first_n_reads {
        sample_max = s;
    }

    let mut records = reader.lazy_records();
    while let Some(record) = records.try_next().await? {
        if let Ok(Some(read_name)) = record.read_name() {
            let name: &str = read_name.as_ref();

            match name.parse::<IlluminaReadName>() {
                Ok(read) => {
                    instrument_names.insert(read.instrument_name);
                    if let Some(fc) = read.flowcell {
                        flowcell_names.insert(fc);
                    }
                }
                Err(_) => {
                    error!(
                        "Could not parse Illumina-formatted query names for read: {}.",
                        name
                    );

                    exit(
                        "Illumina-formatted reads are expected to be colon (:) delimited with \
either five or seven fields. Please see the Wikipedia page on \
FASTQ files (https://en.wikipedia.org/wiki/FASTQ_format#Illumina_sequence_identifiers) \
for more details.",
                        ExitCode::InvalidInputData,
                    );
                }
            }
        }

        if sample_max > 0 {
            samples += 1;
            if samples > sample_max {
                break;
            }
        }
    }

    // (2) Derive the predict instrument results based on these detected
    // instrument names and flowcell names.
    let result = compute::predict(instrument_names, flowcell_names);

    // (3) Print the output to stdout as JSON (more support for different output
    // types may be added in the future, but for now, only JSON).
    let output = serde_json::to_string_pretty(&result).unwrap();
    println!("{}", output);

    Ok(())
}

/// Main function. This sets up an single-threaded, asynchronous runtime via
/// tokio and then runs the `app()` method (which contains the bulk of this
/// subcommand).
pub fn derive(matches: &ArgMatches) -> io::Result<()> {
    let src = matches.value_of("src").unwrap_or_else(|| {
        exit(
            "Could not parse the arguments that were passed in for src.",
            ExitCode::InvalidInputData,
        )
    });

    let first_n_reads = matches.value_of("first-n-reads").map(|s| {
        let num = s.parse::<usize>().unwrap_or_else(|_| {
            exit(
                "--first-n-reads must be specified as a parsable, positive integer.",
                ExitCode::InvalidInputData,
            )
        });

        if num == 0 {
            exit(
                "--first-n-reads must be greater than 0!",
                ExitCode::InvalidInputData,
            );
        }

        num
    });

    info!("Starting derive instrument subcommand.");

    let rt = tokio::runtime::Builder::new_current_thread()
        .enable_all()
        .build()
        .unwrap();

    let app = app(src, first_n_reads);
    rt.block_on(app)
}
