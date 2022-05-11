use futures::TryStreamExt;
use std::{collections::HashSet, io};

use crate::derive::instrument::{compute, reads};

use clap::ArgMatches;
use noodles_bam as bam;
use noodles_sam::{reader::record::Fields, AlignmentRecord};
use tokio::fs::File;
use tracing::info;

/// Runs the bulk of the derive instrument subcommand in an `async` context.
async fn app(src: &str, first_n_reads: Option<usize>) -> io::Result<()> {
    let mut instrument_names = HashSet::new();
    let mut flowcell_names = HashSet::new();

    let mut reader = File::open(src).await.map(bam::AsyncReader::new)?;
    reader.read_header().await?;
    reader.read_reference_sequences().await?;

    // (1) Collect instrument names and flowcell names from reads within the
    // file. Support for sampling only a portion of the reads is provided.
    let mut samples: usize = 0;
    let mut sample_max = 0;

    if let Some(s) = first_n_reads {
        sample_max = s;
    }

    let fields = Fields::READ_NAME;
    let mut records = reader.records_with_fields(fields);
    while let Some(record) = records.try_next().await? {
        if let Some(read_name) = record.read_name() {
            let read = reads::IlluminaReadName::from(read_name.to_string());
            instrument_names.insert(read.instrument_name);
            if let Some(fc) = read.flowcell {
                flowcell_names.insert(fc);
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
    let src = matches
        .value_of("src")
        .expect("Could not parse the arguments that were passed in for src.");

    let first_n_reads = matches.value_of("first_n_reads").map(|s| {
        s.parse::<usize>()
            .expect("Subsample must be specified as a parsable integer")
    });

    info!("Starting derive subcommand.");

    let rt = tokio::runtime::Builder::new_current_thread()
        .enable_all()
        .build()
        .unwrap();

    let app = app(src, first_n_reads);
    rt.block_on(app)
}
