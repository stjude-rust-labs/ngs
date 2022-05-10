use futures::TryStreamExt;
use std::{collections::HashSet, io};

use crate::derive::instrument::{compute, reads};

use clap::ArgMatches;
use noodles_bam as bam;
use noodles_sam::AlignmentRecord;
use tokio::fs::File;
use tracing::info;

/// Main function
async fn app(src: &str) -> io::Result<()> {
    let mut instrument_names = HashSet::new();
    let mut flowcell_names = HashSet::new();

    let mut reader = File::open(src).await.map(bam::AsyncReader::new)?;
    reader.read_header().await?;
    reader.read_reference_sequences().await?;

    let mut samples: usize = 0;
    let mut records = reader.records();
    while let Some(record) = records.try_next().await? {
        if let Some(read_name) = record.read_name() {
            let read = reads::IlluminaReadName::from(read_name.to_string());
            instrument_names.insert(read.instrument_name);
            if let Some(fc) = read.flowcell {
                flowcell_names.insert(fc);
            }
        }

        samples += 1;
        if samples > 10000 {
            break;
        }
    }

    let result = compute::predict(instrument_names, flowcell_names);
    let output = serde_json::to_string_pretty(&result).unwrap();
    println!("{}", output);
    Ok(())
}

pub fn derive(matches: &ArgMatches) -> io::Result<()> {
    let src = matches
        .value_of("src")
        .expect("Could not parse the arguments that were passed in for src.");
    info!("Starting derive subcommand.");

    let rt = tokio::runtime::Builder::new_current_thread()
        .enable_all()
        .build()
        .unwrap();

    let app = app(src);
    rt.block_on(app)
}
