use anyhow::bail;
use futures::TryStreamExt;
use std::{collections::HashSet, path::PathBuf, thread};

use clap::Args;
use noodles_bam as bam;
use tokio::fs::File;
use tracing::info;

use crate::derive::instrument::{compute, reads::IlluminaReadName};

#[derive(Args)]
pub struct Instrument {
    // Source BAM.
    #[arg(value_name = "BAM")]
    src: PathBuf,

    /// Only examine the first n records in the file.
    #[arg(short, long, value_name = "USIZE")]
    num_records: Option<usize>,

    /// Use a specific number of threads.
    #[arg(short, long, value_name = "USIZE")]
    threads: Option<usize>,
}

pub fn derive(args: Instrument) -> anyhow::Result<()> {
    let first_n_reads: Option<usize> = args.num_records;
    let threads = match args.threads {
        Some(t) => t,
        None => thread::available_parallelism().map(usize::from)?,
    };

    info!(
        "Starting derive instrument subcommand with {} threads.",
        threads
    );

    let rt = tokio::runtime::Builder::new_multi_thread()
        .worker_threads(threads)
        .build()?;

    rt.block_on(app(args.src, first_n_reads))
}

async fn app(src: PathBuf, first_n_reads: Option<usize>) -> anyhow::Result<()> {
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
                    bail!(
                        "Could not parse Illumina-formatted query names for read: {}",
                        name
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
    print!("{}", output);

    Ok(())
}
