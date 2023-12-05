//! Functionality relating to the `ngs derive encoding` subcommand itself.

use std::path::PathBuf;

use anyhow::Ok;
use clap::Args;

/// Clap arguments for the `ngs derive encoding` subcommand.
#[derive(Args)]
pub struct DeriveEncodingArgs {
    // Source NGS file (BAM or FASTQ).
    #[arg(value_name = "NGS_FILE")]
    src: PathBuf,

    /// Only examine the first n records in the file.
    #[arg(short, long, value_name = "USIZE")]
    num_records: Option<usize>,
}

/// Main function for the `ngs derive encoding` subcommand.
pub fn derive(args: DeriveEncodingArgs) -> anyhow::Result<()> {
    Ok(())
}
