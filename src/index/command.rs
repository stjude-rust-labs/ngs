//! Functionality related to the `ngs index` command itself.

use std::{ffi::OsStr, path::PathBuf};

use anyhow::bail;
use clap::Args;

//========================//
// Command-line arguments //
//========================//

/// Command line arguments for `ngs index`.
#[derive(Args)]
pub struct IndexArgs {
    /// Path to the file to index.
    #[arg(value_name = "BAM")]
    src: PathBuf,
}

//==============//
// Main command //
//==============//

/// Main method for the `ngs index` subcommand.
pub fn index(args: IndexArgs) -> anyhow::Result<()> {
    let src = args.src;

    match &src.extension().and_then(OsStr::to_str) {
        Some("bam") => super::bam::index(src),
        Some(ext) => {
            bail!(
                "{} files are not supported by this command. This may be \
                because we haven't supported this file format yet or because \
                the file format is not a next-generation sequencing file. \
                If you believe this format should be supported, please search \
                for and upvote the related issue on Github (or file a new one).",
                ext.to_ascii_uppercase()
            )
        }
        None => {
            bail!("Filepaths without an extension are not supported.")
        }
    }?;
    Ok(())
}
