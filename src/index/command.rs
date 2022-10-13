//! Functionality related to the `ngs index` command itself.

use std::path::PathBuf;

use anyhow::bail;
use clap::Args;

use crate::utils::formats::BioinformaticsFileFormat;

//========================//
// Command-line arguments //
//========================//

/// Command line arguments for `ngs index`.
#[derive(Args)]
pub struct IndexArgs {
    /// Path to the file to index.
    #[arg(value_name = "BAM/CRAM/FASTA")]
    src: PathBuf,
}

//==============//
// Main command //
//==============//

/// Main method for the `ngs index` subcommand.
pub fn index(args: IndexArgs) -> anyhow::Result<()> {
    let src = args.src;

    match BioinformaticsFileFormat::try_detect(&src) {
        Some(BioinformaticsFileFormat::BAM) => super::bam::index(src),
        Some(BioinformaticsFileFormat::CRAM) => super::cram::index(src),
        Some(BioinformaticsFileFormat::FASTA) => super::fasta::index(src),
        Some(format) => {
            bail!(
                "{} files are not supported by this command. This may be \
                because we haven't supported this file format yet or because \
                it does not make sense to index a file of this kind. \
                If you believe this format should be supported, please search \
                for and upvote the related issue on Github (or file a new one).",
                format
            )
        }
        None => {
            bail!(
                "Not able to determine bioinformatics file type for path: {}",
                src.display()
            )
        }
    }?;
    Ok(())
}
