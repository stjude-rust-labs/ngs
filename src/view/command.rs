//! Functionality related to the `ngs view` command itself.

use std::path::PathBuf;

use anyhow::bail;
use clap::Args;

use crate::utils::formats::BioinformaticsFileFormat;

//========================//
// Command-line arguments //
//========================//

/// Command line arguments for `ngs index`.
#[derive(Args)]
pub struct ViewArgs {
    /// Path to the file to view.
    #[arg(value_name = "BAM/CRAM")]
    src: PathBuf,

    /// If available, the query location for this view.
    query: Option<String>,

    /// If available, the FASTA reference file used to generate the file.
    #[arg(short, long)]
    reference_fasta: Option<PathBuf>,

    /// Shows the header for the file.
    #[arg(short = 'H')]
    show_header: bool,
}

//==============//
// Main command //
//==============//

/// Main method for the `ngs view` subcommand.
pub fn view(args: ViewArgs) -> anyhow::Result<()> {
    let src = args.src;
    let query = args.query;
    let reference_fasta = args.reference_fasta;
    let show_header = args.show_header;

    match BioinformaticsFileFormat::try_detect(&src) {
        Some(BioinformaticsFileFormat::BAM) => super::bam::view(src, query, show_header),
        Some(BioinformaticsFileFormat::CRAM) => {
            if let Some(reference_fasta) = reference_fasta {
                super::cram::view(src, query, reference_fasta, show_header)
            } else {
                bail!("Reference FASTA is required to view a CRAM file.")
            }
        }
        Some(format) => {
            bail!(
                "{} files are not supported by this command. This may be \
                because we haven't supported this file format yet or because \
                it does not make sense to view a file of this kind. \
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
    }
}
