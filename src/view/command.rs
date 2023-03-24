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
    #[arg(value_name = "FILE")]
    src: PathBuf,

    /// If available, the query region for this view.
    query: Option<String>,

    /// If available, the FASTA reference file used to generate the file.
    #[arg(short, long)]
    reference_fasta: Option<PathBuf>,

    /// Shows either the header ("header-only"), the records ("records-only"),
    /// or both ("full").
    #[arg(short, long, value_enum, default_value = "full")]
    mode: Mode,
}

/// Mode of the view command
#[derive(clap::ValueEnum, PartialEq, Eq, Clone)]
pub enum Mode {
    /// Shows the full contents of the file (header and records).
    Full,

    /// Shows only the header of the file.
    HeaderOnly,

    /// Shows only the records within the file.
    RecordsOnly,
}

//==============//
// Main command //
//==============//

/// Main method for the `ngs view` subcommand.
pub fn view(args: ViewArgs) -> anyhow::Result<()> {
    let src = args.src;
    let query = args.query;
    let reference_fasta = args.reference_fasta;
    let mode = args.mode;

    let rt = tokio::runtime::Builder::new_current_thread()
        .enable_all()
        .build()
        .unwrap();

    match BioinformaticsFileFormat::try_detect(&src) {
        Some(BioinformaticsFileFormat::SAM) => rt.block_on(super::sam::view(src, query, mode)),
        Some(BioinformaticsFileFormat::BAM) => rt.block_on(super::bam::view(src, query, mode)),
        Some(BioinformaticsFileFormat::CRAM) => {
            if let Some(reference_fasta) = reference_fasta {
                rt.block_on(super::cram::view(src, query, reference_fasta, mode))
            } else {
                bail!(
                    "--reference-fasta is a required argument when converting to/from a CRAM file"
                )
            }
        }
        Some(BioinformaticsFileFormat::GFF) | Some(BioinformaticsFileFormat::GFF_GZ) => {
            super::gff::view(src)
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
