//! Functionality related to the `ngs convert` command itself.

use std::path::PathBuf;

use anyhow::bail;
use anyhow::Context;
use clap::arg;
use clap::Args;
use tracing::debug;

use crate::convert::bam;
use crate::convert::cram;
use crate::convert::gff;
use crate::convert::sam;
use crate::utils::args::CompressionStrategy;
use crate::utils::args::NumberOfRecords;
use crate::utils::formats::BioinformaticsFileError;
use crate::utils::formats::BioinformaticsFileFormat;

//========================//
// Command-line arguments //
//========================//

/// Command line arguments for `ngs convert`.
#[derive(Args)]
pub struct ConvertArgs {
    /// Path to the source file from which we are converting.
    from: PathBuf,

    /// Path to the destination file to which we are converting.
    to: PathBuf,

    /// Number of records to process before exiting the conversion.
    #[arg(short = 'n', long, value_name = "USIZE")]
    num_records: Option<usize>,

    /// If available, the FASTA reference file used to generate the file.
    #[arg(short, long)]
    reference_fasta: Option<PathBuf>,

    #[arg(short, long)]
    compression_strategy: Option<CompressionStrategy>,
}

/// Utility struct to join two bioinformatics file formats together as a tuple. Most
/// commonly, this is used to facilitate conversions **from** one bioinformatics file
/// format **to** another bioinformatics file format.
pub struct BioinformaticsFilePair(BioinformaticsFileFormat, BioinformaticsFileFormat);

impl BioinformaticsFilePair {
    /// The from _from_ which we are converting.
    pub fn from(&self) -> &BioinformaticsFileFormat {
        &self.0
    }

    /// The from _to_ which we are converting.
    pub fn to(&self) -> &BioinformaticsFileFormat {
        &self.1
    }
}

/// Runs the main program for the `convert` subcommand.
pub fn convert(args: ConvertArgs) -> anyhow::Result<()> {
    //===============//
    // From Filepath //
    //===============//

    let from = BioinformaticsFileFormat::try_detect(&args.from)
        .ok_or(BioinformaticsFileError::FailedParsing)
        .with_context(|| format!("from input file: {}", &args.from.display()))?;

    //=============//
    // To Filepath //
    //=============//

    let to = BioinformaticsFileFormat::try_detect(&args.to)
        .ok_or(BioinformaticsFileError::FailedParsing)
        .with_context(|| format!("to input file: {}", args.to.display()))?;

    //======================//
    // Compression Strategy //
    //======================//

    let compression_strategy = match args.compression_strategy {
        Some(cs) => cs,
        None => CompressionStrategy::Balanced,
    };

    debug!(
        "If applicable, using the {} compression strategy.",
        compression_strategy
    );

    //===================//
    // Number of Records //
    //===================//

    let num_records = NumberOfRecords::from(args.num_records);

    //==========================//
    // Bioinformatics File Pair //
    //==========================//

    let pair = BioinformaticsFilePair(from, to);

    let rt = tokio::runtime::Builder::new_current_thread()
        .enable_all()
        .build()
        .unwrap();

    match pair {
        BioinformaticsFilePair(BioinformaticsFileFormat::SAM, BioinformaticsFileFormat::BAM) => {
            rt.block_on(sam::to_bam_async(args.from, args.to, num_records))
        }
        BioinformaticsFilePair(BioinformaticsFileFormat::BAM, BioinformaticsFileFormat::SAM) => {
            rt.block_on(bam::to_sam_async(args.from, args.to, num_records))
        }
        BioinformaticsFilePair(BioinformaticsFileFormat::BAM, BioinformaticsFileFormat::CRAM) => {
            todo!()
        }
        BioinformaticsFilePair(BioinformaticsFileFormat::CRAM, BioinformaticsFileFormat::BAM) => {
            todo!()
        }
        BioinformaticsFilePair(BioinformaticsFileFormat::SAM, BioinformaticsFileFormat::CRAM) => {
            todo!()
        }
        BioinformaticsFilePair(BioinformaticsFileFormat::CRAM, BioinformaticsFileFormat::SAM) => {
            let fasta = match args.reference_fasta {
                Some(s) => s,
                None => bail!(
                    "--reference-fasta is a required argument when converting to/from a CRAM file"
                ),
            };

            rt.block_on(cram::to_sam_async(args.from, args.to, fasta, num_records))
        }
        BioinformaticsFilePair(
            BioinformaticsFileFormat::GFF,
            BioinformaticsFileFormat::GFF_BGZ,
        ) => gff::to_block_gzipped_gff(args.from, args.to, num_records, compression_strategy),
        _ => bail!(
            "Conversion from {} to {} is not currently supported",
            pair.from(),
            pair.to()
        ),
    }
}
