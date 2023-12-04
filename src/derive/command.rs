//! Functionality related to the `ngs derive` subcommand itself.

pub mod encoding;
pub mod endedness;
pub mod instrument;
pub mod readlen;

use clap::Args;
use clap::Subcommand;

//===============//
// Command setup //
//===============//

/// Command line arguments for `ngs derive`.
#[derive(Args)]
pub struct DeriveArgs {
    /// The subcommand for `ngs derive`.
    #[command(subcommand)]
    pub subcommand: DeriveSubcommand,
}

/// All possible subcommands for `ngs derive`.
#[derive(Subcommand)]
pub enum DeriveSubcommand {
    /// Derives the quality score encoding used to produce the file.
    Encoding(self::encoding::DeriveEncodingArgs),

    /// Derives the endedness of the file.
    Endedness(self::endedness::DeriveEndednessArgs),

    /// Derives the instrument used to produce the file.
    Instrument(self::instrument::DeriveInstrumentArgs),

    /// Derives the read length of the file.
    Readlen(self::readlen::DeriveReadlenArgs),
}
