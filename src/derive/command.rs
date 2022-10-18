//! Functionality related to the `ngs derive` subcommand itself.

pub mod instrument;

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
    /// Derives the instrument used to produce the file.
    Instrument(self::instrument::DeriveInstrumentArgs),
}
