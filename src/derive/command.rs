pub mod instrument;

use clap::{Args, Subcommand};

//===============//
// Command setup //
//===============//

#[derive(Args)]
pub struct DeriveArgs {
    #[command(subcommand)]
    pub subcommand: DeriveSubcommand,
}

#[derive(Subcommand)]
pub enum DeriveSubcommand {
    /// Derives the instrument used to produce the file.
    Instrument(self::instrument::Instrument),
}
