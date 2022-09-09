pub mod instrument;

use clap::Command;

pub use self::instrument::derive;

pub fn get_command<'a>() -> Command<'a> {
    let derive_instrument_cmd = self::instrument::get_command();
    Command::new("derive")
        .about("Forensic analysis tool useful for backwards computing information from next-generation sequencing data.")
        .subcommand_required(true)
        .subcommand(derive_instrument_cmd)
}
