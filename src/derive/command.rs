pub mod instrument;

use clap::Command;

pub fn get_command() -> Command {
    let derive_instrument_cmd = self::instrument::get_command();
    Command::new("derive")
        .about("Forensic analysis tool for next-generation sequencing data")
        .subcommand_required(true)
        .subcommand(derive_instrument_cmd)
}
