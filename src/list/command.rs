//! Functionality related to the `ngs list` command itself.

use anyhow::bail;
use clap::{builder::PossibleValuesParser, Args};

use prettytable::{row, Table};

use crate::utils::genome::get_all_reference_genomes;

//========================//
// Command-line arguments //
//========================//

/// Command line arguments for `ngs list`.
#[derive(Args)]
pub struct ListArgs {
    /// The subject which you want to list values for.
    #[arg(value_parser = PossibleValuesParser::new(["genomes"]))]
    subject: String,
}

//==============//
// Main command //
//==============//

/// Main method for the `ngs list` subcommand.
pub fn list(args: ListArgs) -> anyhow::Result<()> {
    match args.subject.as_str() {
        "genomes" => {
            let mut table = Table::new();

            table.add_row(row!["Name", "Triplet ID", "Source", "Basis"]);
            for reference in get_all_reference_genomes() {
                table.add_row(row![
                    reference.name(),
                    reference.triplet_id(),
                    reference.source(),
                    reference.basis(),
                ]);
            }

            table.printstd();

            Ok(())
        }
        s => bail!("Unsupported subject: {}", s),
    }
}
