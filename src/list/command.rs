use anyhow::bail;
use clap::{builder::PossibleValuesParser, Args};

use prettytable::{row, Table};

use crate::utils::genome::get_all_reference_genomes;

#[derive(Args)]
pub struct ListArgs {
    /// The subject which you want to list values for.
    #[arg(value_parser = PossibleValuesParser::new(["reference-genomes"]))]
    subject: String,
}

pub fn list(args: ListArgs) -> anyhow::Result<()> {
    match args.subject.as_str() {
        "reference-genomes" => {
            let mut table = Table::new();

            table.add_row(row!["Name", "Source", "Basis"]);
            for reference in get_all_reference_genomes() {
                table.add_row(row![
                    reference.name(),
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
