use anyhow::bail;
use clap::{builder::PossibleValuesParser, Arg, ArgMatches, Command};

use prettytable::{row, Table};

use crate::lib::utils::genome::get_all_reference_genomes;

pub fn get_command<'a>() -> Command<'a> {
    Command::new("list")
        .about("Utility to list various supported items in this command line tool.")
        .arg(
            Arg::new("subject")
                .takes_value(true)
                .help("The subject which you want to list values for.")
                .value_parser(PossibleValuesParser::new(["reference-genomes"]))
                .required(true),
        )
}

pub fn list(matches: &ArgMatches) -> anyhow::Result<()> {
    let subject: &String = matches.get_one("subject").expect("could not parse subject");

    match subject.as_str() {
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
