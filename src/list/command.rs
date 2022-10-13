//! Functionality related to the `ngs list` command itself.

use anyhow::bail;
use clap::{builder::PossibleValuesParser, Args};

use prettytable::{row, Table};

use crate::{
    plot::command::{get_all_cohort_plots, get_all_sample_plots},
    utils::genome::get_all_reference_genomes,
};

//========================//
// Command-line arguments //
//========================//

/// Command line arguments for `ngs list`.
#[derive(Args)]
pub struct ListArgs {
    /// The subject which you want to list values for.
    #[arg(value_parser = PossibleValuesParser::new(["genomes", "plots"]))]
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
        "plots" => {
            let mut sample_table = Table::new();
            sample_table.add_row(row!["Name", "Type", "Description"]);

            for sample_plot in get_all_sample_plots(None)? {
                sample_table.add_row(row![
                    sample_plot.name(),
                    "Sample",
                    sample_plot.description()
                ]);
            }

            println!("Sample Plots:");
            println!();
            sample_table.printstd();
            println!();

            let mut cohort_table = Table::new();
            cohort_table.add_row(row!["Name", "Type", "Description"]);
            for cohort_plot in get_all_cohort_plots(None)? {
                cohort_table.add_row(row![
                    cohort_plot.name(),
                    "Cohort",
                    cohort_plot.description()
                ]);
            }

            println!("Cohort Plots:");
            println!();
            cohort_table.printstd();

            Ok(())
        }
        s => bail!("Unsupported subject: {}", s),
    }
}
