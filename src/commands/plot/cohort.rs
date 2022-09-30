use clap::{ArgMatches, Command};

pub fn get_command<'a>() -> Command<'a> {
    Command::new("cohort").about("Plots cohort-level information produced by the `ngs qc` command.")
}

pub fn plot(matches: &ArgMatches) -> anyhow::Result<()> {
    eprintln!("Hello from plot cohort!");
    Ok(())
}
