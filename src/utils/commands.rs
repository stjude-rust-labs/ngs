use clap::{arg, Command};

/// Utility method to add the verbosity arguments to any subcommands passed to Clap.
///
/// # Arguments
///
/// * `subcommand` — The Clap subcommand to add these arguments to.
pub fn add_verbosity_args(subcommand: Command) -> Command {
    subcommand
        .arg(arg!(-q --quiet "Only errors are printed to the stderr stream."))
        .arg(
            arg!(-v --verbose "All available information, including debug information, is \
                printed to stderr."),
        )
}
