use tracing::error;

pub enum ExitCode {
    /// Indicates that invalid data was supplied to the given subcommand.
    InvalidInputData = 1,
}

pub fn exit<I>(message: I, code: ExitCode) -> !
where
    I: tracing::Value,
{
    error!(message);
    std::process::exit(code as i32);
}
