use tracing::error;

#[allow(dead_code)]
pub enum ExitCode {
    /// Indicates that invalid data was supplied to the given subcommand.
    InvalidInputData = 1,
}

#[allow(dead_code)]
pub fn exit<I>(message: I, code: ExitCode) -> !
where
    I: tracing::Value,
{
    error!(message);
    std::process::exit(code as i32);
}
