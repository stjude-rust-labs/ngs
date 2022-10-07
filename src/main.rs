#![warn(missing_docs)]
#![warn(rust_2018_idioms)]
#![warn(rust_2021_compatibility)]
//! `ngs` is a command line tool written to facilitate the analysis of
//! next-generation sequencing analysis. To learn more, visit the
//! [wiki](https://github.com/stjude-rust-labs/ngs/wiki).

use anyhow::bail;
use clap::Command;

use git_testament::{git_testament, render_testament};
use ngs::{derive, generate, list, plot, qc, utils::commands::add_verbosity_args};

git_testament!(TESTAMENT);

fn main() -> anyhow::Result<()> {
    let derive_cmd = derive::command::get_command();
    let generate_cmd = generate::command::get_command();
    let list_cmd = list::command::get_command();
    let plot_cmd = plot::command::get_command();
    let qc_cmd = qc::command::get_command();

    let matches = Command::new("ngs")
        .version(render_testament!(TESTAMENT))
        .propagate_version(true)
        .subcommand_required(true)
        .subcommand(add_verbosity_args(derive_cmd))
        .subcommand(add_verbosity_args(generate_cmd))
        .subcommand(add_verbosity_args(list_cmd))
        .subcommand(add_verbosity_args(plot_cmd))
        .subcommand(add_verbosity_args(qc_cmd))
        .get_matches();

    if let Some((name, subcommand)) = matches.subcommand() {
        let mut level = tracing::Level::INFO;
        if subcommand.get_flag("quiet") {
            level = tracing::Level::ERROR;
        } else if subcommand.get_flag("verbose") {
            level = tracing::Level::DEBUG;
        }

        let subscriber = tracing_subscriber::fmt::Subscriber::builder()
            .with_max_level(level)
            .with_writer(std::io::stderr)
            .finish();
        let _ = tracing::subscriber::set_global_default(subscriber);

        match name {
            "derive" => {
                if let Some(m) = subcommand.subcommand_matches("instrument") {
                    return derive::command::instrument::derive(m);
                } else {
                    unreachable!();
                }
            }
            "generate" => return generate::command::generate(subcommand),
            "list" => return list::command::list(subcommand),
            "plot" => {
                if let Some(m) = subcommand.subcommand_matches("sample") {
                    return plot::sample::plot(m);
                } else if let Some(m) = subcommand.subcommand_matches("cohort") {
                    return plot::cohort::plot(m);
                } else {
                    unreachable!();
                }
            }
            "qc" => return qc::command::qc(subcommand),
            s => {
                bail!("Unknown subcommand: {}", s);
            }
        }
    }

    unreachable!();
}
