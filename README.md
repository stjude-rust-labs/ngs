<p align="center">
  <h1 align="center">
    ngs
  </h1>

  <p align="center">
    <a href="https://github.com/stjude-rust-labs/ngs/actions/workflows/CI.yml" target="_blank">
      <img alt="CI: Status" src="https://github.com/stjude-rust-labs/ngs/actions/workflows/CI.yml/badge.svg" />
    </a>
    <a href="https://crates.io/crates/ngs" target="_blank">
      <img alt="crates.io version" src="https://img.shields.io/crates/v/ngs">
    </a>
    <img alt="crates.io downloads" src="https://img.shields.io/crates/d/ngs">
    <a href="https://github.com/stjude-rust-labs/ngs/blob/master/LICENSE.md" target="_blank">
      <img alt="License: MIT" src="https://img.shields.io/badge/license-MIT-blue.svg" />
    </a>
  </p>


  <p align="center">
    Command line utility for working with next-generation sequencing files. 
    <br />
    <a href="https://github.com/stjude-rust-labs/ngs/wiki"><strong>Explore the docs »</strong></a>
    <br />
    <br />
    <a href="https://github.com/stjude-rust-labs/ngs/issues/new?assignees=&labels=&template=feature_request.md&title=Descriptive%20Title&labels=enhancement">Request Feature</a>
    ·
    <a href="https://github.com/stjude-rust-labs/ngs/issues/new?assignees=&labels=&template=bug_report.md&title=Descriptive%20Title&labels=bug">Report Bug</a>
    ·
    ⭐ Consider starring the repo! ⭐
    <br />
  </p>

  <p>
    <img src="https://raw.githubusercontent.com/stjude-rust-labs/ngs/main/.github/assets/experimental-warning.png">
  </p>
</p>


## 🎨 Features

* **[`ngs derive`](https://github.com/stjude-rust-labs/ngs/wiki/ngs-derive).** Forensic analysis tool useful for backwards computing information from next-generation sequencing data.
* **[`ngs generate`](https://github.com/stjude-rust-labs/ngs/wiki/ngs-generate).** Tool to generate next-generation sequencing files.
* **[`ngs qc`](https://github.com/stjude-rust-labs/ngs/wiki/ngs-qc).** Provides tools for checking the quality of next-generation sequencing files.

## Guiding Principles

* **Modern, reliable foundation for everyday bioinformatics analysis—written in Rust.** `ngs` aims to package together a fairly comprehensive set of analysis tools and utilities for everyday work in bioinformatics. It is built with modern, multi-core systems in mind and written in Rust. Though we are not there today, we plan to work towards this goal in the future.
* **Runs on readily available hardware/software.** We aim for every subcommand within `ngs` to run within most computing environments without the need for special hardware or software. Practically, this means we've designed `ngs` to run in any UNIX-like environment that has at least four (4) cores and sixteen (16) GB of RAM. Often, tools will run with fewer resources. This design decision is important and sometimes means that `ngs` runs slower than it otherwise could.

## 📚 Getting Started

### Installation

To install the latest released version, you can simply use `cargo`.

```bash
cargo install ngs
```

To install the latest version on `main`, you can use the following command.

```bash
cargo install --locked --git https://github.com/stjude-rust-labs/ngs.git
```

## 🖥️ Development

To bootstrap a development environment, please use the following commands.

```bash
# Clone the repository
git clone git@github.com:stjude-rust-labs/ngs.git
cd ngs

# Run the command line tool using cargo.
cargo run -- -h
```

## 🚧️ Tests

```bash
# Run the project's tests.
cargo test

# Ensure the project doesn't have any linting warnings.
cargo clippy

# Ensure the project passes `cargo fmt`.
cargo fmt --check
```

## Minimum Supported Rust Version (MSRV)

The minimum supported Rust version for this project is 1.64.0.

## 🤝 Contributing

Contributions, issues and feature requests are welcome! Feel free to check
[issues page](https://github.com/stjude-rust-labs/ngs/issues).

## 📝 License

Copyright © 2021-Present [St. Jude Children's Research
Hospital](https://github.com/stjude).

This project is licensed as either [Apache 2.0][license-apache] or
[MIT][license-mit] at your discretion.

[contributing-md]: https://github.com/stjude-rust-labs/ngs/blob/master/CONTRIBUTING.md
[license-apache]: https://github.com/stjude-rust-labs/ngs/blob/master/LICENSE-APACHE
[license-mit]: https://github.com/stjude-rust-labs/ngs/blob/master/LICENSE-MIT
