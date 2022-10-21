//! `ngs` is a command line tool written to facilitate the analysis of
//! next-generation sequencing analysis. This package is composed of both
//! a library crate, as well as a binary crate.
//!
//! This documentation generally refers to the library crate documentation for
//! use by developers of `ngs`. If you're interested in details about using
//! the `ngs` command line tool (much more common), please visit the
//! [wiki page](https://github.com/stjude-rust-labs/ngs/wiki).
#![warn(missing_docs)]
#![warn(rust_2018_idioms)]
#![warn(rust_2021_compatibility)]

pub mod convert;
pub mod derive;
pub mod generate;
pub mod index;
pub mod list;
pub mod plot;
pub mod qc;
pub mod utils;
pub mod view;
