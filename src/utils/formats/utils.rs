//! Utilities related to the next-generation sequencing file formats.

use noodles::sam::Header;

/// Utility struct which contains both the raw header (before being parsed, as a
/// string), and the parsed header (as a `noodles::sam::Header`).
pub struct RawAndParsedHeaders {
    /// The raw, unprocessed (and uncorrected) header [`String`] from the file.
    pub raw: String,
    /// The parsed, processed (and corrected) header [`Header`] from the file.
    pub parsed: Header,
}

#[derive(PartialEq, Eq)]
/// Utility enum to formalize in types whether or not to check for an index.
pub enum IndexCheck {
    /// Opens and parses the associated index to ensure it is well-formed.
    Full,

    /// Opens the index and just checks that the header is well-formed.
    HeaderOnly,

    /// _Does not_ check for an index file when opening and parsing the BAM file.
    None,
}
