//! Utilities related to opening and manipulating Binary Alignment Map (BAM) files.

use std::{
    fs::File,
    io::BufReader,
    path::{Path, PathBuf},
};

use anyhow::bail;
use anyhow::Context;
use indexmap::IndexMap;
use noodles::bam;
use noodles::bam::bai;
use noodles::bgzf;
use noodles::sam::header::record::value::map::Map;
use noodles::sam::header::record::value::map::ReferenceSequence;
use noodles::sam::Header;
use tracing::debug;

use crate::utils::pathbuf::AppendExtension;

use super::BioinformaticsFileFormat;

//==================================//
// Binary Alignment Map (BAM) files //
//==================================//

/// Attempts to open a BAM file from a given source. Note that this file is private
/// because it should never be called by an external module (use [`open_and_parse`]
/// instead).
fn open<P>(src: P) -> anyhow::Result<bam::Reader<bgzf::Reader<BufReader<File>>>>
where
    P: AsRef<Path>,
{
    let path = src.as_ref();
    let file = File::open(path);

    match BioinformaticsFileFormat::try_detect(path) {
        Some(BioinformaticsFileFormat::BAM) => {
            let reader = file
                .map(BufReader::new)
                .with_context(|| "opening src BAM file")?;
            Ok(bam::Reader::new(reader))
        }
        Some(format) => bail!("incompatible formats: required BAM, found {}", format),
        None => {
            let ext = path
                .extension()
                .expect("file extension to exist")
                .to_str()
                .expect("extension to be convertible to &str");
            bail!("Not able to determine filetype for extension: {}", ext)
        }
    }
}

/// Utility struct which contains both the raw header (before being parsed, as a
/// string), and the parsed header (as a `noodles::sam::Header`).
pub struct RawAndParsedHeaders {
    /// The raw, unprocessed (and uncorrected) header [`String`] from the file.
    pub raw: String,
    /// The parsed, processed (and corrected) header [`Header`] from the file.
    pub parsed: Header,
}

/// Contains the BAM file reader, the parsed header from the BAM file, the reference
/// sequences read from the BAM file, and the location of the BAI index file.
pub struct ParsedBAMFile {
    /// A reader for the BAM file.
    pub reader: bam::Reader<bgzf::Reader<BufReader<File>>>,

    /// The raw and processed headers from the file, packaged together for convenience.
    pub header: RawAndParsedHeaders,

    /// The reference sequences read from the BAM file.
    pub reference_sequences: IndexMap<String, Map<ReferenceSequence>>,

    /// The path to the associated BAM index file.
    pub index_path: PathBuf,
}

#[derive(PartialEq, Eq)]
/// Utility enum to formalize in types whether or not to check for an index.
pub enum IndexCheck {
    /// Checks for an index file when opening and parsing the BAM file.
    CheckForIndex,

    /// _Does not_ check for an index file when opening and parsing the BAM file.
    DontCheckForIndex,
}

/// Opens and subsequently parses a BAM file's header. This is useful when opening BAM
/// files when you want the corrections applied by
/// [`super::sam::correct_common_header_mistakes`] to apply.
pub fn open_and_parse<P>(src: P, ensure_index: IndexCheck) -> anyhow::Result<ParsedBAMFile>
where
    P: AsRef<Path>,
{
    // (1) Construct the reader.
    debug!("reading BAM file from disk");
    let mut reader = open(&src)?;

    // (2) Though a BAM index is not always needed, we want, as good practice, to ensure
    // one exists for the times that it is needed (and also to deduplicate code where
    // we're checking for it). Thus, we will ask the user of this API to explicitly opt
    // out if they don't want to check for it.
    let index_path = src
        .as_ref()
        .to_path_buf()
        .append_extension("bai")
        .with_context(|| "constructing the BAM index filepath")?;

    if ensure_index == IndexCheck::CheckForIndex {
        debug!("checking that associated index exists for BAM file");
        bai::read(&index_path).with_context(|| "BAM index")?;
    }

    // (3) Parse the header and reference sequences.
    debug!("parsing the header and reference sequences");
    let raw_header = reader.read_header().with_context(|| "reading header")?;
    let parsed_header = super::sam::parse_header(raw_header.clone());
    let reference_sequences = reader
        .read_reference_sequences()
        .with_context(|| "reading reference sequences")?;

    // (4) Return the result.
    Ok(ParsedBAMFile {
        reader,
        header: RawAndParsedHeaders {
            raw: raw_header,
            parsed: parsed_header,
        },
        reference_sequences,
        index_path,
    })
}
