//! Utilities related to opening and manipulating CRAM files.

use std::io::BufReader;
use std::path::Path;
use std::path::PathBuf;

use anyhow::bail;
use anyhow::Context;
use noodles::cram;
use noodles::cram::crai;
use noodles::cram::FileDefinition;
use tracing::debug;

use crate::utils::formats::utils::IndexCheck;
use crate::utils::formats::utils::RawAndParsedHeaders;
use crate::utils::formats::BioinformaticsFileFormat;
use crate::utils::pathbuf::AppendExtension;

//============//
// CRAM files //
//============//

/// Attempts to open a CRAM file from a given source. Note that this file is private
/// because it should never be called by an external module (use [`open_and_parse`]
/// instead).
fn open<P>(src: P) -> anyhow::Result<cram::Reader<BufReader<std::fs::File>>>
where
    P: AsRef<Path>,
{
    let path = src.as_ref();
    let file = std::fs::File::open(path);

    match BioinformaticsFileFormat::try_detect(path) {
        Some(BioinformaticsFileFormat::CRAM) => {
            let reader = file
                .map(BufReader::new)
                .with_context(|| "opening CRAM file")?;
            Ok(cram::Reader::new(reader))
        }
        Some(format) => bail!("incompatible formats: required CRAM, found {}", format),
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

/// Contains the CRAM file reader, the parsed header from the CRAM file, the file
/// definition for the CRAM file, and the location of the CRAI index file.
pub struct ParsedCRAMFile {
    /// A reader for the CRAM file.
    pub reader: cram::Reader<BufReader<std::fs::File>>,

    /// The raw and processed headers from the file, packaged together for convenience.
    pub header: RawAndParsedHeaders,

    /// File definition
    pub file_definition: FileDefinition,

    /// The path to the associated CRAM index file.
    pub index_path: PathBuf,
}

/// Opens and subsequently parses a CRAM file. This is useful when opening CRAM
/// files when you want the corrections applied by
/// [`super::sam::correct_common_header_mistakes`] to apply.
pub fn open_and_parse<P>(src: P, ensure_index: IndexCheck) -> anyhow::Result<ParsedCRAMFile>
where
    P: AsRef<Path>,
{
    // (1) Construct the reader.
    debug!("reading CRAM file from disk");
    let mut reader = open(&src)?;

    // (2) Checks for an index (or doesn't) based on the user's input.
    let index_path = src
        .as_ref()
        .to_path_buf()
        .append_extension("crai")
        .with_context(|| "constructing CRAM index filepath")?;

    match ensure_index {
        IndexCheck::Full => {
            debug!("checking the index's header and contents");
            crai::read(&index_path).with_context(|| "reading CRAM index")?;
        }
        IndexCheck::HeaderOnly => {
            // A header doesn't really exist for a CRAM file, so this option doesn't
            // make sense in the context of CRAM.
            unimplemented!()
        }
        IndexCheck::None => {}
    }

    // (3) Parse the file definition and the header.
    debug!("parsing the file definition and header");
    let file_definition = reader.read_file_definition()?;
    let raw_header = reader
        .read_file_header()
        .with_context(|| "reading CRAM header")?;
    let parsed_header = raw_header.parse().with_context(|| "parsing CRAM header")?;

    // (4) Return the result.
    Ok(ParsedCRAMFile {
        reader,
        header: RawAndParsedHeaders {
            raw: raw_header,
            parsed: parsed_header,
        },
        file_definition,
        index_path,
    })
}

//==================//
// Async CRAM files //
//==================//

/// Attempts to open a CRAM file from a given source using an `async` function. Note that
/// this file is private because it should never be called by an external module (use
/// [`open_and_parse_async`] instead).
async fn open_async<P>(src: P) -> anyhow::Result<cram::AsyncReader<tokio::fs::File>>
where
    P: AsRef<Path>,
{
    let path = src.as_ref();
    let file = tokio::fs::File::open(path).await;

    match BioinformaticsFileFormat::try_detect(path) {
        Some(BioinformaticsFileFormat::CRAM) => {
            let reader = file
                .map(cram::AsyncReader::new)
                .with_context(|| "opening CRAM file")?;
            Ok(reader)
        }
        Some(format) => bail!("incompatible formats: required CRAM, found {}", format),
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

/// Contains the async CRAM file reader, the parsed header from the CRAM file, the file
/// definition for the CRAM file, and the location of the CRAI index file.
pub struct ParsedAsyncCRAMFile {
    /// An async reader for the CRAM file.
    pub reader: cram::AsyncReader<tokio::fs::File>,

    /// The raw and processed headers from the file, packaged together for convenience.
    pub header: RawAndParsedHeaders,

    /// File definition
    pub file_definition: FileDefinition,

    /// The path to the associated CRAM index file.
    pub index_path: PathBuf,
}

/// Opens and subsequently parses a BAM file in an async fashion. This is useful when
/// opening BAM files when you want the corrections applied by
/// [`super::sam::correct_common_header_mistakes`] to apply.
pub async fn open_and_parse_async<P>(
    src: P,
    ensure_index: IndexCheck,
) -> anyhow::Result<ParsedAsyncCRAMFile>
where
    P: AsRef<Path>,
{
    // (1) Construct the reader.
    debug!("reading CRAM file from disk");
    let mut reader = open_async(&src).await?;

    // (2) Checks for an index (or doesn't) based on the user's input.
    let index_path = src
        .as_ref()
        .to_path_buf()
        .append_extension("crai")
        .with_context(|| "constructing CRAM index filepath")?;

    match ensure_index {
        IndexCheck::Full => {
            debug!("checking the index's header and contents");
            crai::read(&index_path).with_context(|| "reading CRAM index")?;
        }
        IndexCheck::HeaderOnly => {
            // A header doesn't really exist for a CRAM file, so this option doesn't
            // make sense in the context of CRAM.
            unimplemented!()
        }
        IndexCheck::None => {}
    }

    // (3) Parse the file definition and the header.
    debug!("parsing the file definition and header");
    let file_definition = reader.read_file_definition().await?;
    let raw_header = reader
        .read_file_header()
        .await
        .with_context(|| "reading CRAM header")?;
    let parsed_header = raw_header.parse().with_context(|| "parsing CRAM header")?;

    // (4) Return the result.
    Ok(ParsedAsyncCRAMFile {
        reader,
        header: RawAndParsedHeaders {
            raw: raw_header,
            parsed: parsed_header,
        },
        file_definition,
        index_path,
    })
}
