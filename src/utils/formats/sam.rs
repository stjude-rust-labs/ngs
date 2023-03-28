//! Utilities related to opening and manipulating Sequence Alignment Map (SAM) files.

use std::io::BufRead;
use std::path::Path;

use anyhow::bail;
use anyhow::Context;
use noodles::sam;
use tracing::debug;

use crate::utils::formats::utils::RawAndParsedHeaders;

use super::BioinformaticsFileFormat;

//====================================//
// Sequence Alignment Map (SAM) files //
//====================================//

/// Attempts to open a SAM file from a given source. Note that this file is private
/// because it should never be called by an external module (use [`open_and_parse`]
/// instead).
fn open<P>(src: P) -> anyhow::Result<sam::Reader<Box<dyn BufRead>>>
where
    P: AsRef<Path>,
{
    let path = src.as_ref();
    let file = std::fs::File::open(path);

    match BioinformaticsFileFormat::try_detect(path) {
        Some(BioinformaticsFileFormat::SAM) => {
            let reader = file.map(std::io::BufReader::new)?;
            Ok(sam::Reader::new(Box::new(reader)))
        }
        Some(format) => bail!("incompatible formats: required SAM, found {}", format),
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

/// Contains both the opened SAM file and the parsed header from the SAM file.
pub struct ParsedSAMFile {
    /// The reader for the SAM file.
    pub reader: sam::Reader<Box<dyn BufRead>>,

    /// The raw and parsed headers for the SAM file.
    ///
    /// Note that the raw header will be exactly as it is within the SAM file, and the
    /// parsed header will contain common corrections applied by
    /// [`correct_common_header_mistakes`].
    pub header: RawAndParsedHeaders,
}

/// Opens and subsequently parses a SAM file's header. This is useful when opening SAM
/// files when you want the corrections applied by [`correct_common_header_mistakes`] to
/// apply.
pub fn open_and_parse<P>(src: P) -> anyhow::Result<ParsedSAMFile>
where
    P: AsRef<Path>,
{
    // (1) Construct the reader.
    debug!("reading SAM file from disk");
    let mut reader = open(src)?;

    // (2) Parse the header.
    debug!("parsing the header");
    let raw_header = reader.read_header()?;
    let parsed_header = raw_header.parse().with_context(|| "parsing SAM header")?;

    // (3) Return the result.
    Ok(ParsedSAMFile {
        reader,
        header: RawAndParsedHeaders {
            raw: raw_header,
            parsed: parsed_header,
        },
    })
}

//==========================================//
// Async Sequence Alignment Map (SAM) files //
//==========================================//

/// Attempts to open a SAM file from a given source in an asyncronous fashion. Note that
/// this file is private because it should never be called by an external module (use
/// [`open_and_parse_async`] instead).
async fn open_async<P>(
    src: P,
) -> anyhow::Result<sam::AsyncReader<tokio::io::BufReader<tokio::fs::File>>>
where
    P: AsRef<Path>,
{
    let path = src.as_ref();
    let file = tokio::fs::File::open(path).await;

    match BioinformaticsFileFormat::try_detect(path) {
        Some(BioinformaticsFileFormat::SAM) => {
            let reader = file
                .map(tokio::io::BufReader::new)
                .map(sam::AsyncReader::new)?;
            Ok(reader)
        }
        Some(format) => bail!("incompatible formats: required SAM, found {}", format),
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

/// Contains both the opened async SAM file and the parsed header from the SAM file.
pub struct ParsedAsyncSAMFile {
    /// The reader for the SAM file.
    pub reader: sam::AsyncReader<tokio::io::BufReader<tokio::fs::File>>,

    /// The raw and parsed headers for the SAM file.
    ///
    /// Note that the raw header will be exactly as it is within the SAM file, and the
    /// parsed header will contain common corrections applied by
    /// [`correct_common_header_mistakes`].
    pub header: RawAndParsedHeaders,
}

/// Opens and subsequently parses a SAM file in an asyncronous manner. This is useful
/// when opening SAM files when you want the corrections applied by
/// [`correct_common_header_mistakes`] to apply.
pub async fn open_and_parse_async<P>(src: P) -> anyhow::Result<ParsedAsyncSAMFile>
where
    P: AsRef<Path>,
{
    // (1) Construct the reader.
    debug!("reading SAM file from disk");
    let mut reader = open_async(src).await?;

    // (2) Parse the header.
    debug!("parsing the header");
    let raw_header = reader.read_header().await?;
    let parsed_header = raw_header.parse().with_context(|| "parsing SAM header")?;

    // (3) Return the result.
    Ok(ParsedAsyncSAMFile {
        reader,
        header: RawAndParsedHeaders {
            raw: raw_header,
            parsed: parsed_header,
        },
    })
}
