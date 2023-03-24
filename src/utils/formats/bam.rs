//! Utilities related to opening and manipulating Binary Alignment Map (BAM) files.

use std::fs::File;
use std::io::BufReader;
use std::path::Path;
use std::path::PathBuf;

use anyhow::bail;
use anyhow::Context;
use indexmap::IndexMap;
use noodles::bam;
use noodles::bam::bai;
use noodles::bgzf;
use noodles::sam::header::record::value::map::Map;
use noodles::sam::header::record::value::map::ReferenceSequence;
use noodles::sam::record::ReferenceSequenceName;
use tracing::debug;

use crate::utils::formats::utils::IndexCheck;
use crate::utils::formats::utils::RawAndParsedHeaders;
use crate::utils::pathbuf::AppendExtension;

use super::BioinformaticsFileFormat;

//==================================//
// Binary Alignment Map (BAM) files //
//==================================//

/// Attempts to open a BAM file from a given source. Note that this file is private
/// because it should never be called by an external module (use [`open_and_parse`]
/// instead).
fn open<P>(src: P) -> anyhow::Result<bam::Reader<bgzf::Reader<BufReader<std::fs::File>>>>
where
    P: AsRef<Path>,
{
    let path = src.as_ref();
    let file = std::fs::File::open(path);

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

/// Contains the BAM file reader, the parsed header from the BAM file, the reference
/// sequences read from the BAM file, and the location of the BAI index file.
pub struct ParsedBAMFile {
    /// A reader for the BAM file.
    pub reader: bam::Reader<bgzf::Reader<BufReader<std::fs::File>>>,

    /// The raw and processed headers from the file, packaged together for convenience.
    pub header: RawAndParsedHeaders,

    /// The reference sequences read from the BAM file.
    pub reference_sequences: IndexMap<ReferenceSequenceName, Map<ReferenceSequence>>,

    /// The path to the associated BAM index file.
    pub index_path: PathBuf,
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

    // (2) Checks for an index (or doesn't) based on the user's input.
    let index_path = src
        .as_ref()
        .to_path_buf()
        .append_extension("bai")
        .with_context(|| "constructing the BAM index filepath")?;

    match ensure_index {
        IndexCheck::Full => {
            debug!("checking the index's header and contents");
            bai::read(&index_path).with_context(|| "BAM index")?;
        }
        IndexCheck::HeaderOnly => {
            debug!("checking the index's header");
            let mut reader = File::open(&index_path).map(bai::Reader::new)?;
            reader.read_header()?;
        }
        IndexCheck::None => {}
    }

    // (3) Parse the header and reference sequences.
    debug!("parsing the header and reference sequences");
    let raw_header = reader.read_header().with_context(|| "reading header")?;
    let parsed_header =
        super::sam::parse_header(raw_header.clone()).with_context(|| "parsing header")?;
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

//========================================//
// Async Binary Alignment Map (BAM) files //
//========================================//

/// Attempts to open a BAM file from a given source using an `async` function. Note that
/// this file is private because it should never be called by an external module (use
/// [`open_and_parse_async`] instead).
async fn open_async<P>(
    src: P,
) -> anyhow::Result<bam::AsyncReader<bgzf::AsyncReader<tokio::fs::File>>>
where
    P: AsRef<Path>,
{
    let path = src.as_ref();
    let file = tokio::fs::File::open(path).await;

    match BioinformaticsFileFormat::try_detect(path) {
        Some(BioinformaticsFileFormat::BAM) => {
            let reader = file
                .map(bam::AsyncReader::new)
                .with_context(|| "opening src BAM file")?;
            Ok(reader)
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

/// Contains the async BAM file reader, the parsed header from the BAM file, the
/// reference sequences read from the BAM file, and the location of the BAI index file.
pub struct ParsedAsyncBAMFile {
    /// An async reader for the BAM file.
    pub reader: bam::AsyncReader<bgzf::AsyncReader<tokio::fs::File>>,

    /// The raw and processed headers from the file, packaged together for convenience.
    pub header: RawAndParsedHeaders,

    /// The reference sequences read from the BAM file.
    pub reference_sequences: IndexMap<ReferenceSequenceName, Map<ReferenceSequence>>,

    /// The path to the associated BAM index file.
    pub index_path: PathBuf,
}

/// Opens and subsequently parses a BAM file in an async fashion. This is useful when
/// opening BAM files when you want the corrections applied by
/// [`super::sam::correct_common_header_mistakes`] to apply.
pub async fn open_and_parse_async<P>(
    src: P,
    ensure_index: IndexCheck,
) -> anyhow::Result<ParsedAsyncBAMFile>
where
    P: AsRef<Path>,
{
    // (1) Construct the reader.
    debug!("reading BAM file from disk");
    let mut reader = open_async(&src).await?;

    // (2) Though a BAM index is not always needed, we want, as good practice, to ensure
    // one exists for the times that it is needed (and also to deduplicate code where
    // we're checking for it). Thus, we will ask the user of this API to explicitly opt
    // out if they don't want to check for it.
    let index_path = src
        .as_ref()
        .to_path_buf()
        .append_extension("bai")
        .with_context(|| "constructing the BAM index filepath")?;

    if ensure_index == IndexCheck::Full {
        debug!("checking that associated index exists for BAM file");
        bai::read(&index_path).with_context(|| "BAM index")?;
    }

    // (3) Parse the header and reference sequences.
    debug!("parsing the header and reference sequences");
    let raw_header = reader
        .read_header()
        .await
        .with_context(|| "reading header")?;
    let parsed_header =
        super::sam::parse_header(raw_header.clone()).with_context(|| "parsing header")?;
    let reference_sequences = reader
        .read_reference_sequences()
        .await
        .with_context(|| "reading reference sequences")?;

    // (4) Return the result.
    Ok(ParsedAsyncBAMFile {
        reader,
        header: RawAndParsedHeaders {
            raw: raw_header,
            parsed: parsed_header,
        },
        reference_sequences,
        index_path,
    })
}
