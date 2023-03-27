//! CRAM viewing

use std::path::PathBuf;

use anyhow::bail;
use anyhow::Context;
use futures::TryStreamExt;
use noodles::cram;
use noodles::cram::crai;
use noodles::fasta;
use noodles::fasta::repository::adapters::IndexedReader;
use noodles::sam;
use tokio::io;
use tokio::io::AsyncWriteExt;
use tracing::debug;

use crate::utils::formats::sam::parse_header;
use crate::utils::pathbuf::AppendExtension;
use crate::view::command::Mode;

/// Main method for BAM viewing.
pub async fn view(
    src: PathBuf,
    query: Option<String>,
    reference_fasta: PathBuf,
    mode: Mode,
) -> anyhow::Result<()> {
    // (1) Reads the file from disk.
    debug!("reading CRAM file from disk");
    let mut reader = tokio::fs::File::open(&src)
        .await
        .map(cram::AsyncReader::new)
        .with_context(|| "opening CRAM input file")?;

    // (2) Determine the handle with which to write the output.
    let mut handle = io::stdout();

    // (3) Build the FASTA repository and associated index.
    let repository = fasta::indexed_reader::Builder::default()
        .build_from_path(&reference_fasta)
        .map(IndexedReader::new)
        .map(fasta::Repository::new)
        .with_context(|| "reading reference FASTA and associated index")?;

    // TODO: remove in a future version when noodles gives an error message that
    // suggests you should index your FASTA file (as of the time of writing, it
    // just gives an error back saying "unsupported" if you don't have an
    // index).
    let fai_filepath = reference_fasta
        .append_extension("fai")
        .with_context(|| "setting FAI extension")?;
    if !fai_filepath.exists() {
        bail!(
            "couldn't find an index for your reference FASTA: is the FASTA indexed? \
        Run `ngs index [FASTA]` to index the FASTA file."
        )
    }

    // (4) Read the file's definition.
    reader
        .read_file_definition()
        .await
        .with_context(|| "reading CRAM file definition")?;

    // (5) If the user specified to output the header, output the raw header (before
    // applying any corrections).
    let ht = reader
        .read_file_header()
        .await
        .with_context(|| "reading CRAM file header")?;

    if mode == Mode::Full || mode == Mode::HeaderOnly {
        handle
            .write_all(ht.as_bytes())
            .await
            .with_context(|| "writing header to stream")?;
    }

    // (6) If the mode is header-only, nothing left to do, so return.
    if mode == Mode::HeaderOnly {
        return Ok(());
    }

    // (7) Parses the header text.
    let header = parse_header(ht).with_context(|| "parsing CRAM header")?;

    // (8) Writes the records to the output stream.
    let mut writer = sam::AsyncWriter::new(handle);

    if let Some(query) = query {
        // (a) If a query is specified, print just the records that fall within the query.
        let index = crai::r#async::read(src.with_extension("cram.crai"))
            .await
            .with_context(|| "reading CRAM index")?;
        let region = query.parse().with_context(|| "parsing query")?;

        let mut records = reader
            .query(&repository, &header, &index, &region)
            .with_context(|| "querying CRAM file")?;

        while let Some(record) = records
            .try_next()
            .await
            .with_context(|| "reading CRAM record")?
        {
            let record = record
                .try_into_alignment_record(&header)
                .with_context(|| "parsing CRAM record")?;
            writer
                .write_alignment_record(&header, &record)
                .await
                .with_context(|| "writing record to stream")?;
        }
    } else {
        // (b) Else, print all of the records in the file.
        let mut records = reader.records(&repository, &header);

        while let Some(record) = records
            .try_next()
            .await
            .with_context(|| "reading CRAM record")?
        {
            let record = record
                .try_into_alignment_record(&header)
                .with_context(|| "parsing CRAM record")?;

            writer
                .write_alignment_record(&header, &record)
                .await
                .with_context(|| "writing record to stream")?;
        }
    }

    Ok(())
}
