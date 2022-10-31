//! CRAM viewing

use std::path::PathBuf;

use anyhow::bail;
use anyhow::Context;
use futures::TryStreamExt;
use noodles::fasta;
use noodles::fasta::repository::adapters::IndexedReader;
use noodles::sam;
use tokio::io;
use tokio::io::AsyncWriteExt;

use crate::utils::formats::cram::ParsedAsyncCRAMFile;
use crate::utils::formats::utils::IndexCheck;
use crate::utils::pathbuf::AppendExtension;
use crate::view::command::Mode;

/// Main method for CRAM viewing in an asyncronous fashion.
pub async fn view(
    src: PathBuf,
    query: Option<String>,
    reference_fasta: PathBuf,
    mode: Mode,
) -> anyhow::Result<()> {
    // (1) Opens and parses the BAM file.
    let ParsedAsyncCRAMFile {
        mut reader,
        header,
        index_path,
        ..
    } = crate::utils::formats::cram::open_and_parse_async(&src, IndexCheck::Full).await?;

    // (2) Determine the handle with which to write the output. TODO: for now, we just
    // default to stdout, but in the future we will support writing to another path.
    let mut handle = io::stdout();

    // (3) Build the FASTA repository.
    let repository = fasta::indexed_reader::Builder::default()
        .build_from_path(&reference_fasta)
        .map(IndexedReader::new)
        .map(fasta::Repository::new)
        .with_context(|| "building FASTA repository")?;

    // TODO: remove in a future version when noodles gives an error message that
    // suggests you should index your FASTA file (as of the time of writing, it
    // just gives an error back saying "unsupported" if you don't have an
    // index).
    let fai_filepath = reference_fasta.append_extension("fai")?;
    if !fai_filepath.exists() {
        bail!(
            "couldn't find an index for your reference FASTA: is the FASTA indexed? \
        Run `ngs index [FASTA]` to index the FASTA file."
        )
    }

    // (4) If the user specified to output the header, output the header.
    if mode == Mode::Full || mode == Mode::HeaderOnly {
        handle
            .write_all(header.raw.to_string().as_bytes())
            .await
            .with_context(|| "writing header to stream")?;
    }

    // (6) If the mode is header-only, nothing left to do, so return.
    if mode == Mode::HeaderOnly {
        return Ok(());
    }

    // (8) Writes the records to the output stream.
    let mut writer = sam::AsyncWriter::new(handle);

    // if let Some(query) = query {
    //     // (a) If a query is specified, print just the records that fall within the query.
    //     let index = crai::read(index_path).with_context(|| "reading CRAM index")?;
    //     let region = query.parse().with_context(|| "parsing query")?;

    //     let records = reader
    //         .query(&repository, &header, &index, &region)
    //         .with_context(|| "querying CRAM file")?;

    //     for result in records {
    //         let record = result
    //             .and_then(|record| record.try_into_alignment_record(&header))
    //             .with_context(|| "reading record")?;
    //         writer
    //             .write_alignment_record(&header.parsed, &record)
    //             .await?;
    //     }
    // } else {
    // (b) Else, print all of the records in the file.

    // TODO: use `read_record` when this issue is resolved:
    // https://github.com/zaeleus/noodles/issues/126.
    let mut records = reader.records(&repository, &header.parsed);

    while let Some(record) = records.try_next().await? {
        let record = record.try_into_alignment_record(&header.parsed)?;
        writer.write_record(&header.parsed, &record).await?;
    }

    Ok(())
}
