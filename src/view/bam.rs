//! BAM viewing

use std::path::PathBuf;

use anyhow::Context;
use futures::TryStreamExt;
use noodles::bam::bai;
use noodles::sam;
use tokio::io::{self, AsyncWriteExt};

use crate::utils::formats::bam::ParsedAsyncBAMFile;
use crate::utils::formats::utils::IndexCheck;
use crate::view::command::Mode;

/// Main method for BAM viewing.
pub async fn view(src: PathBuf, query: Option<String>, mode: Mode) -> anyhow::Result<()> {
    // (1) Opens and parses the BAM file.
    let ParsedAsyncBAMFile {
        mut reader,
        header,
        index_path: bai_path,
        ..
    } = crate::utils::formats::bam::open_and_parse_async(&src, IndexCheck::HeaderOnly).await?;

    // (2) Determine the handle with which to write the output. TODO: for now, we just
    // default to stdout, but in the future we will support writing to another path.
    let mut handle = io::stdout();

    // (3) If the user specified to output the header, output the header.
    if mode == Mode::Full || mode == Mode::HeaderOnly {
        handle
            .write_all(header.raw.to_string().as_bytes())
            .await
            .with_context(|| "writing header to stream")?;
    }

    // (4) If the mode is header only, nothing left to do, so return.
    if mode == Mode::HeaderOnly {
        return Ok(());
    }

    // (5) Writes the records to the output stream.
    let mut writer = sam::AsyncWriter::new(handle);

    if let Some(query) = query {
        // (a) If a query is specified, print just the records that fall within the query.
        let index = bai::read(bai_path).with_context(|| "reading BAM index")?;
        let region = query.parse().with_context(|| "parsing query")?;

        let mut records = reader
            .query(&header.parsed, &index, &region)
            .with_context(|| "querying BAM file")?;

        while let Some(record) = records.try_next().await? {
            writer
                .write_alignment_record(&header.parsed, &record)
                .await?;
        }
    } else {
        // (b) Else, print all of the records in the file.
        let mut records = reader.records(&header.parsed);
        while let Some(record) = records.try_next().await? {
            writer
                .write_alignment_record(&header.parsed, &record)
                .await?;
        }
    }

    Ok(())
}
