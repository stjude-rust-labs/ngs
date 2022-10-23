//! BAM viewing

use std::io;
use std::io::Write;
use std::path::PathBuf;

use anyhow::Context;
use noodles::bam::bai;
use noodles::sam;
use noodles::sam::AlignmentWriter;

use crate::utils::formats::bam::ParsedBAMFile;
use crate::utils::formats::utils::IndexCheck;
use crate::view::command::Mode;

/// Main method for BAM viewing.
pub fn view(src: PathBuf, query: Option<String>, mode: Mode) -> anyhow::Result<()> {
    // (1) Opens and parses the BAM file.
    let ParsedBAMFile {
        mut reader,
        header,
        index_path: bai_path,
        ..
    } = crate::utils::formats::bam::open_and_parse(&src, IndexCheck::HeaderOnly)?;

    // (2) Determine the handle with which to write the output. TODO: for now, we just
    // default to stdout, but in the future we will support writing to another path.
    let stdout = io::stdout();
    let mut handle = stdout.lock();

    // (3) If the user specified to output the header, output the header.
    if mode == Mode::Full || mode == Mode::HeaderOnly {
        write!(handle, "{}", header.raw).with_context(|| "writing header to stream")?;
    }

    // (4) If the mode is header only, nothing left to do, so return.
    if mode == Mode::HeaderOnly {
        return Ok(());
    }

    // (5) Writes the records to the output stream.
    let mut writer = sam::Writer::new(handle);

    if let Some(query) = query {
        // (a) If a query is specified, print just the records that fall within the query.
        let index = bai::read(bai_path).with_context(|| "reading BAM index")?;
        let region = query.parse().with_context(|| "parsing query")?;

        let records = reader
            .query(header.parsed.reference_sequences(), &index, &region)
            .with_context(|| "querying BAM file")?;

        for result in records {
            let record = result?;
            writer.write_alignment_record(&header.parsed, &record)?;
        }
    } else {
        // (b) Else, print all of the records in the file.
        for result in reader.records() {
            let record = result?;
            writer.write_alignment_record(&header.parsed, &record)?;
        }
    }

    writer.finish(&header.parsed)?;
    Ok(())
}
