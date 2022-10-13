//! BAM viewing

use std::{
    fs::File,
    io::{self as io, Write},
    path::PathBuf,
};

use anyhow::Context;
use noodles::{
    bam::{self, bai},
    sam::{self, AlignmentWriter},
};
use tracing::debug;

use crate::utils::formats::sam::parse_header;

/// Main method for BAM viewing.
pub fn view(src: PathBuf, query: Option<String>, show_header: bool) -> anyhow::Result<()> {
    // (1) Reads the file from disk.
    debug!("reading BAM file from disk");
    let mut reader = File::open(&src)
        .map(bam::Reader::new)
        .with_context(|| "opening src file")?;

    // (2) Determine the handle with which to write the output.
    let stdout = io::stdout();
    let mut handle = stdout.lock();

    // (3) If the user specified to output the header, output the header.
    let ht = reader.read_header().with_context(|| "reading BAM header")?;
    if show_header {
        write!(handle, "{}", ht).with_context(|| "writing header to stream")?;
    }

    // (4) Parse the header and reference sequences.
    let header = parse_header(ht);
    reader
        .read_reference_sequences()
        .with_context(|| "reading reference sequences")?;

    let mut writer = sam::Writer::new(handle);
    if let Some(query) = query {
        // (5a) If a query is specified, print just the records that fall within the query.
        let index =
            bai::read(src.with_extension("bam.bai")).with_context(|| "reading BAM index")?;
        let region = query.parse().with_context(|| "parsing query")?;

        let records = reader
            .query(header.reference_sequences(), &index, &region)
            .with_context(|| "querying BAM file")?;

        for result in records {
            let record = result?;
            writer.write_alignment_record(&header, &record)?;
        }
    } else {
        // (5b) Else, print all of the records in the file.
        for result in reader.records() {
            let record = result?;
            writer.write_alignment_record(&header, &record)?;
        }
    }

    writer.finish(&header)?;
    Ok(())
}
