//! SAM viewing

use std::io;
use std::io::Write;
use std::path::PathBuf;

use anyhow::bail;
use anyhow::Context;
use noodles::sam;
use noodles::sam::AlignmentWriter;

use crate::utils::formats::sam::ParsedSAMFile;
use crate::view::command::Mode;

/// Main method for SAM viewing.
pub fn view(src: PathBuf, query: Option<String>, mode: Mode) -> anyhow::Result<()> {
    // (1) Check if the user provided a query. If they did, we do not support any sort
    //     tabix-like indexing and we would highly recommend the user take advantage of
    //     BAM/CRAM files. If anyone stumbles across this comment and sees a reason we
    //     _should_ support it, please file an issue.
    if query.is_some() {
        bail!(
            "querying is not supported for SAM files. Please convert your SAM file \
             to a BAM/CRAM file and then index it you'd like to query a region within \
             the file."
        )
    }

    // (2) Opens and parses the SAM file.
    let ParsedSAMFile {
        mut reader, header, ..
    } = crate::utils::formats::sam::open_and_parse(&src)?;

    // (3) Determine the handle with which to write the output. TODO: for now, we just
    // default to stdout, but in the future we will support writing to another path.
    let stdout = io::stdout();
    let mut handle = stdout.lock();

    // (4) If the user specified to output the header, output the header.
    if mode == Mode::Full || mode == Mode::HeaderOnly {
        write!(handle, "{}", header.raw).with_context(|| "writing header to stream")?;
    }

    // (5) If the mode is header only, nothing left to do, so return.
    if mode == Mode::HeaderOnly {
        return Ok(());
    }

    // (6) Writes the records to the output stream.
    let mut writer = sam::Writer::new(handle);

    for result in reader.records(&header.parsed) {
        let record = result?;
        writer.write_alignment_record(&header.parsed, &record)?;
    }

    writer.finish(&header.parsed)?;
    Ok(())
}
