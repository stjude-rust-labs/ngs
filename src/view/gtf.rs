//! GTF viewing

use std::io;
use std::path::PathBuf;

use anyhow::Context;
use noodles::gtf;

/// Main method for GTF viewing.
pub fn view(src: PathBuf) -> anyhow::Result<()> {
    // (1) Opens and parses the GTF file.
    let mut gtf =
        crate::utils::formats::gtf::open(src).with_context(|| "reading GTF input file")?;

    // (2) Determine the handle with which to write the output. TODO: for now, we just
    // default to stdout, but in the future we will support writing to another path.
    let stdout = io::stdout();
    let handle = stdout.lock();

    // (3) Writes the lines to the stream.
    let mut writer = gtf::Writer::new(handle);
    for result in gtf.lines() {
        let line = result.with_context(|| "reading GTF line")?;
        writer
            .write_line(&line)
            .with_context(|| "writing GTF line to stream")?;
    }

    Ok(())
}
