//! GFF viewing

use std::io;
use std::path::PathBuf;

use anyhow::Context;
use noodles::gff;

/// Main method for GFF viewing.
pub fn view(src: PathBuf) -> anyhow::Result<()> {
    // (1) Opens and parses the GFF file.
    let mut gff =
        crate::utils::formats::gff::open(src).with_context(|| "reading GFF input file")?;

    // (2) Determine the handle with which to write the output. TODO: for now, we just
    // default to stdout, but in the future we will support writing to another path.
    let stdout = io::stdout();
    let handle = stdout.lock();

    // (3) Writes the lines to the stream.
    let mut writer = gff::Writer::new(handle);
    for result in gff.lines() {
        let line = result.with_context(|| "reading GFF line")?;
        writer
            .write_line(&line)
            .with_context(|| "writing GFF line to stream")?;
    }

    Ok(())
}
