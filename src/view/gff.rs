//! GFF viewing

use std::io;
use std::path::PathBuf;

use noodles::gff;

/// Main method for GFF viewing.
pub fn view(src: PathBuf) -> anyhow::Result<()> {
    // (1) Opens and parses the GFF file.
    let mut gff = crate::utils::formats::gff::open(src)?;

    // (2) Determine the handle with which to write the output. TODO: for now, we just
    // default to stdout, but in the future we will support writing to another path.
    let stdout = io::stdout();
    let handle = stdout.lock();

    // (3) Writes the lines to the stream.
    let mut writer = gff::Writer::new(handle);
    for result in gff.lines() {
        let line = result?;
        writer.write_line(&line)?;
    }

    Ok(())
}
