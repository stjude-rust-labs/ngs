//! CRAM indexing

use std::{fs::File, path::PathBuf};

use anyhow::{bail, Context};
use noodles::cram::{self as cram, crai};
use tracing::{debug, info};

use crate::utils::pathbuf::AppendExtension;

/// Main method for CRAM indexing.
pub fn index(src: PathBuf) -> anyhow::Result<()> {
    info!("indexing CRAM file");

    // (1) Reads the file from disk.
    debug!("reading CRAM file from disk");
    let index = cram::index(&src).with_context(|| "building CRAM index")?;

    // (2) Calculate where the CRAM index should go and check if a file is
    // already there. Error out if so.
    let crai = src
        .append_extension("crai")
        .with_context(|| "constructing CRAM index filepath")?;

    if crai.exists() {
        bail!(
            "refusing to overwrite existing index file: {}. Please delete \
                and rerun if you'd like to replace it.",
            crai.display()
        );
    }

    // (3) Write the index to disk.
    debug!("writing the CRAM index to disk");
    let mut writer = File::create(crai)
        .map(crai::Writer::new)
        .with_context(|| "creating CRAM index output file")?;

    writer.write_index(&index)?;

    Ok(())
}
