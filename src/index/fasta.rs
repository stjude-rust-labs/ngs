//! FASTA indexing
//!
//! NOTICE: this was taken almost verbatim from @zaeleus's excellent example in
//! noodles. You can find that source code in the `fasta_index.rs` example of the
//! `noodles_fasta` crate at commit 44b7d19.

use crate::utils::pathbuf::AppendExtension;
use std::{fs::File, path::PathBuf};

use anyhow::{bail, Context};
use noodles::fasta::{self, fai};
use tracing::{debug, info};

/// Indexes a FASTA file stored at `src`.
pub fn index(src: PathBuf) -> anyhow::Result<()> {
    info!("indexing FASTA file");

    // (1) Reads the file from disk.
    debug!("reading FASTA file from disk");
    let index = fasta::index(&src).with_context(|| "building FASTA index")?;

    // (2) Calculate where the FASTA index should go and check if a file is
    // already there. Error out if so.
    let fai = src
        .append_extension("fai")
        .with_context(|| "building filepath for FASTA index")?;

    if fai.exists() {
        bail!(
            "refusing to overwrite existing index file: {}. Please delete \
                and rerun if you'd like to replace it.",
            fai.display()
        );
    }

    // (3) Write the index to disk.
    debug!("writing the FASTA index to disk");
    let mut writer = File::create(fai)
        .map(fai::Writer::new)
        .with_context(|| "initializing the new FASTA index file")?;
    writer
        .write_index(&index)
        .with_context(|| "writing the FASTA index file")?;

    Ok(())
}
