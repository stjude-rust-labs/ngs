//! FASTA indexing
//!
//! NOTICE: this was taken almost verbatim from @zaeleus's excellent example in
//! noodles. You can find that source code in the `fasta_index.rs` example of the
//! `noodles_fasta` crate at commit 44b7d19.

use crate::utils::path::AppendExtension;
use std::{fs::File, path::PathBuf};

use anyhow::{bail, Context};
use noodles::fasta::{self, fai};

/// Indexes a FASTA file stored at `src`.
pub fn index(src: PathBuf) -> anyhow::Result<()> {
    let index = fasta::index(&src).with_context(|| "building FASTA index")?;
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

    let mut writer = File::create(fai)
        .map(fai::Writer::new)
        .with_context(|| "initializing the new FASTA index file")?;
    writer
        .write_index(&index)
        .with_context(|| "writing the FASTA index file")?;

    Ok(())
}
