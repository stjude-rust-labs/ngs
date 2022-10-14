//! Utilities related to opening and manipulating FASTA files.

use std::{
    fs::File,
    io::{BufRead, BufReader, Seek},
    path::{Path, PathBuf},
};

use anyhow::bail;
use noodles::fasta::{self, repository::adapters::IndexedReader};

use crate::utils::pathbuf::AppendExtension;

use super::BioinformaticsFileFormat;

mod phred;

trait BufReadSeek: BufRead + Seek {}

/// Attempts to open a FASTA file from a given source.
pub fn open<P>(src: P) -> anyhow::Result<fasta::Reader<Box<dyn BufRead>>>
where
    P: AsRef<Path>,
{
    let path = src.as_ref();
    let file = File::open(path);

    match BioinformaticsFileFormat::try_detect(path) {
        Some(BioinformaticsFileFormat::FASTA_GZ) => bail!(
            "This command does not yet support gzipped FASTA \
            files. Please unzip your FASTA file and try again."
        ),
        Some(BioinformaticsFileFormat::FASTA) => {
            let reader = fasta::reader::Builder::default()
                .build_from_path(src)
                .map(Box::new)?;
            Ok(reader)
        }
        Some(format) => bail!("incompatible formats: required FASTA, found {}", format),
        None => {
            let ext = path
                .extension()
                .expect("file extension to exist")
                .to_str()
                .expect("extension to be convertible to &str");
            bail!("Not able to determine filetype for extension: {}", ext)
        }
    }
}

pub fn open_repository(src: PathBuf) -> anyhow::Result<()> {
    // TODO: remove in a future version when noodles gives an error message that
    // suggests you should index your FASTA file (as of the time of writing, it
    // just gives an error back saying "unsupported" if you don't have an
    // index).

    // (1) Open up the FASTA file and map it to a Repository.
    let fasta = open(&src)?;

    // (2) Check to make sure that the FASTA index actually exists. If it does
    // not then error out.
    let fai_filepath = src.append_extension("fai")?;
    if !fai_filepath.exists() {
        bail!(
            "couldn't find an index for your reference FASTA: is the FASTA indexed? \
        Run `ngs index [FASTA]` to index the FASTA file."
        )
    }

    let l = IndexedReader::new(fasta);
    // .map(fasta::Repository::new)
    // .with_context(|| "building FASTA repository")?;

    Ok(())
}
