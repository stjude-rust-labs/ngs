//! Utilities related to opening and manipulating FASTA files.

use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use anyhow::bail;
use noodles::fasta;

use super::BioinformaticsFileFormat;

mod phred;

/// Attempts to open a FASTA file from a given source.
pub fn open<P>(src: P) -> anyhow::Result<fasta::Reader<BufReader<File>>>
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
            Ok(file.map(BufReader::new).map(fasta::Reader::new)?)
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
