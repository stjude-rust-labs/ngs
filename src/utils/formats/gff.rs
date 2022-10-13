//! Utilities related to opening and manipulating GFF files.

use anyhow::bail;
use flate2::read::MultiGzDecoder;
use noodles::gff;
use std::{
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
};

use super::BioinformaticsFileFormat;

/// Attempts to open a GFF file from a given source.
pub fn open<P>(src: P) -> anyhow::Result<gff::Reader<Box<dyn BufRead>>>
where
    P: AsRef<Path>,
{
    let path = src.as_ref();
    let file = File::open(path);

    match BioinformaticsFileFormat::try_detect(path) {
        Some(BioinformaticsFileFormat::GFF_GZ) => {
            let reader = file.map(MultiGzDecoder::new).map(BufReader::new)?;
            Ok(gff::Reader::new(Box::new(reader)))
        }
        Some(BioinformaticsFileFormat::GFF) => {
            let reader = file.map(BufReader::new)?;
            Ok(gff::Reader::new(Box::new(reader)))
        }
        Some(format) => bail!("incompatible formats: required GFF, found {}", format),
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
