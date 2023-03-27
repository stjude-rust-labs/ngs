//! Utilities related to opening and manipulating GTF files.

use anyhow::bail;
use flate2::read::MultiGzDecoder;
use noodles::gtf;
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::path::Path;

use super::BioinformaticsFileFormat;

//============================//
// Gene Transfer Format Files //
//============================//

/// Attempts to open a GTF file from a given source.
pub fn open<P>(src: P) -> anyhow::Result<gtf::Reader<Box<dyn BufRead>>>
where
    P: AsRef<Path>,
{
    let path = src.as_ref();
    let file = File::open(path);

    match BioinformaticsFileFormat::try_detect(path) {
        Some(BioinformaticsFileFormat::GTF_GZ) => {
            let reader = file.map(MultiGzDecoder::new).map(BufReader::new)?;
            Ok(gtf::Reader::new(Box::new(reader)))
        }
        Some(BioinformaticsFileFormat::GTF) => {
            let reader = file.map(BufReader::new)?;
            Ok(gtf::Reader::new(Box::new(reader)))
        }
        Some(format) => bail!("incompatible formats: required GTF, found {}", format),
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
