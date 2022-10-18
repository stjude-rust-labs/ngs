//! Utilities related to opening and manipulating FASTQ files.

use std::fs::File;
use std::io::BufWriter;
use std::io::Write;
use std::path::Path;

use anyhow::bail;
use flate2::write::GzEncoder;
use flate2::Compression;
use noodles::fastq;

use super::BioinformaticsFileFormat;

/// Attempts to open a FASTQ file from a given source.
pub fn writer<P>(src: P) -> anyhow::Result<fastq::Writer<Box<dyn Write>>>
where
    P: AsRef<Path>,
{
    let path = src.as_ref();
    let file = File::create(path);

    match BioinformaticsFileFormat::try_detect(path) {
        Some(BioinformaticsFileFormat::FASTQ_GZ) => {
            let writer = file
                .map(|f| GzEncoder::new(f, Compression::default()))
                .map(BufWriter::new)?;
            Ok(fastq::Writer::new(Box::new(writer)))
        }
        Some(BioinformaticsFileFormat::FASTQ) => {
            let writer = file.map(BufWriter::new)?;
            Ok(fastq::Writer::new(Box::new(writer)))
        }
        Some(format) => bail!("incompatible formats: required FASTQ, found {}", format),
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
