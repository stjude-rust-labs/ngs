//! Utilities related to opening and manipulating FASTQ files.

use std::{
    fs::File,
    io::{BufWriter, Write},
    path::Path,
};

use anyhow::bail;
use flate2::{write::GzEncoder, Compression};
use noodles::fastq;

/// Attempts to open a FASTQ file from a given source.
pub fn writer<P>(src: P) -> anyhow::Result<fastq::Writer<Box<dyn Write>>>
where
    P: AsRef<Path>,
{
    let path = src.as_ref();
    let file = File::create(path);
    match path.extension().and_then(|x| x.to_str()) {
        Some("gz") => {
            let writer = file
                .map(|f| GzEncoder::new(f, Compression::default()))
                .map(BufWriter::new)?;
            Ok(fastq::Writer::new(Box::new(writer)))
        }
        Some("fq") | Some("fastq") => {
            let writer = file.map(BufWriter::new)?;
            Ok(fastq::Writer::new(Box::new(writer)))
        }
        _ => bail!("Unknown extension for FASTQ file."),
    }
}
