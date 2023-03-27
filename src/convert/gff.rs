//! Conversions from a GFF file to other next-generation sequencing file formats.

use std::fs::File;
use std::path::PathBuf;

use anyhow::Context;
use noodles::bgzf;
use noodles::bgzf::writer::CompressionLevel;
use noodles::gff;
use noodles::gff::Line::Record;

use crate::utils::args::CompressionStrategy;
use crate::utils::args::NumberOfRecords;
use crate::utils::display::RecordCounter;
use crate::utils::formats;

/// Converts a GFF file to a Block-gzipped GFF file.
pub fn to_block_gzipped_gff(
    from: PathBuf,
    to: PathBuf,
    max_records: NumberOfRecords,
    compression_strategy: CompressionStrategy,
) -> anyhow::Result<()> {
    // (1) Open the GFF file.
    let mut reader = formats::gff::open(from).with_context(|| "opening GFF input file")?;

    // (2) Determine the compression level.
    let compression_level: CompressionLevel = compression_strategy.into();

    // (3) Open the output file.
    let mut writer = File::create(to)
        .map(|f| {
            bgzf::writer::Builder::default()
                .set_compression_level(compression_level)
                .build_with_writer(f)
        })
        .map(gff::Writer::new)
        .with_context(|| "opening bgzipped GFF output file")?;

    // (4) Write the Bgzipped GFF.
    let mut counter = RecordCounter::default();

    for result in reader.lines() {
        let line = result.with_context(|| "reading GFF line")?;
        writer
            .write_line(&line)
            .with_context(|| "writing bgzipped GFF line")?;

        if let Record(_) = line {
            counter.inc();
            if counter.time_to_break(&max_records) {
                break;
            }
        }
    }
    Ok(())
}
