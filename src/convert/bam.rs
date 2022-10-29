//! Conversions from a BAM file to other next-generation sequencing file formats.

use std::path::PathBuf;

use noodles::sam;
use noodles::sam::alignment::Record;
use tokio::fs::File;

use crate::utils::args::NumberOfRecords;
use crate::utils::display::RecordCounter;
use crate::utils::formats;
use crate::utils::formats::bam::ParsedAsyncBAMFile;
use crate::utils::formats::utils::IndexCheck;

/// Converts a BAM file to a SAM file in an asyncronous fashion.
pub async fn to_sam_async(
    from: PathBuf,
    to: PathBuf,
    max_records: NumberOfRecords,
) -> anyhow::Result<()> {
    // (1) Open and parse the BAM file.
    let ParsedAsyncBAMFile {
        mut reader, header, ..
    } = formats::bam::open_and_parse_async(from, IndexCheck::None).await?;

    // (2) Open the SAM file writer.
    let handle = File::create(to).await?;
    let mut writer = sam::AsyncWriter::new(handle);

    // (3) Write the header.
    writer.write_header(&header.parsed).await?;

    let mut counter = RecordCounter::new();
    let mut record = Record::default();

    // (4) Write each record in the BAM file to the SAM file.
    while reader.read_record(&mut record).await? != 0 {
        writer
            .write_alignment_record(&header.parsed, &record)
            .await?;

        counter.inc();

        if counter.time_to_break(&max_records) {
            break;
        }
    }

    Ok(())
}
