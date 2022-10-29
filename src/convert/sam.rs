//! Conversions from a SAM file to other next-generation sequencing file formats.

use std::path::PathBuf;

use noodles::bam;
use noodles::sam::alignment::Record;
use tokio::fs::File;

use crate::utils::args::NumberOfRecords;
use crate::utils::display::RecordCounter;
use crate::utils::formats;
use crate::utils::formats::sam::ParsedAsyncSAMFile;

/// Converts a SAM file to a BAM file in an asyncronous fashion.
pub async fn to_bam_async(
    from: PathBuf,
    to: PathBuf,
    max_records: NumberOfRecords,
) -> anyhow::Result<()> {
    // (1) Open and parse the BAM file.
    let ParsedAsyncSAMFile { mut reader, header } =
        formats::sam::open_and_parse_async(from).await?;

    // (2) Open the SAM file writer.
    let handle = File::create(to).await?;
    let mut writer = bam::AsyncWriter::new(handle);

    // (3) Write the header and the reference sequences.
    writer.write_header(&header.parsed).await?;
    writer
        .write_reference_sequences(header.parsed.reference_sequences())
        .await?;

    let mut counter = RecordCounter::new();
    let mut record = Record::default();

    // (4) Write each record in the BAM file to the SAM file.
    while reader.read_record(&header.parsed, &mut record).await? != 0 {
        writer
            .write_alignment_record(&header.parsed, &record)
            .await?;

        counter.inc();

        if counter.time_to_break(&max_records) {
            break;
        }
    }

    // (5) Shutdown the async writer.
    writer.shutdown().await?;

    Ok(())
}
