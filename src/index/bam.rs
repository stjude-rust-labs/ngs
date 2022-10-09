//! BAM indexing
//!
//! NOTICE: this was taken almost verbatim from @zaeleus's excellent example in
//! noodles. You can find that source code in the `bam_index.rs` example of the
//! `noodles_bam` crate at commit 0a709087934d.

use anyhow::{bail, Context};
use noodles::bam::{self as bam, bai};
use noodles::csi::index::reference_sequence::bin::Chunk;
use noodles::sam::{self as sam, alignment::Record, header::record::value::map::header::SortOrder};
use num_format::{Locale, ToFormattedString};
use std::{fs::File, io, path::PathBuf};
use tracing::{debug, info};

use crate::utils::formats::sam::parse_header;
use crate::utils::pathbuf::AppendExtension;

//==================================//
// Individual indexing methods: BAM //
//==================================//

/// Checks if a SAM header indicates the file is coordinate sorted.
fn is_coordinate_sorted(header: &sam::Header) -> bool {
    if let Some(hdr) = header.header() {
        if let Some(sort_order) = hdr.sort_order() {
            return sort_order == SortOrder::Coordinate;
        }
    }

    false
}

/// Main method for BAM indexing.
pub fn index(src: PathBuf) -> anyhow::Result<()> {
    info!("indexing BAM file");

    // (1) Reads the file from disk.
    debug!("reading BAM file from disk");
    let mut reader = File::open(&src)
        .map(bam::Reader::new)
        .with_context(|| "opening src file")?;

    // (2) Calculate where the BAM index should go and check if a file is
    // already there. Error out if so.
    let bai = src
        .append_extension("bai")
        .with_context(|| "constructing the BAM index filepath")?;
    if bai.exists() {
        bail!(
            "refusing to overwrite existing index file: {}. Please delete \
                and rerun if you'd like to replace it.",
            bai.display()
        );
    }

    // (3) Read the header and reference sequences.
    let ht = reader.read_header().with_context(|| "reading header")?;
    let header = parse_header(ht);
    reader
        .read_reference_sequences()
        .with_context(|| "reading reference sequences")?;

    // (4) Check if the BAM is coordinate sorted. If it isn't, then return an error.
    debug!("checking the BAM is coordinate sorted");
    if !is_coordinate_sorted(&header) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "the input BAM must be coordinate-sorted to be indexed",
        )
        .into());
    }

    // (5) Build the BAM index.
    debug!("building the BAM index");
    let mut record = Record::default();
    let mut builder = bai::Index::builder();
    let mut start_position = reader.virtual_position();
    let mut records_processed = 0;

    loop {
        match reader.read_record(&mut record) {
            Ok(0) => break,
            Ok(_) => {}
            Err(e) => return Err(e.into()),
        }

        let end_position = reader.virtual_position();
        let chunk = Chunk::new(start_position, end_position);

        builder
            .add_record(&record, chunk)
            .with_context(|| "building BAM index")?;

        start_position = end_position;
        records_processed += 1;

        if records_processed % 1_000_000 == 0 {
            debug!(
                "  [*] Processed {} records.",
                records_processed.to_formatted_string(&Locale::en),
            );
        }
    }

    let index = builder.build(header.reference_sequences().len());

    // (6) Write the index to disk.
    debug!("writing the BAM index to disk");
    let mut writer = File::create(bai)
        .map(bai::Writer::new)
        .with_context(|| "creating BAM index output file")?;

    writer
        .write_header()
        .with_context(|| "writing BAM header")?;
    writer
        .write_index(&index)
        .with_context(|| "writing BAM index")?;

    Ok(())
}
