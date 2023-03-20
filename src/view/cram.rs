//! CRAM viewing

use std::fs::File;
use std::io;
use std::io::Write;
use std::path::PathBuf;

use anyhow::bail;
use anyhow::Context;
use noodles::cram;
use noodles::cram::crai;
use noodles::fasta;
use noodles::fasta::repository::adapters::IndexedReader;
use noodles::sam;
use noodles::sam::AlignmentWriter;
use tracing::debug;

use crate::utils::formats::sam::parse_header;
use crate::utils::pathbuf::AppendExtension;
use crate::view::command::Mode;

/// Main method for BAM viewing.
pub fn view(
    src: PathBuf,
    query: Option<String>,
    reference_fasta: PathBuf,
    mode: Mode,
) -> anyhow::Result<()> {
    // (1) Reads the file from disk.
    debug!("reading CRAM file from disk");
    let mut reader = File::open(&src)
        .map(cram::Reader::new)
        .with_context(|| "opening src file")?;

    // (2) Determine the handle with which to write the output.
    let stdout = io::stdout();
    let mut handle = stdout.lock();

    // (3) Build the FASTA repository and associated index.
    let repository = fasta::indexed_reader::Builder::default()
        .build_from_path(&reference_fasta)
        .map(IndexedReader::new)
        .map(fasta::Repository::new)
        .with_context(|| "building FASTA repository and associated index")?;

    // TODO: remove in a future version when noodles gives an error message that
    // suggests you should index your FASTA file (as of the time of writing, it
    // just gives an error back saying "unsupported" if you don't have an
    // index).
    let fai_filepath = reference_fasta.append_extension("fai")?;
    if !fai_filepath.exists() {
        bail!(
            "couldn't find an index for your reference FASTA: is the FASTA indexed? \
        Run `ngs index [FASTA]` to index the FASTA file."
        )
    }

    // (4) Read the file's definition.
    reader.read_file_definition()?;

    // (5) If the user specified to output the header, output the raw header (before
    // applying any corrections).
    let ht = reader
        .read_file_header()
        .with_context(|| "reading CRAM header")?;

    if mode == Mode::Full || mode == Mode::HeaderOnly {
        write!(handle, "{}", ht).with_context(|| "writing header to stream")?;
    }

    // (6) If the mode is header-only, nothing left to do, so return.
    if mode == Mode::HeaderOnly {
        return Ok(());
    }

    // (7) Parse the header and apply corrections.
    let header = parse_header(ht).with_context(|| "parsing header")?;

    // (8) Writes the records to the output stream.
    let mut writer = sam::Writer::new(handle);

    if let Some(query) = query {
        // (a) If a query is specified, print just the records that fall within the query.
        let index =
            crai::read(src.with_extension("cram.crai")).with_context(|| "reading CRAM index")?;
        let region = query.parse().with_context(|| "parsing query")?;

        let records = reader
            .query(&repository, &header, &index, &region)
            .with_context(|| "querying CRAM file")?;

        for result in records {
            let record = result
                .and_then(|record| record.try_into_alignment_record(&header))
                .with_context(|| "reading record")?;
            writer.write_alignment_record(&header, &record)?;
        }
    } else {
        // (b) Else, print all of the records in the file.
        for result in reader.records(&repository, &header) {
            let record = result
                .and_then(|record| record.try_into_alignment_record(&header))
                .with_context(|| "reading record")?;
            writer.write_alignment_record(&header, &record)?;
        }
    }

    writer.finish(&header)?;
    Ok(())
}
