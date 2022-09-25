use std::{fs::File, io::BufReader, path::PathBuf};

use anyhow::{bail, Context};
use fasta::record::Sequence;
use noodles_fasta as fasta;
use noodles_sam::Header;
use serde::Serialize;

use crate::lib::utils::{
    alignment::ReferenceRecordStrepthrough, formats, histogram::SimpleHistogram,
};

use super::SequenceBasedQualityCheckFacet;

#[derive(Clone, Debug, Default, Serialize)]
pub struct EditMetricsSummary {
    pub mean_edits_read_one: f64,
    pub mean_edits_read_two: f64,
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct EditMetrics {
    pub read_one_edits: SimpleHistogram,
    pub read_two_edits: SimpleHistogram,
    pub summary: Option<EditMetricsSummary>,
}

pub struct EditsFacet<'a> {
    pub metrics: EditMetrics,
    pub fasta: fasta::Reader<BufReader<File>>,
    pub current_sequence: Option<Sequence>,
    pub header: &'a Header,
}

impl<'a> EditsFacet<'a> {
    pub fn try_from(reference_fasta: &PathBuf, header: &'a Header) -> anyhow::Result<Self> {
        let fasta = formats::fasta::open(reference_fasta).with_context(|| {
            format!(
                "Error opening reference FASTA file: {}.",
                reference_fasta.display()
            )
        })?;

        Ok(EditsFacet {
            metrics: EditMetrics::default(),
            fasta,
            current_sequence: None,
            header,
        })
    }
}

impl<'a> SequenceBasedQualityCheckFacet<'a> for EditsFacet<'a> {
    fn name(&self) -> &'static str {
        "Edit Metrics"
    }

    fn computational_load(&self) -> super::ComputationalLoad {
        super::ComputationalLoad::Heavy
    }

    fn supports_sequence_name(&self, _: &str) -> bool {
        true
    }

    fn setup_sequence(
        &mut self,
        seq: &noodles_sam::header::ReferenceSequence,
    ) -> anyhow::Result<()> {
        let seq_name = seq.name().as_str();

        for result in self.fasta.records() {
            let record = result?;
            if seq_name == record.name() {
                self.current_sequence = Some(record.sequence().clone());
                return Ok(());
            }
        }

        bail!("Sequence {} not found in reference FASTA.", seq_name)
    }

    fn process_record<'b>(
        &mut self,
        _: &'b noodles_sam::header::ReferenceSequence,
        record: &noodles_sam::alignment::Record,
    ) -> anyhow::Result<()>
    where
        'b: 'a,
    {
        // (1) First, if the read is unmapped, we need to ignore it for this
        // analysis because there is no reference to compare it to.
        if record.flags().is_unmapped() {
            return Ok(());
        }

        // (2) Parse the read name.
        let read_name = match record.read_name() {
            Some(name) => name,
            _ => bail!("Could not parse read name"),
        };

        // (3) Parse the reference sequence id from the record
        // let id = match record.reference_sequence_id() {
        //     Some(id) => id,
        //     _ => {
        //         return Err(Error::new(format!(
        //             "Could not parse reference sequence id for read: {}",
        //             read_name
        //         )))
        //     }
        // };

        // (4) Lookup this record's reference sequence name
        // let seq_name = match self
        //     .header
        //     .reference_sequences()
        //     .get_index(id)
        //     .map(|(_, rs)| Some(rs.name().as_str()))
        // {
        //     Some(Some(name)) => name,
        //     _ => {
        //         return Err(Error::new(format!(
        //             "Could not map reference sequence id to header for read: {}",
        //             read_name
        //         )))
        //     }
        // };

        // (4e) Do the actual calculation
        let reference_start = record.alignment_start().unwrap();
        let cigar = record.cigar();
        let alignment_span = cigar.alignment_span();
        let reference_end = reference_start.checked_add(alignment_span).unwrap();

        if let Some(current_sequence) = &self.current_sequence {
            let reference_seq = match current_sequence.get(reference_start..reference_end) {
                Some(result) => result,
                None => bail!(
                    "Could not lookup reference sequence for read: {}",
                    read_name
                ),
            };
            let record_seq_sequence = record.sequence();
            let record_seq = record_seq_sequence.as_ref();

            let rrs = ReferenceRecordStrepthrough::new(reference_seq, record_seq, cigar.clone());
            let edits = rrs.edits();

            if record.flags().is_first_segment() {
                self.metrics.read_one_edits.increment(edits).unwrap();
            } else {
                self.metrics.read_two_edits.increment(edits).unwrap();
            }
        }

        Ok(())
    }

    fn teardown_sequence(
        &mut self,
        _: &noodles_sam::header::ReferenceSequence,
    ) -> anyhow::Result<()> {
        self.current_sequence = None;
        Ok(())
    }

    fn aggregate_results(&mut self, results: &mut super::results::Results) {
        self.metrics.summary = Some(EditMetricsSummary {
            mean_edits_read_one: self.metrics.read_one_edits.mean(),
            mean_edits_read_two: self.metrics.read_two_edits.mean(),
        });

        results.set_edits(self.metrics.clone());
    }
}