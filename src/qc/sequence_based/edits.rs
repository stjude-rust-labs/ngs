//! Functionality related to the Edits quality control facet.

use std::{fs::File, io::BufReader, path::PathBuf};

use anyhow::{bail, Context};
use fasta::record::Sequence;
use noodles::fasta;
use noodles::sam::{
    alignment::Record,
    header::record::value::{map::ReferenceSequence, Map},
};
use serde::{Deserialize, Serialize};

use crate::{
    qc::{results, ComputationalLoad, SequenceBasedQualityControlFacet},
    utils::{alignment::ReferenceRecordStepThrough, formats, histogram::Histogram},
};

//=========//
// Metrics //
//=========//

/// Summary statistics for the Edits quality control facet.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct EditMetricsSummary {
    /// Mean number of edits for the read ones in the file.
    pub mean_edits_read_one: f64,

    /// Mean number of edits for the read twos in the file.
    pub mean_edits_read_two: f64,
}

/// Primary metrics struct that is comprised of all of the minor metrics structs
/// for this quality control facet.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct EditMetrics {
    /// The distribution of edit counts for all read ones in the file.
    pub read_one_edits: Histogram,

    /// The distribution of edit counts for all read twos in the file.
    pub read_two_edits: Histogram,

    /// Summary statistics for the Edits quality control facet.
    pub summary: Option<EditMetricsSummary>,
}

/// Primary struct used to compile stats regarding edits.
pub struct EditsFacet {
    /// Metrics related to the Edits quality control facet.
    pub metrics: EditMetrics,

    /// The FASTA reader, which is used to cache the current sequence being
    /// reviewed as processing occurs (recall that this is a sequence-based
    /// quailty check facet, so all of the sequences cannot be held in memory at
    /// the same time).
    pub fasta: fasta::Reader<BufReader<File>>,

    /// The sequence currently being processed by the quality control facet.
    /// This is updated over time as the
    /// [`setup`](../../trait.SequenceBasedQualityControlFacet.html#tymethod.setup)
    /// is called.
    pub current_sequence: Option<Sequence>,
}

impl EditsFacet {
    /// Tries to create an [`EditsFacet`] from a reference FASTA file.
    pub fn try_from(reference_fasta: PathBuf) -> anyhow::Result<Self> {
        let fasta = formats::fasta::open(&reference_fasta).with_context(|| {
            format!(
                "Error opening reference FASTA file: {}.",
                reference_fasta.display()
            )
        })?;

        Ok(EditsFacet {
            metrics: EditMetrics::default(),
            fasta,
            current_sequence: None,
        })
    }
}

impl SequenceBasedQualityControlFacet for EditsFacet {
    fn name(&self) -> &'static str {
        "Edit Metrics"
    }

    fn computational_load(&self) -> ComputationalLoad {
        ComputationalLoad::Heavy
    }

    fn supports_sequence_name(&self, _: &str) -> bool {
        true
    }

    fn setup(&mut self, sequence: &Map<ReferenceSequence>) -> anyhow::Result<()> {
        let seq_name = sequence.name().as_str();

        for result in self.fasta.records() {
            let record = result?;
            if seq_name == record.name() {
                self.current_sequence = Some(record.sequence().clone());
                return Ok(());
            }
        }

        bail!("Sequence {} not found in reference FASTA.", seq_name)
    }

    fn process<'b>(
        &mut self,
        _: &'b Map<ReferenceSequence>,
        record: &Record,
    ) -> anyhow::Result<()> {
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

        // (3) Do the actual calculation
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

            let rrs = ReferenceRecordStepThrough::new(reference_seq, record_seq, cigar.clone());
            let edits = rrs.edits();

            if record.flags().is_first_segment() {
                self.metrics.read_one_edits.increment(edits).unwrap();
            } else {
                self.metrics.read_two_edits.increment(edits).unwrap();
            }
        }

        Ok(())
    }

    fn teardown(&mut self, _: &Map<ReferenceSequence>) -> anyhow::Result<()> {
        self.current_sequence = None;
        Ok(())
    }

    fn aggregate(&mut self, results: &mut results::Results) {
        self.metrics.summary = Some(EditMetricsSummary {
            mean_edits_read_one: self.metrics.read_one_edits.mean(),
            mean_edits_read_two: self.metrics.read_two_edits.mean(),
        });

        results.edits = Some(self.metrics.clone());
    }
}
