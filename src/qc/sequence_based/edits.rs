//! Functionality related to the Edits quality control facet.

use std::fs::File;
use std::io::Write;
use std::path::PathBuf;

use anyhow::bail;
use anyhow::Context;
use fasta::record::Sequence;
use noodles::fasta;
use noodles::sam::alignment::Record;
use noodles::sam::header::record::value::{map::ReferenceSequence, Map};
use noodles::sam::record::cigar::op::Kind;
use noodles::sam::record::sequence::base::TryFromCharError;
use noodles::sam::record::sequence::Base;
use serde::Deserialize;
use serde::Serialize;

use crate::qc::results;
use crate::qc::ComputationalLoad;
use crate::qc::SequenceBasedQualityControlFacet;
use crate::utils::alignment::ReferenceRecordStepThrough;
use crate::utils::formats;
use crate::utils::histogram::Histogram;

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
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct EditMetrics {
    /// The distribution of edit counts for all read ones in the file.
    pub read_one_edits: Histogram,

    /// The distribution of edit counts for all read twos in the file.
    pub read_two_edits: Histogram,

    /// Histogram representing the distribution of variant allele fractions (VAF) across
    /// all variants observed.
    pub vaf_histogram: Histogram,

    /// Summary statistics for the Edits quality control facet.
    pub summary: Option<EditMetricsSummary>,
}

impl Default for EditMetrics {
    fn default() -> Self {
        Self {
            read_one_edits: Histogram::default(),
            read_two_edits: Histogram::default(),
            vaf_histogram: Histogram::zero_based_with_capacity(100),
            summary: None,
        }
    }
}

/// Primary struct related to the Edits quality control facet.
pub struct EditsFacet {
    /// Metrics related to the Edits quality control facet.
    ///
    /// This is what will ultimately store the output of this quality control facet.
    pub metrics: EditMetrics,

    /// Refs per position within the current sequence being evaluated.
    ///
    /// This is structured as a Histogram that covers all positions in the current
    /// sequence being evaluated. The value in each bin represents how many records
    /// supporting the reference allele at that position were detected.
    pub refs_per_position: Histogram,

    /// Alts per position within the current sequence being evaluated.
    ///
    /// This is structured as a Histogram that covers all positions in the current
    /// sequence being evaluated. The value in each bin represents how many records
    /// supporting an alternative allele at that position were detected. No distinction
    /// is made for differing alternative alleles (bi-allelic sites, etc).
    pub alts_per_position: Histogram,

    /// The FASTA file path.
    ///
    /// A path to the FASTA  is retained so that sequences can be retrieved on
    /// demand. This is essentially a space optimization: we only go to the file and
    /// pull the sequence into memory once it is that sequence's turn to be processed.
    pub reference_fasta_file_path: PathBuf,

    /// The sequence currently being processed by the quality control facet.
    ///
    /// This is the other part of the space optimization listed in the `fasta`
    /// description. This field is updated over time as the
    /// [`setup`](../../trait.SequenceBasedQualityControlFacet.html#tymethod.setup)
    /// method is called on each sequence to ensure only one sequence is in memory at a
    /// time.
    pub current_sequence: Option<Sequence>,

    /// File which contains a variant allele fraction (VAF) for every position in the
    /// genome.
    ///
    /// If `None` is provided, then no VAF file will be output. If `Some(handle)` is
    /// provided, then the VAF contents will be written to `handle`.
    pub vaf_file: Option<File>,
}

impl EditsFacet {
    /// Tries to create an [`EditsFacet`] from a reference FASTA file
    /// (`reference_fasta`). If a `vaf_file_path` is provided, then the Edits facet will
    /// also output a VAF file for every position covered by a record.
    pub fn try_from(
        reference_fasta: &PathBuf,
        vaf_file_path: Option<PathBuf>,
    ) -> anyhow::Result<Self> {
        // (1) Tries to open the reference FASTA file. We do this once without storing
        // the result to make sure that all is well with the reference FASTA file before
        // proceeding.
        formats::fasta::open(reference_fasta).with_context(|| {
            format!(
                "opening reference FASTA file: {}.",
                reference_fasta.display()
            )
        })?;

        // (2) Tries to prepare the `vaf_file_path`, if it exists.
        let vaf_file = match vaf_file_path {
            Some(file_path) => {
                if file_path.exists() {
                    bail!(
                        "refusing to overwrite existing VAF file: {}. \
                        Please delete and rerun if you'd like to replace it.",
                        file_path.display()
                    )
                }

                let mut f = File::create(file_path).with_context(|| "creating VAF file")?;
                writeln!(f, "Sequence\tPosition\tVAF")
                    .with_context(|| "writing VAF file header")?;
                Some(f)
            }
            None => None,
        };

        Ok(EditsFacet {
            metrics: EditMetrics::default(),
            refs_per_position: Histogram::default(),
            alts_per_position: Histogram::default(),
            reference_fasta_file_path: reference_fasta.clone(),
            current_sequence: None,
            vaf_file,
        })
    }
}

impl SequenceBasedQualityControlFacet for EditsFacet {
    fn name(&self) -> &'static str {
        "Edits"
    }

    fn computational_load(&self) -> ComputationalLoad {
        ComputationalLoad::Heavy
    }

    fn supports_sequence_name(&self, _: &str) -> bool {
        true
    }

    fn setup(&mut self, sequence: &Map<ReferenceSequence>) -> anyhow::Result<()> {
        let seq_name = sequence.name().as_str();

        // (1) Load the reference sequence we're currently evaluating into memory.
        let mut fasta =
            formats::fasta::open(&self.reference_fasta_file_path).with_context(|| {
                format!(
                    "opening reference FASTA file: {}.",
                    self.reference_fasta_file_path.display()
                )
            })?;

        // (2) Hook the reference genome sequence that we are currently process it and
        // store it for later use.
        for result in fasta.records() {
            let record = result?;
            if seq_name == record.name() {
                self.current_sequence = Some(record.sequence().clone());
                break;
            }
        }

        if self.current_sequence.is_none() {
            bail!("sequence {} not found in reference FASTA.", seq_name)
        }

        // (3) Set up the `refs_per_position` and `alts_per_position` histograms for
        // tracking during processing.

        let seq_length = usize::from(sequence.length());
        self.refs_per_position = Histogram::zero_based_with_capacity(seq_length);
        self.alts_per_position = Histogram::zero_based_with_capacity(seq_length);

        Ok(())
    }

    fn process<'b>(
        &mut self,
        _: &'b Map<ReferenceSequence>,
        record: &Record,
    ) -> anyhow::Result<()> {
        // (1) First, if the read is unmapped, we need to ignore it for this
        // analysis because there is no reference to compare it to. Further, we don't
        // care about duplicate reads for this analysis (and they will only serve to
        // mess up our calculations).
        if record.flags().is_unmapped() || record.flags().is_duplicate() {
            return Ok(());
        }

        // (2) Parse the record name. This is useful for reporting an error with
        // the record if appropriate.
        let read_name = match record.read_name() {
            Some(name) => name,
            _ => bail!("Could not parse read name"),
        };

        // (3) Step-through each record in the file, counting up the number of edits as
        // well as the refs and alts at each position (to eventually calculate the VAFs
        // in the teardown step).
        let cigar = record.cigar();
        let reference_start = record.alignment_start().unwrap();
        let reference_end = reference_start.checked_add(cigar.alignment_span()).unwrap();

        let current_sequence = match &self.current_sequence {
            Some(s) => s,
            None => bail!(
                "could not lookup reference sequence for read: {}",
                read_name
            ),
        };

        // Map the reference sequences retrieved to `Base`s so they can be easily
        // compared. TODO: there may be a better way to do this without so many
        // conversions and copies, which would, in turn, speed up the facet and require
        // less memory.
        let reference_seq_vec: Result<Vec<Base>, TryFromCharError> = current_sequence
            .get(reference_start..reference_end)
            .map(|x| x.iter().copied().map(Base::try_from))
            .unwrap()
            .collect();

        let reference_seq = reference_seq_vec?;
        let reference_seq = reference_seq.as_ref();
        let record_seq = record.sequence().as_ref();

        // Starts the reference-record step through.
        let rrs = ReferenceRecordStepThrough::new(reference_seq, record_seq, cigar.clone());

        // TODO: this definition of edits, whereby an edit is only counted as one when
        // (a) the CIGAR string says it's an `M` and (b) the bases don't match could
        // probably be improved to include the edit distance of deletes, insertions,
        // etc.
        let mut edits = 0;
        rrs.stepthrough(|cigar, reference_base, reference_ptr, record_base, _| {
            if cigar == Kind::Match {
                let reference_position = usize::from(reference_start)
                    .checked_add(reference_ptr)
                    .unwrap();

                if reference_base != record_base {
                    edits += 1;
                    self.alts_per_position
                        .increment(reference_position)
                        .unwrap();
                } else {
                    self.refs_per_position
                        .increment(reference_position)
                        .unwrap();
                }
            }

            Ok(())
        })?;

        if record.flags().is_first_segment() {
            self.metrics.read_one_edits.increment(edits).unwrap();
        } else {
            self.metrics.read_two_edits.increment(edits).unwrap();
        }

        Ok(())
    }

    fn teardown(&mut self, seq: &Map<ReferenceSequence>) -> anyhow::Result<()> {
        let seq_name = seq.name().to_string();

        // (1) Resets the current sequence to None so that there can be no
        // confusion as to whether the sequence we were looking at is still
        // around in the next iteration.
        self.current_sequence = None;

        // (3) Calculate up the VAFs at every position for this sequence. If the VAF
        // file option was provided (and thus the `vaf_file` key is `Some`), then write
        // every VAF as a tab-delimited file at that location.
        for i in self.refs_per_position.range_start()..=self.refs_per_position.range_stop() {
            let refs_at_this_position = self.refs_per_position.get(i);
            let alts_at_this_position = self.alts_per_position.get(i);
            let total_at_this_position = refs_at_this_position + alts_at_this_position;

            if total_at_this_position == 0 {
                // No records covered this position, simply continue to the next
                // position.
                continue;
            }

            let vaf_at_this_position = alts_at_this_position as f32 / total_at_this_position as f32;
            self.metrics
                .vaf_histogram
                .increment((vaf_at_this_position * 100.0) as usize)
                .unwrap();

            if let Some(f) = &mut self.vaf_file {
                writeln!(f, "{}\t{}\t{}", seq_name, i, vaf_at_this_position)
                    .with_context(|| "writing VAF file")?;
            }
        }

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
