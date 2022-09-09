use std::{
    fs::File,
    io::{self, BufReader, ErrorKind},
    ops::Range,
    path::PathBuf,
};

use super::utils::{PairedRead, SeqLen};

use noodles_core::Position;
use noodles_fasta::{self as fasta, repository::adapters::IndexedReader};
use noodles_fastq as fastq;
use rand::{
    distributions::{Uniform, WeightedIndex},
    prelude::*,
    seq::SliceRandom,
};
use tracing::debug;

/// Generic trait that represents anything that can provide a read pair to add
/// to the generated files. The simplest `SequenceProvider` is the implemented
/// `ReferenceGenomeSequenceProvider`, which takes a reference fasta and
/// simulates reads read from that genome. Other sequence providers could be
/// contaminates from various sources.
pub trait SequenceProvider {
    fn generate_read_pair(&mut self, read_prefix: String) -> PairedRead;
}

/// The premier struct that implements the `SequenceProvider` trait. The
/// `ReferenceGenomeSequenceProvider` takes a reference fasta, as well as other
/// parameters related to sequencing, and generates random read pairs resulting
/// from that genome. This simulates sequencing of a genome in question.
pub struct ReferenceGenomeSequenceProvider {
    error_freq: usize,
    inner_distance_dist: Uniform<i64>,
    repository: fasta::Repository,
    read_length: usize,
    rng: ThreadRng,
    // Key is the sequence name, value is the length of the sequence.
    sequence_names: Vec<String>,
    sequence_lengths: Vec<usize>,
    sequence_dist: WeightedIndex<usize>,
    total_size: usize,
}

impl ReferenceGenomeSequenceProvider {
    /// Creates a new reference genome sequence provider.
    ///
    /// # Arguments
    ///
    /// * `path` — Path to the reference FASTA file.
    /// * `read_length` — Read length to generate for each read.
    /// * `inner_distance` — A negative to positive range that represents the offset
    ///   that the inner distance causes to the read pair. Negative offsets represent
    ///   the read pairs overlapping and positive offsets represent space between the
    ///   reads in the read pair.
    /// * `error_freq` — The provider will simulate a sequencing error for one in every
    ///   N nucleotides sequenced.
    pub fn new(
        path: PathBuf,
        read_length: usize,
        inner_distance: Range<i64>,
        error_freq: usize,
    ) -> io::Result<Self> {
        // (1) Parses the extension for the reference file.
        let extension = path
            .extension()
            .expect("Could not parse extension from reference file.")
            .to_str()
            .expect("Could not convert extension to string.");

        // (2) Instantiates a reader that is appropriate based on the reference
        // file's type. For the moment, we only parse uncompressed FASTA files,
        // but we may expand this in the future (thus, the match statement).
        let mut reader = match extension {
            "fa" | "fna" | "fasta" => File::open(&path)
                .map(BufReader::new)
                .map(fasta::Reader::new)?,
            ext => {
                return Err(io::Error::new(
                    ErrorKind::InvalidInput,
                    format!("Could not parse reference genome with extension {}. Supported extension types are: fa, fna, fasta.", ext),
                ))
            }
        };

        // (3) Identify and load in the associated index file for our reference
        // file.
        let fai_ext = format!("{}.fai", extension);
        let index = match fasta::fai::read(path.with_extension(fai_ext)) {
            Ok(i) => Ok(i),
            _ => Err(io::Error::new(
                ErrorKind::InvalidInput,
                "Could not find or read associated \
            index file for the reference file. Please use `samtools faidx` on \
            the reference fasta to generate one.",
            )),
        }?;

        // (4) Load all of the sequences and their lengths into a cache for
        // later use.
        let mut sequence_names = Vec::new();

        for result in reader.records() {
            let record = result?;
            let name = record.name().to_owned();
            sequence_names.push(name);
        }

        let repository = fasta::Repository::new(IndexedReader::new(reader, index));

        let mut total_size: usize = 0;
        let mut sequence_lengths = Vec::new();
        for seq in &sequence_names {
            // precaching and length lookup
            debug!("Reading sequence for {}", &seq);
            let len = repository.get(seq).unwrap().unwrap().len();
            sequence_lengths.push(len);
            total_size += len;
        }

        let sequence_dist = WeightedIndex::new(&sequence_lengths).unwrap();

        // (5) Finally, return the initialized object.
        Ok(ReferenceGenomeSequenceProvider {
            error_freq,
            inner_distance_dist: Uniform::from(inner_distance),
            read_length,
            repository,
            rng: thread_rng(),
            sequence_names,
            sequence_lengths,
            sequence_dist,
            total_size,
        })
    }

    /// Selects a random sequence from the reference FASTA. This includes only the
    /// sequence name and the length of the sequence (things such as the exact position
    /// within the sequence are computed later).
    pub fn random_sequence(&mut self) -> SeqLen {
        let index = self.sequence_dist.sample(&mut self.rng);
        SeqLen(
            self.sequence_names[index].clone(),
            self.sequence_lengths[index],
        )
    }

    /// Calculates the number of reads needed for a mean coverage of `coverage`. The
    /// calculation is simple and somewhat naive. First, we determine how many reads
    /// would be needed end to end to covere the whole genome one time. Then, we
    /// multiply that by the coverage we'd like. Of course, this does not guarentee that
    /// you will see `coverage`x coverage at any given site of interest.
    ///
    /// # Arguments
    ///
    /// * `coverage` — The desired mean coverage for this generated file.
    pub fn reads_needed_for_coverage(&self, coverage: usize) -> usize {
        coverage * (self.total_size / self.read_length)
    }
}

impl SequenceProvider for ReferenceGenomeSequenceProvider {
    /// Generates a random read pair from a given reference genome.
    ///
    /// # Arguments
    ///
    /// * `read_prefix` — The prefixed name to assign to the read (reads are automatically
    ///   interleaved for convenience).
    fn generate_read_pair(&mut self, read_prefix: String) -> PairedRead {
        let mut forward_sequence: Option<Vec<u8>> = None;
        let mut reverse_sequence: Option<Vec<u8>> = None;

        // Generating a read pair can run into a variety of issues. Thus, we loop here
        // until a proper read pair is generated and allow any errors to just reset the
        // loop. This is a very practical strategy to avoiding mysterious errors. Note
        // that the vast majority of calls will see this loop satisfied in one
        // iteration.
        while forward_sequence.is_none() || reverse_sequence.is_none() {
            let random_seq = self.random_sequence();

            let seq = random_seq.get_seq_name();
            let len = random_seq.get_seq_len();

            let min_start = 0usize;
            // at most, read to the end of the chromosome.
            let max_start = len - (self.read_length * 2);
            let start = self.rng.gen_range(min_start..max_start);
            let inner_distance_offset = self.inner_distance_dist.sample(&mut self.rng);

            if let Ok(start_pos) = Position::try_from(start as usize) {
                if let Ok(end_as_isize) = i64::try_from(start + (self.read_length * 2)) {
                    if let Ok(end) = u64::try_from(end_as_isize + inner_distance_offset) {
                        if let Ok(end_pos) = Position::try_from(end as usize) {
                            if let Some(Ok(chr)) = self.repository.get(seq) {
                                if let Some(sequence) = chr.get(start_pos..end_pos) {
                                    let sequence_as_vec = sequence.to_vec();
                                    reverse_sequence =
                                        super::utils::reverse_compliment(&sequence_as_vec);
                                    forward_sequence = Some(sequence_as_vec);
                                }
                            }
                        }
                    }
                }
            }
        }

        let mut read_name_one = read_prefix.clone();
        let mut read_name_two = read_prefix;
        read_name_one.push_str("/1");
        read_name_two.push_str("/2");

        let mut fwd_vec = forward_sequence.unwrap();
        let fwd = fwd_vec.get(0..self.read_length as usize).unwrap();
        let mut rev_vec = reverse_sequence.unwrap();
        let rev = rev_vec.get(0..self.read_length as usize).unwrap();

        fwd_vec = simulate_errors(fwd, self.error_freq, &mut self.rng);
        rev_vec = simulate_errors(rev, self.error_freq, &mut self.rng);

        PairedRead(
            fastq::Record::new(
                read_name_one,
                fwd_vec,
                "J".repeat(self.read_length as usize),
            ),
            fastq::Record::new(
                read_name_two,
                rev_vec,
                "J".repeat(self.read_length as usize),
            ),
        )
    }
}

/// Simulates sequencing errors in a very naive way.
///
/// # Arguments
///
/// * `seq` — The byte string containing the existing sequence.
/// * `error_freq` — The frequency with which to simulate sequencing errors.
/// * `rng` — The random number generator to use.
///
/// TODO: this function could be greatly improved!
pub fn simulate_errors(seq: &[u8], error_freq: usize, rng: &mut ThreadRng) -> Vec<u8> {
    // Represents A, C, G, T in byte form.
    let bases: Vec<u8> = vec![0x41, 0x43, 0x47, 0x54];

    // Loops through the sequence, if our random number generator hits, we keep spinning
    // the slot machine until we find a base not equal to the true base in the genome.
    let mut result = Vec::new();
    for mut current in seq.iter().cloned() {
        if rng.gen_ratio(1, error_freq as u32) {
            let mut new = current;
            while current == new {
                new = *bases.choose(rng).unwrap();
            }
            current = new;
        }
        result.push(current);
    }

    result
}
