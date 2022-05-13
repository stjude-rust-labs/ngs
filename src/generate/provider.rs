use std::{
    fs::File,
    io::{self, BufReader, ErrorKind},
    ops::Range,
    path::PathBuf,
};

use noodles_core::Position;
use noodles_fasta::{self as fasta, repository::adapters::IndexedReader};
use noodles_fastq as fastq;
use rand::{distributions::Uniform, prelude::*};
use rand::{distributions::WeightedIndex, prelude::Distribution};
use tracing::{debug, info};

#[derive(Debug)]
pub struct SeqLen(String, usize);

impl SeqLen {
    pub fn get_seq_name(&self) -> &String {
        &self.0
    }

    pub fn get_seq_len(&self) -> usize {
        self.1
    }
}

#[derive(Debug)]
pub struct PairedRead(fastq::Record, fastq::Record);

impl PairedRead {
    pub fn get_forward_read(&self) -> &fastq::Record {
        &self.0
    }

    pub fn get_reverse_read(&self) -> &fastq::Record {
        &self.1
    }
}

pub trait SequenceProvider {
    fn generate_read_pair(&mut self, read_name: String) -> PairedRead;
}

pub struct ReferenceGenomeSequenceProvider {
    repository: fasta::Repository,
    inner_distance_dist: Uniform<isize>,
    read_length: usize,
    rng: ThreadRng,
    // Key is the sequence name, value is the length of the sequence.
    sequence_names: Vec<String>,
    sequence_lengths: Vec<usize>,
    sequence_dist: WeightedIndex<usize>,
    total_size: usize,
}

impl ReferenceGenomeSequenceProvider {
    pub fn new(
        path: PathBuf,
        read_length: usize,
        inner_distance: Range<isize>,
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
            "fa" => File::open(&path)
                .map(BufReader::new)
                .map(fasta::Reader::new)?,
            ext => {
                return Err(io::Error::new(
                    ErrorKind::InvalidInput,
                    format!("Could not parse reference genome with extension {}.", ext),
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
        let mut total_size: usize = 0;
        let mut sequence_names = Vec::new();
        let mut sequence_lengths = Vec::new();

        for result in reader.records() {
            let record = result?;
            debug!("Reading sequence for {}", record.name());

            let name = record.name().to_owned();
            let len = record.sequence().len();

            sequence_names.push(name);
            sequence_lengths.push(len);
            total_size += len;
        }

        let sequence_dist = WeightedIndex::new(&sequence_lengths).unwrap();

        let repository = fasta::Repository::new(IndexedReader::new(reader, index));
        // (5) Finally, return the initialized object.
        Ok(ReferenceGenomeSequenceProvider {
            repository,
            inner_distance_dist: Uniform::from(inner_distance),
            read_length,
            rng: thread_rng(),
            sequence_names,
            sequence_lengths,
            sequence_dist,
            total_size,
        })
    }

    pub fn precache(&mut self) {
        info!("Precaching sequences...");
        for name in &self.sequence_names {
            debug!("  [*] Caching sequences for {}", name);
            self.repository.get(name);
        }
    }

    pub fn random_sequence(&mut self) -> SeqLen {
        let index = self.sequence_dist.sample(&mut self.rng);
        SeqLen(
            self.sequence_names[index].clone(),
            self.sequence_lengths[index],
        )
    }
}

impl SequenceProvider for ReferenceGenomeSequenceProvider {
    fn generate_read_pair(&mut self, read_name: String) -> PairedRead {
        let mut forward_sequence: Option<Vec<u8>> = None;
        let mut reverse_sequence: Option<Vec<u8>> = None;

        while forward_sequence.is_none() || reverse_sequence.is_none() {
            let random_seq = self.random_sequence();

            let seq = random_seq.get_seq_name();
            let len = random_seq.get_seq_len();

            let min_start = 0usize;
            // at most, read to the end of the chromosome.
            let max_start = len - (self.read_length * 2);
            let start = self.rng.gen_range(min_start..max_start);
            let inner_distance_offset = self.inner_distance_dist.sample(&mut self.rng);

            if let Ok(start_pos) = Position::try_from(start) {
                if let Ok(end_as_isize) = isize::try_from(start + (self.read_length * 2)) {
                    if let Ok(end) = usize::try_from(end_as_isize + inner_distance_offset) {
                        if let Ok(end_pos) = Position::try_from(end) {
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

        let mut read_name_one = read_name.clone();
        let mut read_name_two = read_name;
        read_name_one.push_str("/1");
        read_name_two.push_str("/2");

        PairedRead(
            fastq::Record::new(
                read_name_one,
                forward_sequence.unwrap()[0..self.read_length].to_vec(),
                "J".repeat(self.read_length),
            ),
            fastq::Record::new(
                read_name_two,
                reverse_sequence.unwrap()[0..self.read_length].to_vec(),
                "J".repeat(self.read_length),
            ),
        )
    }
}
