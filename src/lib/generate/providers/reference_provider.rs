//! Functionality related to the `ReferenceGenomeSequenceProvider`. This provider
//! generates records based on a provided reference genome (typically in the
//! FASTA format). In general, the goal is to simulate the sequencing of the
//! provided reference genome.

use std::{
    collections::HashMap,
    fmt::Display,
    path::{Path, PathBuf},
    str::FromStr,
};

use fasta::record::Sequence;
use noodles_core::Position;
use noodles_fasta::{self as fasta};
use noodles_fastq::{self as fastq};
use rand::{rngs::ThreadRng, Rng};
use rand_distr::{Distribution, Normal, WeightedIndex};

use crate::lib::{
    generate::utils::{self, NormalDistributionParams, PairedRead, SeqLen},
    utils::formats,
};

use super::{simulate_errors, SequenceProvider};

//========================================//
// Primary implementation and Constructor //
//========================================//

/// The premier struct that implements the `SequenceProvider` trait. Currently,
/// sequence names and lengths are stored redundantly: once in
/// `sequence_names`/`sequence_lengths` and again in `sequences`. This is to
/// simplify the code when picking a random sequence and was an intentional
/// design decision.
pub struct ReferenceGenomeSequenceProvider {
    /// Filename, just to keep track of which provider we are looking at.
    name: String,
    /// Sequence names stored as a one to one mapping with sequence lengths.
    sequence_names: Vec<String>,
    /// Sequence names stored as a one to one mapping with sequence lengths.
    sequence_lengths: Vec<usize>,
    /// A hashmap containing a mapping between the sequence name and the
    /// sequence itself.
    sequences: HashMap<String, Sequence>,
    /// Distribution to pull sequences from for this reference genome.
    sequence_distribution: WeightedIndex<usize>,

    /// Simulate generating a sequencing error every N nucleotides.
    error_frequency: usize,
    /// Inner distances are pulled from this defined normal distribution.  
    inner_distance_distribution: Normal<f64>,
    /// Assigned read length for the provider.
    read_length: usize,
    /// Total size of the reference genome.
    total_size: usize,

    /// Weight of this reference genome provider. This will be used in
    /// conjunction with the weights of other reference genome providers to
    /// facilitate random sampling.
    weight: usize,
}

impl ReferenceGenomeSequenceProvider {
    /// Tries to create a new `ReferenceGenomeSequenceProvider` from the
    /// provided arguments.
    ///
    /// # Arguments
    ///
    /// * `src` — Path to the reference FASTA for this provider.
    /// * `error_frequency` — Simulate sequencing errors one in every N nucleobases.
    /// * `inner_distance` — Distribution to pull the inner distances from.
    /// * `read_length` — Read length to generate for this provider.
    /// * `weight` — Represents how often this `ReferenceGenomeSequenceProvider`
    ///   should be used to generate reads with respect to the other
    ///   `ReferenceGenomeSequenceProvider`s passed into the command.
    pub fn try_new<P>(
        src: P,
        error_frequency: usize,
        inner_distance: NormalDistributionParams,
        read_length: usize,
        weight: usize,
    ) -> anyhow::Result<Self>
    where
        P: AsRef<Path>,
    {
        // (1) Instantiates the appropriate FASTA reader based on suffix.
        let mut reader = formats::fasta::open(&src)?;
        let filename = src
            .as_ref()
            .file_name()
            .unwrap() // safe, as we expect this must have a filename.
            .to_string_lossy()
            .to_string();

        // (2) Reads the FASTA sequences into memory for fast lookups and
        // calculates some other useful statistics along the way.

        let mut sequence_names = Vec::new();
        let mut sequence_lengths = Vec::new();
        let mut sequences = HashMap::new();
        let mut total_size = 0usize;

        for result in reader.records() {
            let record = result?;
            let name = String::from(record.name());
            let sequence = record.sequence().clone();
            let sequence_length = sequence.len();

            sequence_names.push(name.clone());
            sequence_lengths.push(sequence_length);
            sequences.insert(name, sequence);
            total_size += sequence_length;
        }

        // (3) Generates a sequence distribution from which to pull a random
        // sequence based on the lengths of the sequences read. This is
        // effectively a uniform distribution across the whole genome.
        let sequence_distribution = WeightedIndex::new(&sequence_lengths).unwrap();

        Ok(Self {
            name: filename,
            sequence_names,
            sequence_lengths,
            sequences,
            sequence_distribution,
            error_frequency,
            inner_distance_distribution: Normal::from(inner_distance),
            read_length,
            total_size,
            weight,
        })
    }

    pub fn name(&self) -> &str {
        self.name.as_ref()
    }

    pub fn weight(&self) -> usize {
        self.weight
    }

    /// Selects a random sequence from the reference FASTA. This includes only the
    /// sequence name and the length of the sequence (things such as the exact position
    /// within the sequence are computed later).
    ///
    /// # Arguments
    ///
    /// * `rng` — A `ThreadRng` used to sample the sequence distribution.
    pub fn random_sequence(&self, rng: &mut ThreadRng) -> SeqLen {
        let index = self.sequence_distribution.sample(rng);
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

//===========================//
// Implementation of Display //
//===========================//

impl Display for ReferenceGenomeSequenceProvider {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{{ ReferenceGenomeSequenceProvider: {{ ")?;
        write!(f, "Name: {}, ", self.name)?;
        write!(f, "Sequences: {}, ", self.sequence_names.len())?;
        write!(f, "Total Size: {}, ", self.total_size)?;
        write!(f, "Weight: {} ", self.weight)?;
        write!(f, "}} }}")
    }
}

//===========================//
// Implementation of FromStr //
//===========================//

/// This method is pulled out from the actual `FromStr` impl to make the whole
/// system more easily testable.
fn reference_provider_params_from_str(
    s: &str,
) -> Result<
    (PathBuf, usize, NormalDistributionParams, usize, usize),
    <ReferenceGenomeSequenceProvider as FromStr>::Err,
> {
    let parts: Vec<_> = s.split(':').collect();

    if parts.len() != 6 {
        return Err("Invalid format for reference genome sequence provider, \
            please check the docs."
            .into());
    }

    let src = PathBuf::from(parts[0]);
    let error_frequency = match parts[1].parse::<usize>() {
        Ok(result) => result,
        Err(_) => {
            return Err(format!(
                "Could not parse the error frequency for reference provider: {}.",
                s
            ))
        }
    };

    let mu = match parts[2].parse::<f64>() {
        Ok(result) => result,
        Err(_) => {
            return Err(format!(
            "Could not parse the mean for inner distance distribution for reference provider: {}.",
            s
        ))
        }
    };

    let sigma = match parts[3].parse::<f64>() {
            Ok(result) => result,
            Err(_) => return Err(format!(
                "Could not parse the std deviation for inner distance distribution for reference provider: {}.",
                s
            ))
        };

    let read_length = match parts[4].parse::<usize>() {
        Ok(result) => result,
        Err(_) => {
            return Err(format!(
                "Could not parse the read length for reference provider: {}.",
                s
            ))
        }
    };

    let weight = match parts[5].parse::<usize>() {
        Ok(result) => result,
        Err(_) => {
            return Err(format!(
                "Could not parse the weight for reference provider: {}.",
                s
            ))
        }
    };

    let params = NormalDistributionParams::new(mu, sigma);
    Ok((src, error_frequency, params, read_length, weight))
}

impl FromStr for ReferenceGenomeSequenceProvider {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let (src, error_frequency, params, read_length, weight) =
            reference_provider_params_from_str(s)?;
        match ReferenceGenomeSequenceProvider::try_new(
            src,
            error_frequency,
            params,
            read_length,
            weight,
        ) {
            Ok(rsp) => Ok(rsp),
            Err(e) => Err(e.to_string()),
        }
    }
}

//====================================//
// Implementation of SequenceProvider //
//====================================//

impl SequenceProvider for ReferenceGenomeSequenceProvider {
    /// Generates a random read pair from a given reference genome.
    ///
    /// # Arguments
    ///
    /// * `read_prefix` — The prefixed name to assign to the read (reads are automatically
    ///   interleaved for convenience).
    fn generate_read_pair(&self, read_prefix: String, read_number: String) -> PairedRead {
        let mut rng = ThreadRng::default();
        let mut forward_sequence: Option<Vec<u8>> = None;
        let mut reverse_sequence: Option<Vec<u8>> = None;

        let mut selected_sequence: Option<String> = None;
        let mut forward_start_position: Option<Position> = None;
        let mut reverse_start_position: Option<Position> = None;

        // Generating a read pair can run into a variety of issues. Thus, we loop here
        // until a proper read pair is generated and allow any errors to just reset the
        // loop. This is a very practical strategy to avoiding mysterious errors. Note
        // that the vast majority of calls will see this loop satisfied in one
        // iteration.
        while forward_sequence.is_none() || reverse_sequence.is_none() {
            let random_seq = self.random_sequence(&mut rng);

            let seq = random_seq.get_seq_name();
            let len = random_seq.get_seq_len();

            let min_start = 0usize;
            // at most, read to the end of the chromosome.
            let max_start = len - (self.read_length * 2);
            let start = rng.gen_range(min_start..max_start);
            let mut inner_distance_offset =
                self.inner_distance_distribution.sample(&mut rng).round() as i64;

            // Clamp inner distances so that it can never be less than (mean - 3
            // * std) and never be more than (mean + 3 * std). This helps with
            //   stability and predictability of fragment sizes.

            let lower_bound = (self.inner_distance_distribution.mean()
                - (3.0 * self.inner_distance_distribution.std_dev()).floor())
                as i64;
            let upper_bound = (self.inner_distance_distribution.mean()
                + (3.0 * self.inner_distance_distribution.std_dev()).ceil())
                as i64;

            inner_distance_offset = i64::clamp(inner_distance_offset, lower_bound, upper_bound);

            if let Ok(start_pos) = Position::try_from(start as usize) {
                if let Ok(end_as_isize) = i64::try_from(start + (self.read_length * 2)) {
                    if let Ok(end) = u64::try_from(end_as_isize + inner_distance_offset) {
                        if let Ok(end_pos) = Position::try_from(end as usize) {
                            if let Some(chr) = self.sequences.get(seq) {
                                if let Some(sequence) = chr.get(start_pos..end_pos) {
                                    let sequence_as_vec = sequence.to_vec();
                                    reverse_sequence = utils::reverse_compliment(&sequence_as_vec);
                                    forward_sequence = Some(sequence_as_vec);

                                    selected_sequence = Some(seq.clone());
                                    forward_start_position = Some(start_pos);
                                    reverse_start_position =
                                        Position::new(usize::from(end_pos) - self.read_length);
                                }
                            }
                        }
                    }
                }
            }
        }

        let read_name_one = format!(
            "{}:{}:{}:{}/1",
            &read_prefix,
            selected_sequence.clone().unwrap(),
            forward_start_position.unwrap(),
            &read_number
        );
        let read_name_two = format!(
            "{}:{}:{}:{}/2",
            &read_prefix,
            selected_sequence.unwrap(),
            reverse_start_position.unwrap(),
            &read_number
        );

        let mut fwd_vec = forward_sequence.unwrap();
        let fwd = fwd_vec
            .get(0..self.read_length as usize)
            .unwrap_or_else(|| {
                panic!(
                    "Forward read fragment is too short for the specified read \
             length. This usually means you need to increase the specified \
             inner distance or reduce the standard deviation for genome {} \
             such that fragments this short cannot be generated.",
                    self.name()
                )
            });
        let mut rev_vec = reverse_sequence.unwrap();
        let rev = rev_vec
            .get(0..self.read_length as usize)
            .unwrap_or_else(|| {
                panic!(
                    "Reverse read fragment is too short for the specified read \
             length. This usually means you need to increase the specified \
             inner distance or reduce the standard deviation for genome {} \
             such that fragments this short cannot be generated.",
                    self.name()
                )
            });

        fwd_vec = simulate_errors(fwd, self.error_frequency, &mut rng);
        rev_vec = simulate_errors(rev, self.error_frequency, &mut rng);

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

//=======//
// Tests //
//=======//

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    pub fn test_valid_str() {
        let result = reference_provider_params_from_str("path:20:5.0:6.0:150:100");
        assert!(result.is_ok());

        let (src, error_frequency, params, read_length, weight) = result.unwrap();
        assert_eq!(src, PathBuf::from("path"));
        assert_eq!(error_frequency, 20);
        assert_eq!(params.get_mu(), 5.0);
        assert_eq!(params.get_sigma(), 6.0);
        assert_eq!(read_length, 150);
        assert_eq!(weight, 100);
    }

    #[test]
    pub fn test_invalid_str_malformed() {
        let result = reference_provider_params_from_str("path:20:5.");
        assert!(result.is_err());
        assert_eq!(
            result.err().unwrap(),
            "Invalid format for reference genome sequence provider, please check the docs."
        );
    }

    #[test]
    pub fn test_invalid_str_bad_error_frequency() {
        let result = reference_provider_params_from_str("path:20.0:5.0:6.0:150:100");
        assert!(result.is_err());
        assert_eq!(
            result.err().unwrap(),
            "Could not parse the error frequency for reference provider: path:20.0:5.0:6.0:150:100."
        );
    }

    #[test]
    pub fn test_invalid_str_bad_mean() {
        let result = reference_provider_params_from_str("path:20:err:6.0:150:100");
        assert!(result.is_err());
        assert_eq!(
            result.err().unwrap(),
            "Could not parse the mean for inner distance distribution for reference provider: path:20:err:6.0:150:100."
        );
    }

    #[test]
    pub fn test_invalid_str_bad_std_dev() {
        let result = reference_provider_params_from_str("path:20:5.0:err:150:100");

        assert!(result.is_err());
        assert_eq!(
            result.err().unwrap(),
            "Could not parse the std deviation for inner distance distribution for reference provider: path:20:5.0:err:150:100."
        );
    }

    #[test]
    pub fn test_invalid_str_bad_read_length() {
        let result = reference_provider_params_from_str("path:20:5.0:6.0:false:100");
        assert!(result.is_err());
        assert_eq!(
            result.err().unwrap(),
            "Could not parse the read length for reference provider: path:20:5.0:6.0:false:100."
        );
    }

    #[test]
    pub fn test_invalid_str_bad_weight() {
        let result = reference_provider_params_from_str("path:20:5.0:6.0:150:false");
        assert!(result.is_err());
        assert_eq!(
            result.err().unwrap(),
            "Could not parse the weight for reference provider: path:20:5.0:6.0:150:false."
        );
    }
}
