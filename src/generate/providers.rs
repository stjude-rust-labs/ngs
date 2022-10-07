use super::utils::PairedRead;

use rand::{prelude::*, seq::SliceRandom};

pub mod reference_provider;

/// Generic trait that represents anything that can provide a read pair to add
/// to the generated files. The simplest `SequenceProvider` is the implemented
/// `ReferenceGenomeSequenceProvider`, which takes a reference fasta and
/// simulates reads read from that genome. Other sequence providers could be
/// contaminates from various sources.
pub trait SequenceProvider {
    fn generate_read_pair(&self, read_prefix: String, read_number: String) -> PairedRead;
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
