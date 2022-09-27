pub mod ncbi;
pub mod one_thousand_genomes;

use std::{fmt::Debug, rc::Rc, str::FromStr};

use self::ncbi::grch38_no_alt::GRCh38NoAlt;
use self::one_thousand_genomes::hs37d5::HS37D5;

//=================//
// Utility methods //
//=================//

/// Gets all of the supported reference genomes for the tool. When new reference
/// genomes are added, this needs to be updated.
pub fn get_all_reference_genomes() -> Vec<Box<dyn ReferenceGenome>> {
    vec![Box::new(HS37D5), Box::new(GRCh38NoAlt)]
}

/// Utility method to map a string (generally passed on the command line) to a
/// reference genome. I considered implementing FromStr for ReferenceGenome, but
/// that solution felt a lot less clean considering we have the names stored in
/// `.name()` and want to do the lowercase matching.
pub fn get_reference_genome(s: &str) -> Option<Box<dyn ReferenceGenome>> {
    get_all_reference_genomes()
        .into_iter()
        .find(|genome| s.eq_ignore_ascii_case(genome.name()))
}

/// Gets the primary assembly for a given reference genome (as defined by the
/// autosomes, the sex chromosomes, the unlocalized sequences, and the unplaced
/// sequences).
pub fn get_primary_assembly(reference_genome: Rc<Box<dyn ReferenceGenome>>) -> Vec<Sequence> {
    let mut primary: Vec<Sequence> = Vec::new();

    if let Some(autosomes) = reference_genome.autosomes() {
        primary.extend(autosomes);
    }

    if let Some(sex_chromosomes) = reference_genome.sex_chromosomes() {
        primary.extend(sex_chromosomes);
    }

    if let Some(unlocalized) = reference_genome.unlocalized_sequences() {
        primary.extend(unlocalized);
    }

    if let Some(unplaced) = reference_genome.unplaced_sequences() {
        primary.extend(unplaced);
    }

    primary
}

/// Gets all sequences for the given assembly.
pub fn get_all_sequences(reference_genome: Rc<Box<dyn ReferenceGenome>>) -> Vec<Sequence> {
    let mut all: Vec<Sequence> = Vec::new();

    if let Some(autosomes) = reference_genome.autosomes() {
        all.extend(autosomes);
    }

    if let Some(sex_chromosomes) = reference_genome.sex_chromosomes() {
        all.extend(sex_chromosomes);
    }

    if let Some(mitochondrion_chromosome) = reference_genome.mitochondrion_chromosome() {
        all.push(mitochondrion_chromosome);
    }

    if let Some(ebv_chromosome) = reference_genome.ebv_chromosome() {
        all.push(ebv_chromosome);
    }

    if let Some(unlocalized_sequences) = reference_genome.unlocalized_sequences() {
        all.extend(unlocalized_sequences);
    }

    if let Some(unplaced_sequences) = reference_genome.unplaced_sequences() {
        all.extend(unplaced_sequences);
    }

    if let Some(decoy_sequences) = reference_genome.decoy_sequences() {
        all.extend(decoy_sequences);
    }

    all
}

//====================//
// Types of sequences //
//====================//

#[derive(Clone)]
pub enum SequenceKind {
    Chromosome,
    Mitochondrion,
    EpsteinBarrVirus,
    Unlocalized,
    Unplaced,
    Decoy,
}

impl FromStr for SequenceKind {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_ref() {
            "chromosome" => Ok(SequenceKind::Chromosome),
            "mitochondrion" => Ok(SequenceKind::Mitochondrion),
            "ebv" => Ok(SequenceKind::EpsteinBarrVirus),
            "unlocalized" => Ok(SequenceKind::Unlocalized),
            "unplaced" => Ok(SequenceKind::Unplaced),
            "decoy" => Ok(SequenceKind::Decoy),
            s => Err(format!("Unknown sequence kind: {}", s)),
        }
    }
}

//===========//
// Sequences //
//===========//

#[derive(Clone)]
pub struct Sequence {
    name: &'static str,
    #[allow(dead_code)]
    kind: SequenceKind,
}

impl Sequence {
    pub fn new(name: &'static str, kind: SequenceKind) -> Self {
        Self { name, kind }
    }

    pub fn name(&self) -> &str {
        self.name
    }
}

/// Expands the provided arguments into a new Sequence.
#[macro_export]
macro_rules! sequence {
    ($name:expr,$kind:expr) => {
        Sequence::new($name, $kind.parse().unwrap())
    };
}

//===================//
// Reference Genomes //
//===================//

pub trait ReferenceGenome: Debug {
    fn name(&self) -> &'static str;
    fn autosomes(&self) -> Option<Vec<Sequence>>;
    fn sex_chromosomes(&self) -> Option<Vec<Sequence>>;
    fn mitochondrion_chromosome(&self) -> Option<Sequence>;
    fn ebv_chromosome(&self) -> Option<Sequence>;
    fn unlocalized_sequences(&self) -> Option<Vec<Sequence>>;
    fn unplaced_sequences(&self) -> Option<Vec<Sequence>>;
    fn decoy_sequences(&self) -> Option<Vec<Sequence>>;
}
