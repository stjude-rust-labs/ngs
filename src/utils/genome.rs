//! Utilities related to reference genomes.

#[cfg(feature = "extended-reference-genomes")]
pub mod microsoft;
pub mod ncbi;
#[cfg(feature = "extended-reference-genomes")]
pub mod one_thousand_genomes;
#[cfg(feature = "extended-reference-genomes")]
pub mod t2t_consortium;

use std::fmt;
use std::fmt::Debug;
use std::rc::Rc;
use std::str::FromStr;

#[cfg(feature = "extended-reference-genomes")]
use self::microsoft::hg38m1x::HG38M1X;
use self::ncbi::grch38_no_alt::GRCh38NoAltAnalysisSet;
#[cfg(feature = "extended-reference-genomes")]
use self::one_thousand_genomes::grch38_full_analysis_set_with_decoy_hla::GRCh38FullAnalysisSetWithDecoyHLA;
#[cfg(feature = "extended-reference-genomes")]
use self::one_thousand_genomes::hs37d5::HS37D5;
#[cfg(feature = "extended-reference-genomes")]
use self::t2t_consortium::t2t_chm13::T2T_CHM13;

//=================//
// Utility methods //
//=================//

/// Gets all of the supported reference genomes for the tool. When new reference
/// genomes are added, this needs to be updated.
pub fn get_all_reference_genomes() -> Vec<Box<dyn ReferenceGenome>> {
    vec![
        #[cfg(feature = "extended-reference-genomes")]
        Box::new(HS37D5),
        #[cfg(feature = "extended-reference-genomes")]
        Box::new(GRCh38FullAnalysisSetWithDecoyHLA),
        Box::new(GRCh38NoAltAnalysisSet),
        #[cfg(feature = "extended-reference-genomes")]
        Box::new(HG38M1X),
        #[cfg(feature = "extended-reference-genomes")]
        Box::new(T2T_CHM13),
    ]
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

/// Gets the primary assembly for a given reference genome (as defined by the autosomes,
/// the sex chromosomes, the alternative_contigs, the unlocalized sequences, and the
/// unplaced sequences).
pub fn get_primary_assembly(reference_genome: Rc<Box<dyn ReferenceGenome>>) -> Vec<Sequence> {
    let mut primary: Vec<Sequence> = Vec::new();

    if let Some(autosomes) = reference_genome.autosomes() {
        primary.extend(autosomes);
    }

    if let Some(sex_chromosomes) = reference_genome.sex_chromosomes() {
        primary.extend(sex_chromosomes);
    }

    if let Some(alternative_contig_sequences) = reference_genome.alternative_contig_sequences() {
        primary.extend(alternative_contig_sequences);
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

    if let Some(alternative_contig_sequences) = reference_genome.alternative_contig_sequences() {
        all.extend(alternative_contig_sequences);
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

    if let Some(other_sequences) = reference_genome.other_sequences() {
        all.extend(other_sequences);
    }

    all
}

//====================//
// Types of sequences //
//====================//

/// A kind of sequence that can be included in a reference genome.
#[derive(Clone)]
pub enum SequenceKind {
    /// A canonical sequence containing species-specific DNA.
    Chromosome,

    /// A sequence representing mitochondrial DNA.
    Mitochondrion,

    /// An alternative contig.
    AlternativeContig,

    /// A sequence reprsenting the Epstein-Barr virus.
    EpsteinBarrVirus,

    /// A sequence representing an unlocalized contig. It is generally known
    /// which chromosome these belong to, but not their exact position within
    /// the chromosome.
    Unlocalized,

    /// A sequence representing an unplaced contig. These are generally known to
    /// be part of a species's DNA, but not their exact position or chromosome.
    Unplaced,

    /// A sequence that is not part of the species canonical DNA. These are
    /// typically inserted to clean up results in variant calling and remove
    /// contaminants.
    Decoy,

    /// A sequence that doesn't fall neatly into the categories above.
    Other,
}

impl FromStr for SequenceKind {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_ref() {
            "chromosome" => Ok(SequenceKind::Chromosome),
            "mitochondrion" => Ok(SequenceKind::Mitochondrion),
            "alt" => Ok(SequenceKind::AlternativeContig),
            "ebv" => Ok(SequenceKind::EpsteinBarrVirus),
            "unlocalized" => Ok(SequenceKind::Unlocalized),
            "unplaced" => Ok(SequenceKind::Unplaced),
            "decoy" => Ok(SequenceKind::Decoy),
            "other" => Ok(SequenceKind::Other),
            s => Err(format!("Unknown sequence kind: {}", s)),
        }
    }
}

//===========//
// Sequences //
//===========//

/// A sequence contained within a reference genome, including the name and the
/// kind of the sequence.
#[derive(Clone)]
pub struct Sequence {
    name: &'static str,
    #[allow(dead_code)]
    kind: SequenceKind,
}

impl Sequence {
    /// Creates a new [`Sequence`].
    pub fn new(name: &'static str, kind: SequenceKind) -> Self {
        Self { name, kind }
    }

    /// Gives the name of the [`Sequence`].
    pub fn name(&self) -> &str {
        self.name
    }
}

/// Expands the provided arguments into a new Sequence. This is provided for
/// convenience when specifying large reference genomes.
#[macro_export]
macro_rules! sequence {
    ($name:expr,$kind:expr) => {
        Sequence::new($name, $kind.parse().unwrap())
    };
}

//==============//
// Genome Basis //
//==============//

/// The basis upon which a reference genome is created. These are generally
/// revisions of the Genome Reference Consortiums regular releases of the
/// reference genome.
#[derive(Debug)]
pub enum GenomeBasis {
    /// GRCh37-based reference genomes
    GRCh37,

    /// GRCh38-based reference genomes
    GRCh38,

    /// T2T-based reference genomes
    T2tChm13,
}

impl fmt::Display for GenomeBasis {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::GRCh37 => write!(f, "GRCh37"),
            Self::GRCh38 => write!(f, "GRCh38"),
            Self::T2tChm13 => write!(f, "T2tChm13"),
        }
    }
}

//===================//
// Reference Genomes //
//===================//

/// A struct is [`ReferenceGenome`] if it represents a reference genome that is
/// supported by this tool.
pub trait ReferenceGenome: Debug {
    /// Name of the reference genome.
    fn name(&self) -> &'static str;

    /// If available, the version of the reference genome.
    fn version(&self) -> Option<&'static str>;

    /// Center that produced this reference genome.
    fn source(&self) -> &'static str;

    /// Build upon which this reference genome is based.
    fn basis(&self) -> GenomeBasis;

    /// If available, the patch of the reference genome upon which this genome is based.
    fn patch(&self) -> Option<usize>;

    /// If available, url where this reference genome is available.
    fn url(&self) -> Option<&'static str>;

    /// A triplet id uniquely identifies a reference genome with respect to its
    /// colloquial name, the version of the release, and the reference genome it
    /// is based off of.
    ///
    /// # Examples
    ///
    /// ```
    /// use ngs::utils::genome::ncbi::grch38_no_alt::GRCh38NoAltAnalysisSet;
    /// // Trait must be in scope to use the `triplet_id` method.
    /// use crate::ngs::utils::genome::ReferenceGenome;
    /// let reference = GRCh38NoAltAnalysisSet;
    /// assert_eq!(reference.triplet_id(), "GRCh38_no_alt_AnalysisSet-none-GRCh38.p0")
    /// ```
    fn triplet_id(&self) -> String {
        // (1) Colloquial name
        let mut result = String::from(self.name());
        result.push('-');

        // (2) Version, if available
        let version = self.version().unwrap_or("none");
        result.push_str(version);
        result.push('-');

        // (3) Basis, with patch if available
        result.push_str(&self.basis().to_string());
        if let Some(patch) = self.patch() {
            result.push_str(".p");
            result.push_str(&patch.to_string());
        }

        result
    }

    /// If available, the autosomes included in this reference genome.
    fn autosomes(&self) -> Option<Vec<Sequence>>;

    /// If available, the sex chromosomes included in this reference genome.
    fn sex_chromosomes(&self) -> Option<Vec<Sequence>>;

    /// If available, the mitochondrial DNA included in this reference genome.
    fn mitochondrion_chromosome(&self) -> Option<Sequence>;

    /// Alternative contigs that were included in this reference genome.
    fn alternative_contig_sequences(&self) -> Option<Vec<Sequence>>;

    /// If available, the epstein-barr virus DNA included in this reference
    /// genome.
    fn ebv_chromosome(&self) -> Option<Sequence>;

    /// If available, a number of unlocalized sequences included in this
    /// reference genome.
    fn unlocalized_sequences(&self) -> Option<Vec<Sequence>>;

    /// If available, a number of unplaced sequences included in this reference
    /// genome.
    fn unplaced_sequences(&self) -> Option<Vec<Sequence>>;

    /// If available, any decoy sequences included in this reference genome.
    fn decoy_sequences(&self) -> Option<Vec<Sequence>>;

    /// Sequences that don't fit neatly into the categories included above but that
    /// should still be considered as part of the reference genome.
    fn other_sequences(&self) -> Option<Vec<Sequence>>;
}
