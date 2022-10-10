//! Version 2.0 of the T2T-CHM13 reference genome. This reference genome is
//! described in detail at
//! <https://www.ncbi.nlm.nih.gov/assembly/GCF_009914755.1/>.
//!
//! Link:
//! <https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz>

use crate::sequence;
use crate::utils::genome::ReferenceGenome;
use crate::utils::genome::Sequence;

/// Main struct for the T2T-CHM13 reference genome.
#[derive(Debug)]
#[allow(non_camel_case_types)]
pub struct T2T_CHM13;

impl ReferenceGenome for T2T_CHM13 {
    fn name(&self) -> &'static str {
        "T2T-CHM13"
    }

    fn version(&self) -> Option<&'static str> {
        Some("v2.0")
    }

    fn source(&self) -> &'static str {
        "T2T Consortium"
    }

    fn basis(&self) -> crate::utils::genome::GenomeBasis {
        crate::utils::genome::GenomeBasis::T2tChm13
    }

    fn patch(&self) -> Option<usize> {
        None
    }

    fn url(&self) -> Option<&'static str> {
        Some("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz")
    }

    fn autosomes(&self) -> Option<Vec<crate::utils::genome::Sequence>> {
        Some(vec![
            sequence!("NC_060925.1", "chromosome"), // chr1
            sequence!("NC_060926.1", "chromosome"), // chr2
            sequence!("NC_060927.1", "chromosome"), // chr3
            sequence!("NC_060928.1", "chromosome"), // chr4
            sequence!("NC_060929.1", "chromosome"), // chr5
            sequence!("NC_060930.1", "chromosome"), // chr6
            sequence!("NC_060931.1", "chromosome"), // chr7
            sequence!("NC_060932.1", "chromosome"), // chr8
            sequence!("NC_060933.1", "chromosome"), // chr9
            sequence!("NC_060934.1", "chromosome"), // chr10
            sequence!("NC_060935.1", "chromosome"), // chr11
            sequence!("NC_060936.1", "chromosome"), // chr12
            sequence!("NC_060937.1", "chromosome"), // chr13
            sequence!("NC_060938.1", "chromosome"), // chr14
            sequence!("NC_060939.1", "chromosome"), // chr15
            sequence!("NC_060940.1", "chromosome"), // chr16
            sequence!("NC_060941.1", "chromosome"), // chr17
            sequence!("NC_060942.1", "chromosome"), // chr18
            sequence!("NC_060943.1", "chromosome"), // chr19
            sequence!("NC_060944.1", "chromosome"), // chr20
            sequence!("NC_060945.1", "chromosome"), // chr21
            sequence!("NC_060946.1", "chromosome"), // chr22
        ])
    }

    fn sex_chromosomes(&self) -> Option<Vec<crate::utils::genome::Sequence>> {
        Some(vec![
            sequence!("NC_060947.1", "chromosome"), // chrX
            sequence!("NC_060948.1", "chromosome"), // chrY
        ])
    }

    fn mitochondrion_chromosome(&self) -> Option<crate::utils::genome::Sequence> {
        None
    }

    fn ebv_chromosome(&self) -> Option<crate::utils::genome::Sequence> {
        None
    }

    fn unlocalized_sequences(&self) -> Option<Vec<crate::utils::genome::Sequence>> {
        None
    }

    fn unplaced_sequences(&self) -> Option<Vec<crate::utils::genome::Sequence>> {
        None
    }

    fn decoy_sequences(&self) -> Option<Vec<crate::utils::genome::Sequence>> {
        None
    }
}

#[cfg(test)]
mod tests {

    use std::rc::Rc;

    use crate::utils::genome::get_primary_assembly;

    use super::*;

    #[test]
    pub fn it_has_the_correct_number_of_autosomes() {
        let grch38_no_alt = T2T_CHM13;
        assert_eq!(grch38_no_alt.autosomes().unwrap().len(), 22);
    }

    #[test]
    pub fn it_has_the_correct_number_of_sex_chromosomes() {
        let grch38_no_alt = T2T_CHM13;
        assert_eq!(grch38_no_alt.sex_chromosomes().unwrap().len(), 2);
    }

    #[test]
    pub fn it_has_the_correct_number_of_sequence_in_the_primary_assembly() {
        let grch38_no_alt: Rc<Box<dyn ReferenceGenome>> = Rc::new(Box::new(T2T_CHM13));
        assert_eq!(get_primary_assembly(grch38_no_alt).len(), 24);
    }

    #[test]
    pub fn it_contains_the_mitochodrion_chromosome() {
        let grch38_no_alt = T2T_CHM13;
        assert!(grch38_no_alt.mitochondrion_chromosome().is_none());
    }

    #[test]
    pub fn it_contains_the_epstein_barr_virus() {
        let grch38_no_alt = T2T_CHM13;
        assert!(grch38_no_alt.ebv_chromosome().is_none());
    }

    #[test]
    pub fn it_has_the_correct_number_of_unlocalized_sequences() {
        let grch38_no_alt = T2T_CHM13;
        assert!(grch38_no_alt.unlocalized_sequences().is_none());
    }

    #[test]
    pub fn it_has_the_correct_number_of_unplaced_sequences() {
        let grch38_no_alt = T2T_CHM13;
        assert!(grch38_no_alt.unplaced_sequences().is_none());
    }

    #[test]
    pub fn it_has_the_correct_number_of_decoy_sequences() {
        let grch38_no_alt = T2T_CHM13;
        assert!(grch38_no_alt.decoy_sequences().is_none());
    }
}
