//! This reference genome is the reference genome used for most of phase 2 +
//! phase 3 of the 1000 genomes project. It contains the GRCh37 primary assembly
//! (chromosomal plus unlocalized and unplaced contigs), the rCRS mitochondrial
//! sequence (AC: NC_012920), the Human herpesvirus 4 type 1 (AC:NC_007605)
//! [otherwise known as the Epstein-Barr Virus], and a set of concatenated decoy
//! sequences.
//!
//! For more information, see the README for this reference genome at:
//! http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/README_human_reference_20110707
//!
//! Link:
//! https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz

use crate::{
    sequence,
    utils::genome::{GenomeBasis, ReferenceGenome, Sequence},
};

#[derive(Debug)]
pub struct HS37D5;

impl ReferenceGenome for HS37D5 {
    fn name(&self) -> &'static str {
        "hs37d5"
    }

    fn source(&self) -> &'static str {
        "1000 Genomes Project"
    }

    fn basis(&self) -> GenomeBasis {
        GenomeBasis::GRCh37
    }

    fn url(&self) -> &'static str {
        "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz"
    }

    fn autosomes(&self) -> Option<Vec<Sequence>> {
        Some(vec![
            sequence!("1", "chromosome"),
            sequence!("2", "chromosome"),
            sequence!("3", "chromosome"),
            sequence!("4", "chromosome"),
            sequence!("5", "chromosome"),
            sequence!("6", "chromosome"),
            sequence!("7", "chromosome"),
            sequence!("8", "chromosome"),
            sequence!("9", "chromosome"),
            sequence!("10", "chromosome"),
            sequence!("11", "chromosome"),
            sequence!("12", "chromosome"),
            sequence!("13", "chromosome"),
            sequence!("14", "chromosome"),
            sequence!("15", "chromosome"),
            sequence!("16", "chromosome"),
            sequence!("17", "chromosome"),
            sequence!("18", "chromosome"),
            sequence!("19", "chromosome"),
            sequence!("20", "chromosome"),
            sequence!("21", "chromosome"),
            sequence!("22", "chromosome"),
        ])
    }

    fn sex_chromosomes(&self) -> Option<Vec<Sequence>> {
        Some(vec![
            sequence!("X", "chromosome"),
            sequence!("Y", "chromosome"),
        ])
    }

    fn mitochondrion_chromosome(&self) -> Option<Sequence> {
        Some(sequence!("MT", "mitochondrion"))
    }

    fn ebv_chromosome(&self) -> Option<Sequence> {
        Some(sequence!("NC_007605", "ebv"))
    }

    fn unlocalized_sequences(&self) -> Option<Vec<Sequence>> {
        Some(vec![
            sequence!("GL000207.1", "unlocalized"), // chr18_gl000207_random
            sequence!("GL000210.1", "unlocalized"), // chr21_gl000210_random
            sequence!("GL000201.1", "unlocalized"), // chr9_gl000201_random
            sequence!("GL000197.1", "unlocalized"), // chr8_gl000197_random
            sequence!("GL000203.1", "unlocalized"), // chr17_gl000203_random
            sequence!("GL000196.1", "unlocalized"), // chr8_gl000196_random
            sequence!("GL000202.1", "unlocalized"), // chr11_gl000202_random
            sequence!("GL000206.1", "unlocalized"), // chr17_gl000206_random
            sequence!("GL000204.1", "unlocalized"), // chr17_gl000204_random
            sequence!("GL000198.1", "unlocalized"), // chr9_gl000198_random
            sequence!("GL000208.1", "unlocalized"), // chr19_gl000208_random
            sequence!("GL000191.1", "unlocalized"), // chr1_gl000191_random
            sequence!("GL000209.1", "unlocalized"), // chr19_gl000209_random
            sequence!("GL000199.1", "unlocalized"), // chr9_gl000199_random
            sequence!("GL000205.1", "unlocalized"), // chr17_gl000205_random
            sequence!("GL000195.1", "unlocalized"), // chr7_gl000195_random
            sequence!("GL000200.1", "unlocalized"), // chr9_gl000200_random
            sequence!("GL000193.1", "unlocalized"), // chr4_gl000193_random
            sequence!("GL000194.1", "unlocalized"), // chr4_gl000194_random
            sequence!("GL000192.1", "unlocalized"), // chr1_gl000192_random
        ])
    }

    fn unplaced_sequences(&self) -> Option<Vec<Sequence>> {
        Some(vec![
            sequence!("GL000226.1", "unplaced"), // chrUn_gl000226
            sequence!("GL000229.1", "unplaced"), // chrUn_gl000229
            sequence!("GL000231.1", "unplaced"), // chrUn_gl000231
            sequence!("GL000239.1", "unplaced"), // chrUn_gl000239
            sequence!("GL000235.1", "unplaced"), // chrUn_gl000235
            sequence!("GL000247.1", "unplaced"), // chrUn_gl000247
            sequence!("GL000245.1", "unplaced"), // chrUn_gl000245
            sequence!("GL000246.1", "unplaced"), // chrUn_gl000246
            sequence!("GL000249.1", "unplaced"), // chrUn_gl000249
            sequence!("GL000248.1", "unplaced"), // chrUn_gl000248
            sequence!("GL000244.1", "unplaced"), // chrUn_gl000244
            sequence!("GL000238.1", "unplaced"), // chrUn_gl000238
            sequence!("GL000234.1", "unplaced"), // chrUn_gl000234
            sequence!("GL000232.1", "unplaced"), // chrUn_gl000232
            sequence!("GL000240.1", "unplaced"), // chrUn_gl000240
            sequence!("GL000236.1", "unplaced"), // chrUn_gl000236
            sequence!("GL000241.1", "unplaced"), // chrUn_gl000241
            sequence!("GL000243.1", "unplaced"), // chrUn_gl000243
            sequence!("GL000242.1", "unplaced"), // chrUn_gl000242
            sequence!("GL000230.1", "unplaced"), // chrUn_gl000230
            sequence!("GL000237.1", "unplaced"), // chrUn_gl000237
            sequence!("GL000233.1", "unplaced"), // chrUn_gl000233
            sequence!("GL000227.1", "unplaced"), // chrUn_gl000227
            sequence!("GL000228.1", "unplaced"), // chrUn_gl000228
            sequence!("GL000214.1", "unplaced"), // chrUn_gl000214
            sequence!("GL000221.1", "unplaced"), // chrUn_gl000221
            sequence!("GL000218.1", "unplaced"), // chrUn_gl000218
            sequence!("GL000220.1", "unplaced"), // chrUn_gl000220
            sequence!("GL000213.1", "unplaced"), // chrUn_gl000213
            sequence!("GL000211.1", "unplaced"), // chrUn_gl000211
            sequence!("GL000217.1", "unplaced"), // chrUn_gl000217
            sequence!("GL000216.1", "unplaced"), // chrUn_gl000216
            sequence!("GL000215.1", "unplaced"), // chrUn_gl000215
            sequence!("GL000219.1", "unplaced"), // chrUn_gl000219
            sequence!("GL000224.1", "unplaced"), // chrUn_gl000224
            sequence!("GL000223.1", "unplaced"), // chrUn_gl000223
            sequence!("GL000212.1", "unplaced"), // chrUn_gl000212
            sequence!("GL000222.1", "unplaced"), // chrUn_gl000222
            sequence!("GL000225.1", "unplaced"), // chrUn_gl000225
        ])
    }

    fn decoy_sequences(&self) -> Option<Vec<Sequence>> {
        Some(vec![sequence!("hs37d5", "decoy")])
    }
}

#[cfg(test)]
mod tests {

    use std::rc::Rc;

    use crate::utils::genome::get_primary_assembly;

    use super::*;

    #[test]
    pub fn it_has_the_correct_number_of_autosomes() {
        let hs37d5 = HS37D5;
        assert_eq!(hs37d5.autosomes().unwrap().len(), 22);
    }

    #[test]
    pub fn it_has_the_correct_number_of_sex_chromosomes() {
        let hs37d5 = HS37D5;
        assert_eq!(hs37d5.sex_chromosomes().unwrap().len(), 2);
    }

    #[test]
    pub fn it_has_the_correct_number_of_sequence_in_the_primary_assembly() {
        let hs37d5: Rc<Box<dyn ReferenceGenome>> = Rc::new(Box::new(HS37D5));
        assert_eq!(get_primary_assembly(hs37d5).len(), 83);
    }

    #[test]
    pub fn it_contains_the_mitochodrion_chromosome() {
        let hs37d5 = HS37D5;
        assert!(hs37d5.mitochondrion_chromosome().is_some());
    }

    #[test]
    pub fn it_contains_the_epstein_barr_virus() {
        let hs37d5 = HS37D5;
        assert!(hs37d5.ebv_chromosome().is_some());
    }

    #[test]
    pub fn it_has_the_correct_number_of_unlocalized_sequences() {
        let hs37d5 = HS37D5;
        assert_eq!(hs37d5.unlocalized_sequences().unwrap().len(), 20);
    }

    #[test]
    pub fn it_has_the_correct_number_of_unplaced_sequences() {
        let hs37d5 = HS37D5;
        assert_eq!(hs37d5.unplaced_sequences().unwrap().len(), 39);
    }

    #[test]
    pub fn it_has_the_correct_number_of_decoy_sequences() {
        let hs37d5 = HS37D5;
        assert_eq!(hs37d5.decoy_sequences().unwrap().len(), 1);
    }
}
