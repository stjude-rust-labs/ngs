//! This reference genome represents the hg38 no alt analysis set provided by
//! Microsoft Genomics Service. Currently, this reference genome is not publicly
//! available for download at any link—you must email Microsoft to get a copy of
//! it.

use crate::sequence;
use crate::utils::genome::GenomeBasis;
use crate::utils::genome::ReferenceGenome;
use crate::utils::genome::Sequence;

/// Main struct for the hg38m1x reference genome.
#[derive(Debug)]
pub struct HG38M1X;

impl ReferenceGenome for HG38M1X {
    fn name(&self) -> &'static str {
        "hg38m1x"
    }

    fn version(&self) -> Option<&'static str> {
        None
    }

    fn source(&self) -> &'static str {
        "Microsoft Genomics"
    }

    fn basis(&self) -> GenomeBasis {
        GenomeBasis::GRCh38
    }

    fn patch(&self) -> Option<usize> {
        // Based on the NCBI's GRCh38_no_alt_AnalysisSet, which is based on the
        // first patch (patch 0) of GRCh38.
        Some(0)
    }

    fn url(&self) -> Option<&'static str> {
        // Download from the Azure Open Datasets bucket. You can learn more
        // about which reference genomes they host at
        // https://learn.microsoft.com/en-us/azure/open-datasets/dataset-human-genomes?tabs=azure-storage
        Some("https://datasetmsgen.blob.core.windows.net/dataset/hg38m1x/hg38m1x.fa?sv=2020-10-02&si=prod&sr=c&sig=rIK0VRz7ZbYdcNhvwMMc3M1hr%2BgI0zgEG0Shbblc3nc%3D")
    }

    fn autosomes(&self) -> Option<Vec<Sequence>> {
        Some(vec![
            sequence!("chr1", "chromosome"),
            sequence!("chr2", "chromosome"),
            sequence!("chr3", "chromosome"),
            sequence!("chr4", "chromosome"),
            sequence!("chr5", "chromosome"),
            sequence!("chr6", "chromosome"),
            sequence!("chr7", "chromosome"),
            sequence!("chr8", "chromosome"),
            sequence!("chr9", "chromosome"),
            sequence!("chr10", "chromosome"),
            sequence!("chr11", "chromosome"),
            sequence!("chr12", "chromosome"),
            sequence!("chr13", "chromosome"),
            sequence!("chr14", "chromosome"),
            sequence!("chr15", "chromosome"),
            sequence!("chr16", "chromosome"),
            sequence!("chr17", "chromosome"),
            sequence!("chr18", "chromosome"),
            sequence!("chr19", "chromosome"),
            sequence!("chr20", "chromosome"),
            sequence!("chr21", "chromosome"),
            sequence!("chr22", "chromosome"),
        ])
    }

    fn sex_chromosomes(&self) -> Option<Vec<Sequence>> {
        Some(vec![
            sequence!("chrX", "chromosome"),
            sequence!("chrY", "chromosome"),
        ])
    }

    fn mitochondrion_chromosome(&self) -> Option<Sequence> {
        Some(sequence!("chrM", "mitochondrion"))
    }

    fn alternative_contig_sequences(&self) -> Option<Vec<Sequence>> {
        None
    }

    fn ebv_chromosome(&self) -> Option<Sequence> {
        None
    }

    fn unlocalized_sequences(&self) -> Option<Vec<Sequence>> {
        Some(vec![
            sequence!("chr1_KI270706v1_random", "unlocalized"),
            sequence!("chr1_KI270707v1_random", "unlocalized"),
            sequence!("chr1_KI270708v1_random", "unlocalized"),
            sequence!("chr1_KI270709v1_random", "unlocalized"),
            sequence!("chr1_KI270710v1_random", "unlocalized"),
            sequence!("chr1_KI270711v1_random", "unlocalized"),
            sequence!("chr1_KI270712v1_random", "unlocalized"),
            sequence!("chr1_KI270713v1_random", "unlocalized"),
            sequence!("chr1_KI270714v1_random", "unlocalized"),
            sequence!("chr2_KI270715v1_random", "unlocalized"),
            sequence!("chr2_KI270716v1_random", "unlocalized"),
            sequence!("chr3_GL000221v1_random", "unlocalized"),
            sequence!("chr4_GL000008v2_random", "unlocalized"),
            sequence!("chr5_GL000208v1_random", "unlocalized"),
            sequence!("chr9_KI270717v1_random", "unlocalized"),
            sequence!("chr9_KI270718v1_random", "unlocalized"),
            sequence!("chr9_KI270719v1_random", "unlocalized"),
            sequence!("chr9_KI270720v1_random", "unlocalized"),
            sequence!("chr11_KI270721v1_random", "unlocalized"),
            sequence!("chr14_GL000009v2_random", "unlocalized"),
            sequence!("chr14_GL000225v1_random", "unlocalized"),
            sequence!("chr14_KI270722v1_random", "unlocalized"),
            sequence!("chr14_GL000194v1_random", "unlocalized"),
            sequence!("chr14_KI270723v1_random", "unlocalized"),
            sequence!("chr14_KI270724v1_random", "unlocalized"),
            sequence!("chr14_KI270725v1_random", "unlocalized"),
            sequence!("chr14_KI270726v1_random", "unlocalized"),
            sequence!("chr15_KI270727v1_random", "unlocalized"),
            sequence!("chr16_KI270728v1_random", "unlocalized"),
            sequence!("chr17_GL000205v2_random", "unlocalized"),
            sequence!("chr17_KI270729v1_random", "unlocalized"),
            sequence!("chr17_KI270730v1_random", "unlocalized"),
            sequence!("chr22_KI270731v1_random", "unlocalized"),
            sequence!("chr22_KI270732v1_random", "unlocalized"),
            sequence!("chr22_KI270733v1_random", "unlocalized"),
            sequence!("chr22_KI270734v1_random", "unlocalized"),
            sequence!("chr22_KI270735v1_random", "unlocalized"),
            sequence!("chr22_KI270736v1_random", "unlocalized"),
            sequence!("chr22_KI270737v1_random", "unlocalized"),
            sequence!("chr22_KI270738v1_random", "unlocalized"),
            sequence!("chr22_KI270739v1_random", "unlocalized"),
            sequence!("chrY_KI270740v1_random", "unlocalized"),
        ])
    }

    fn unplaced_sequences(&self) -> Option<Vec<Sequence>> {
        Some(vec![
            sequence!("chrUn_KI270302v1", "unplaced"),
            sequence!("chrUn_KI270304v1", "unplaced"),
            sequence!("chrUn_KI270303v1", "unplaced"),
            sequence!("chrUn_KI270305v1", "unplaced"),
            sequence!("chrUn_KI270322v1", "unplaced"),
            sequence!("chrUn_KI270320v1", "unplaced"),
            sequence!("chrUn_KI270310v1", "unplaced"),
            sequence!("chrUn_KI270316v1", "unplaced"),
            sequence!("chrUn_KI270315v1", "unplaced"),
            sequence!("chrUn_KI270312v1", "unplaced"),
            sequence!("chrUn_KI270311v1", "unplaced"),
            sequence!("chrUn_KI270317v1", "unplaced"),
            sequence!("chrUn_KI270412v1", "unplaced"),
            sequence!("chrUn_KI270411v1", "unplaced"),
            sequence!("chrUn_KI270414v1", "unplaced"),
            sequence!("chrUn_KI270419v1", "unplaced"),
            sequence!("chrUn_KI270418v1", "unplaced"),
            sequence!("chrUn_KI270420v1", "unplaced"),
            sequence!("chrUn_KI270424v1", "unplaced"),
            sequence!("chrUn_KI270417v1", "unplaced"),
            sequence!("chrUn_KI270422v1", "unplaced"),
            sequence!("chrUn_KI270423v1", "unplaced"),
            sequence!("chrUn_KI270425v1", "unplaced"),
            sequence!("chrUn_KI270429v1", "unplaced"),
            sequence!("chrUn_KI270442v1", "unplaced"),
            sequence!("chrUn_KI270466v1", "unplaced"),
            sequence!("chrUn_KI270465v1", "unplaced"),
            sequence!("chrUn_KI270467v1", "unplaced"),
            sequence!("chrUn_KI270435v1", "unplaced"),
            sequence!("chrUn_KI270438v1", "unplaced"),
            sequence!("chrUn_KI270468v1", "unplaced"),
            sequence!("chrUn_KI270510v1", "unplaced"),
            sequence!("chrUn_KI270509v1", "unplaced"),
            sequence!("chrUn_KI270518v1", "unplaced"),
            sequence!("chrUn_KI270508v1", "unplaced"),
            sequence!("chrUn_KI270516v1", "unplaced"),
            sequence!("chrUn_KI270512v1", "unplaced"),
            sequence!("chrUn_KI270519v1", "unplaced"),
            sequence!("chrUn_KI270522v1", "unplaced"),
            sequence!("chrUn_KI270511v1", "unplaced"),
            sequence!("chrUn_KI270515v1", "unplaced"),
            sequence!("chrUn_KI270507v1", "unplaced"),
            sequence!("chrUn_KI270517v1", "unplaced"),
            sequence!("chrUn_KI270529v1", "unplaced"),
            sequence!("chrUn_KI270528v1", "unplaced"),
            sequence!("chrUn_KI270530v1", "unplaced"),
            sequence!("chrUn_KI270539v1", "unplaced"),
            sequence!("chrUn_KI270538v1", "unplaced"),
            sequence!("chrUn_KI270544v1", "unplaced"),
            sequence!("chrUn_KI270548v1", "unplaced"),
            sequence!("chrUn_KI270583v1", "unplaced"),
            sequence!("chrUn_KI270587v1", "unplaced"),
            sequence!("chrUn_KI270580v1", "unplaced"),
            sequence!("chrUn_KI270581v1", "unplaced"),
            sequence!("chrUn_KI270579v1", "unplaced"),
            sequence!("chrUn_KI270589v1", "unplaced"),
            sequence!("chrUn_KI270590v1", "unplaced"),
            sequence!("chrUn_KI270584v1", "unplaced"),
            sequence!("chrUn_KI270582v1", "unplaced"),
            sequence!("chrUn_KI270588v1", "unplaced"),
            sequence!("chrUn_KI270593v1", "unplaced"),
            sequence!("chrUn_KI270591v1", "unplaced"),
            sequence!("chrUn_KI270330v1", "unplaced"),
            sequence!("chrUn_KI270329v1", "unplaced"),
            sequence!("chrUn_KI270334v1", "unplaced"),
            sequence!("chrUn_KI270333v1", "unplaced"),
            sequence!("chrUn_KI270335v1", "unplaced"),
            sequence!("chrUn_KI270338v1", "unplaced"),
            sequence!("chrUn_KI270340v1", "unplaced"),
            sequence!("chrUn_KI270336v1", "unplaced"),
            sequence!("chrUn_KI270337v1", "unplaced"),
            sequence!("chrUn_KI270363v1", "unplaced"),
            sequence!("chrUn_KI270364v1", "unplaced"),
            sequence!("chrUn_KI270362v1", "unplaced"),
            sequence!("chrUn_KI270366v1", "unplaced"),
            sequence!("chrUn_KI270378v1", "unplaced"),
            sequence!("chrUn_KI270379v1", "unplaced"),
            sequence!("chrUn_KI270389v1", "unplaced"),
            sequence!("chrUn_KI270390v1", "unplaced"),
            sequence!("chrUn_KI270387v1", "unplaced"),
            sequence!("chrUn_KI270395v1", "unplaced"),
            sequence!("chrUn_KI270396v1", "unplaced"),
            sequence!("chrUn_KI270388v1", "unplaced"),
            sequence!("chrUn_KI270394v1", "unplaced"),
            sequence!("chrUn_KI270386v1", "unplaced"),
            sequence!("chrUn_KI270391v1", "unplaced"),
            sequence!("chrUn_KI270383v1", "unplaced"),
            sequence!("chrUn_KI270393v1", "unplaced"),
            sequence!("chrUn_KI270384v1", "unplaced"),
            sequence!("chrUn_KI270392v1", "unplaced"),
            sequence!("chrUn_KI270381v1", "unplaced"),
            sequence!("chrUn_KI270385v1", "unplaced"),
            sequence!("chrUn_KI270382v1", "unplaced"),
            sequence!("chrUn_KI270376v1", "unplaced"),
            sequence!("chrUn_KI270374v1", "unplaced"),
            sequence!("chrUn_KI270372v1", "unplaced"),
            sequence!("chrUn_KI270373v1", "unplaced"),
            sequence!("chrUn_KI270375v1", "unplaced"),
            sequence!("chrUn_KI270371v1", "unplaced"),
            sequence!("chrUn_KI270448v1", "unplaced"),
            sequence!("chrUn_KI270521v1", "unplaced"),
            sequence!("chrUn_GL000195v1", "unplaced"),
            sequence!("chrUn_GL000219v1", "unplaced"),
            sequence!("chrUn_GL000220v1", "unplaced"),
            sequence!("chrUn_GL000224v1", "unplaced"),
            sequence!("chrUn_KI270741v1", "unplaced"),
            sequence!("chrUn_GL000226v1", "unplaced"),
            sequence!("chrUn_GL000213v1", "unplaced"),
            sequence!("chrUn_KI270743v1", "unplaced"),
            sequence!("chrUn_KI270744v1", "unplaced"),
            sequence!("chrUn_KI270745v1", "unplaced"),
            sequence!("chrUn_KI270746v1", "unplaced"),
            sequence!("chrUn_KI270747v1", "unplaced"),
            sequence!("chrUn_KI270748v1", "unplaced"),
            sequence!("chrUn_KI270749v1", "unplaced"),
            sequence!("chrUn_KI270750v1", "unplaced"),
            sequence!("chrUn_KI270751v1", "unplaced"),
            sequence!("chrUn_KI270752v1", "unplaced"),
            sequence!("chrUn_KI270753v1", "unplaced"),
            sequence!("chrUn_KI270754v1", "unplaced"),
            sequence!("chrUn_KI270755v1", "unplaced"),
            sequence!("chrUn_KI270756v1", "unplaced"),
            sequence!("chrUn_KI270757v1", "unplaced"),
            sequence!("chrUn_GL000214v1", "unplaced"),
            sequence!("chrUn_KI270742v1", "unplaced"),
            sequence!("chrUn_GL000216v2", "unplaced"),
            sequence!("chrUn_GL000218v1", "unplaced"),
        ])
    }

    fn decoy_sequences(&self) -> Option<Vec<Sequence>> {
        None
    }

    fn other_sequences(&self) -> Option<Vec<Sequence>> {
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
        let hg38m1x = HG38M1X;
        assert_eq!(hg38m1x.autosomes().unwrap().len(), 22);
    }

    #[test]
    pub fn it_has_the_correct_number_of_sex_chromosomes() {
        let hg38m1x = HG38M1X;
        assert_eq!(hg38m1x.sex_chromosomes().unwrap().len(), 2);
    }

    #[test]
    pub fn it_has_the_correct_number_of_sequence_in_the_primary_assembly() {
        let hg38m1x: Rc<Box<dyn ReferenceGenome>> = Rc::new(Box::new(HG38M1X));
        assert_eq!(get_primary_assembly(hg38m1x).len(), 193);
    }

    #[test]
    pub fn it_contains_the_mitochodrion_chromosome() {
        let hg38m1x = HG38M1X;
        assert!(hg38m1x.mitochondrion_chromosome().is_some());
    }

    #[test]
    pub fn it_does_not_contain_any_alternative_contigs() {
        let hg38m1x = HG38M1X;
        assert!(hg38m1x.alternative_contig_sequences().is_none());
    }

    #[test]
    pub fn it_does_not_contain_the_epstein_barr_virus() {
        let hg38m1x = HG38M1X;
        assert!(hg38m1x.ebv_chromosome().is_none());
    }

    #[test]
    pub fn it_has_the_correct_number_of_unlocalized_sequences() {
        let hg38m1x = HG38M1X;
        assert_eq!(hg38m1x.unlocalized_sequences().unwrap().len(), 42);
    }

    #[test]
    pub fn it_has_the_correct_number_of_unplaced_sequences() {
        let hg38m1x = HG38M1X;
        assert_eq!(hg38m1x.unplaced_sequences().unwrap().len(), 127);
    }

    #[test]
    pub fn it_has_the_correct_number_of_decoy_sequences() {
        let hg38m1x = HG38M1X;
        assert!(hg38m1x.decoy_sequences().is_none());
    }

    #[test]
    pub fn it_does_not_contain_any_other_sequences() {
        let hg38m1x = HG38M1X;
        assert!(hg38m1x.other_sequences().is_none());
    }
}
