//! Convenient, extended containers for bioinformatics formats supported by `noodles`.

use std::fmt::Display;
use std::path::PathBuf;

pub mod bam;
pub mod fasta;
pub mod fastq;
pub mod gff;
pub mod sam;

/// Represents all of the supported bioinformatics file formats that can be
/// detected by the extension of the filename.
#[derive(Debug, PartialEq, Eq)]
#[allow(non_camel_case_types)]
pub enum BioinformaticsFileFormat {
    /// A FASTA file.
    FASTA,

    /// A Gzipped FASTA file.
    FASTA_GZ,

    /// A FASTQ file.
    FASTQ,

    /// A Gzipped FASTQ file.
    FASTQ_GZ,

    /// A SAM file.
    SAM,

    /// A BAM file.
    BAM,

    /// A CRAM file.
    CRAM,

    /// A VCF file.
    VCF,

    /// A Gzipped VCF file.
    VCF_GZ,

    /// A BCF file.
    BCF,

    /// A GFF file.
    GFF,

    /// A Gzipped GFF file.
    GFF_GZ,

    /// A GTF file.
    GTF,

    /// A Gzipped GTF file.
    GTF_GZ,

    /// A BED file.
    BED,
}

impl Display for BioinformaticsFileFormat {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            BioinformaticsFileFormat::FASTA => write!(f, "FASTA"),
            BioinformaticsFileFormat::FASTA_GZ => write!(f, "Gzipped FASTA"),
            BioinformaticsFileFormat::FASTQ => write!(f, "FASTQ"),
            BioinformaticsFileFormat::FASTQ_GZ => write!(f, "Gzipped FASTQ"),
            BioinformaticsFileFormat::SAM => write!(f, "SAM"),
            BioinformaticsFileFormat::BAM => write!(f, "BAM"),
            BioinformaticsFileFormat::CRAM => write!(f, "CRAM"),
            BioinformaticsFileFormat::VCF => write!(f, "VCF"),
            BioinformaticsFileFormat::VCF_GZ => write!(f, "Gzipped VCF"),
            BioinformaticsFileFormat::BCF => write!(f, "BCF"),
            BioinformaticsFileFormat::GFF => write!(f, "GFF"),
            BioinformaticsFileFormat::GFF_GZ => write!(f, "Gzipped GFF"),
            BioinformaticsFileFormat::GTF => write!(f, "GTF"),
            BioinformaticsFileFormat::GTF_GZ => write!(f, "Gzipped GTF"),
            BioinformaticsFileFormat::BED => write!(f, "BED"),
        }
    }
}

impl BioinformaticsFileFormat {
    /// Tries to detect a bioinformatics file format from the extension of the
    /// provided `filepath`.
    pub fn try_detect(filepath: impl Into<PathBuf>) -> Option<BioinformaticsFileFormat> {
        let path: PathBuf = filepath.into();

        if let Some(ext) = path.extension() {
            // (1) First, check if the extension is "gz". If it is, then we need
            // a different set of matches.

            if ext.eq_ignore_ascii_case("gz") {
                let path_as_str = path.to_str().expect("path to be convertible to &str");

                if path_as_str.ends_with("fasta.gz")
                    || path_as_str.ends_with("fna.gz")
                    || path_as_str.ends_with("fa.gz")
                {
                    return Some(Self::FASTA_GZ);
                } else if path_as_str.ends_with("fq.gz") || path_as_str.ends_with("fastq.gz") {
                    return Some(Self::FASTQ_GZ);
                } else if path_as_str.ends_with("vcf.gz") {
                    return Some(Self::VCF_GZ);
                } else if path_as_str.ends_with("gff.gz") || path_as_str.ends_with("gff3.gz") {
                    return Some(Self::GFF_GZ);
                } else if path_as_str.ends_with("gtf.gz") {
                    return Some(Self::GTF_GZ);
                } else {
                    return None;
                }
            }

            // (2) Next, if the extension is _not_ "gz", then we can match using
            // a simple `match` statement.

            match ext
                .to_ascii_lowercase()
                .to_str()
                .expect("extension to be convertible to &str")
            {
                "fasta" | "fna" | "fa" => Some(Self::FASTA),
                "fastq" | "fq" => Some(Self::FASTQ),
                "sam" => Some(Self::SAM),
                "bam" => Some(Self::BAM),
                "cram" => Some(Self::CRAM),
                "vcf" => Some(Self::VCF),
                "bcf" => Some(Self::BCF),
                "gff" | "gff3" => Some(Self::GFF),
                "gtf" => Some(Self::GTF),
                "bed" => Some(Self::BED),
                _ => None,
            }
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::BioinformaticsFileFormat;

    #[test]
    fn it_correctly_identifies_fasta_files() {
        assert_eq!(
            BioinformaticsFileFormat::try_detect("sample.fasta"),
            Some(BioinformaticsFileFormat::FASTA)
        );
        assert_eq!(
            BioinformaticsFileFormat::try_detect("sample.fna"),
            Some(BioinformaticsFileFormat::FASTA)
        );
        assert_eq!(
            BioinformaticsFileFormat::try_detect("sample.fa"),
            Some(BioinformaticsFileFormat::FASTA)
        );
        assert_eq!(
            BioinformaticsFileFormat::try_detect("sample.fasta.gz"),
            Some(BioinformaticsFileFormat::FASTA_GZ)
        );
        assert_eq!(
            BioinformaticsFileFormat::try_detect("sample.fna.gz"),
            Some(BioinformaticsFileFormat::FASTA_GZ)
        );
        assert_eq!(
            BioinformaticsFileFormat::try_detect("sample.fa.gz"),
            Some(BioinformaticsFileFormat::FASTA_GZ)
        );
    }

    #[test]
    fn it_correctly_identifies_fastq_files() {
        assert_eq!(
            BioinformaticsFileFormat::try_detect("sample.fastq"),
            Some(BioinformaticsFileFormat::FASTQ)
        );
        assert_eq!(
            BioinformaticsFileFormat::try_detect("sample.fq"),
            Some(BioinformaticsFileFormat::FASTQ)
        );
        assert_eq!(
            BioinformaticsFileFormat::try_detect("sample.fastq.gz"),
            Some(BioinformaticsFileFormat::FASTQ_GZ)
        );
        assert_eq!(
            BioinformaticsFileFormat::try_detect("sample.fq.gz"),
            Some(BioinformaticsFileFormat::FASTQ_GZ)
        );
    }

    #[test]
    fn it_correctly_identifies_sam_files() {
        assert_eq!(
            BioinformaticsFileFormat::try_detect("sample.sam"),
            Some(BioinformaticsFileFormat::SAM)
        );
    }

    #[test]
    fn it_correctly_identifies_bam_files() {
        assert_eq!(
            BioinformaticsFileFormat::try_detect("sample.bam"),
            Some(BioinformaticsFileFormat::BAM)
        );
    }

    #[test]
    fn it_correctly_identifies_cram_files() {
        assert_eq!(
            BioinformaticsFileFormat::try_detect("sample.cram"),
            Some(BioinformaticsFileFormat::CRAM)
        );
    }

    #[test]
    fn it_correctly_identifies_vcf_files() {
        assert_eq!(
            BioinformaticsFileFormat::try_detect("sample.vcf"),
            Some(BioinformaticsFileFormat::VCF)
        );
        assert_eq!(
            BioinformaticsFileFormat::try_detect("sample.vcf.gz"),
            Some(BioinformaticsFileFormat::VCF_GZ)
        );
    }

    #[test]
    fn it_correctly_identifies_bcf_files() {
        assert_eq!(
            BioinformaticsFileFormat::try_detect("sample.bcf"),
            Some(BioinformaticsFileFormat::BCF)
        );
    }

    #[test]
    fn it_correctly_identifies_gff_files() {
        assert_eq!(
            BioinformaticsFileFormat::try_detect("sample.gff"),
            Some(BioinformaticsFileFormat::GFF)
        );
        assert_eq!(
            BioinformaticsFileFormat::try_detect("sample.gff.gz"),
            Some(BioinformaticsFileFormat::GFF_GZ)
        );
    }

    #[test]
    fn it_correctly_identifies_gtf_files() {
        assert_eq!(
            BioinformaticsFileFormat::try_detect("sample.gtf"),
            Some(BioinformaticsFileFormat::GTF)
        );
        assert_eq!(
            BioinformaticsFileFormat::try_detect("sample.gtf.gz"),
            Some(BioinformaticsFileFormat::GTF_GZ)
        );
    }

    #[test]
    fn it_correctly_identifies_bed_files() {
        assert_eq!(
            BioinformaticsFileFormat::try_detect("sample.bed"),
            Some(BioinformaticsFileFormat::BED)
        );
    }
}
