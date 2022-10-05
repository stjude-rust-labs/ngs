use std::{
    collections::HashMap,
    sync::{
        atomic::{AtomicUsize, Ordering},
        Mutex,
    },
};

use noodles_fasta::record::Sequence;
use serde::{Deserialize, Serialize};

use crate::lib::utils::hashmap::MutexBackedHashMap;

//==========================//
// Read Designation Metrics //
//==========================//

#[derive(Debug, Default)]
pub struct AtomicReadDesignationMetrics {
    pub primary: AtomicUsize,
    pub secondary: AtomicUsize,
    pub supplementary: AtomicUsize,
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct ReadDesignationMetrics {
    primary: usize,
    secondary: usize,
    supplementary: usize,
}

impl AtomicReadDesignationMetrics {
    pub fn clone_nonatomic(&self) -> ReadDesignationMetrics {
        ReadDesignationMetrics {
            primary: self.primary.load(Ordering::SeqCst),
            secondary: self.secondary.load(Ordering::SeqCst),
            supplementary: self.supplementary.load(Ordering::SeqCst),
        }
    }
}

//================//
// Record Metrics //
//================//

#[derive(Debug, Default)]
pub struct AtomicRecordMetrics {
    pub total: AtomicUsize,
    pub unmapped: AtomicUsize,
    pub duplicate: AtomicUsize,
    pub designation: AtomicReadDesignationMetrics,
    pub primary_mapped: AtomicUsize,
    pub primary_duplicate: AtomicUsize,
    pub paired: AtomicUsize,
    pub read_1: AtomicUsize,
    pub read_2: AtomicUsize,
    pub proper_pair: AtomicUsize,
    pub singleton: AtomicUsize,
    pub mate_mapped: AtomicUsize,
    pub mate_reference_sequence_id_mismatch: AtomicUsize,
    pub mate_reference_sequence_id_mismatch_hq: AtomicUsize,
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct RecordMetrics {
    total: usize,
    unmapped: usize,
    duplicate: usize,
    designation: ReadDesignationMetrics,
    primary_mapped: usize,
    primary_duplicate: usize,
    paired: usize,
    read_1: usize,
    read_2: usize,
    proper_pair: usize,
    singleton: usize,
    mate_mapped: usize,
    mate_reference_sequence_id_mismatch: usize,
    mate_reference_sequence_id_mismatch_hq: usize,
}

impl AtomicRecordMetrics {
    pub fn clone_nonatomic(&self) -> RecordMetrics {
        RecordMetrics {
            total: self.total.load(Ordering::SeqCst),
            unmapped: self.unmapped.load(Ordering::SeqCst),
            duplicate: self.duplicate.load(Ordering::SeqCst),
            designation: self.designation.clone_nonatomic(),
            primary_mapped: self.primary_mapped.load(Ordering::SeqCst),
            primary_duplicate: self.primary_duplicate.load(Ordering::SeqCst),
            paired: self.paired.load(Ordering::SeqCst),
            read_1: self.read_1.load(Ordering::SeqCst),
            read_2: self.read_2.load(Ordering::SeqCst),
            proper_pair: self.proper_pair.load(Ordering::SeqCst),
            singleton: self.singleton.load(Ordering::SeqCst),
            mate_mapped: self.mate_mapped.load(Ordering::SeqCst),
            mate_reference_sequence_id_mismatch: self
                .mate_reference_sequence_id_mismatch
                .load(Ordering::SeqCst),
            mate_reference_sequence_id_mismatch_hq: self
                .mate_reference_sequence_id_mismatch_hq
                .load(Ordering::SeqCst),
        }
    }
}

//===============//
// Cigar Metrics //
//===============//

#[derive(Debug, Default)]
pub struct AtomicCigarMetrics {
    // Cigar operation pileups
    pub read_one_cigar_ops: Mutex<HashMap<String, usize>>,
    pub read_two_cigar_ops: Mutex<HashMap<String, usize>>,
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct CigarMetrics {
    // Cigar operation pileups
    read_one_cigar_ops: HashMap<String, usize>,
    read_two_cigar_ops: HashMap<String, usize>,
}

impl AtomicCigarMetrics {
    pub fn clone_nonatomic(&self) -> CigarMetrics {
        CigarMetrics {
            read_one_cigar_ops: self.read_one_cigar_ops.lock().unwrap().clone(),
            read_two_cigar_ops: self.read_two_cigar_ops.lock().unwrap().clone(),
        }
    }
}

//=================//
// Summary Metrics //
//=================//

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct SummaryMetrics {
    pub duplication_pct: f64,
    pub unmapped_pct: f64,
    pub mate_reference_sequence_id_mismatch_pct: f64,
    pub mate_reference_sequence_id_mismatch_hq_pct: f64,
}

//==============//
// Main Metrics //
//==============//

#[derive(Debug, Default)]
pub struct AtomicGeneralMetrics {
    pub records: AtomicRecordMetrics,
    pub cigar: AtomicCigarMetrics,
    pub summary: Option<SummaryMetrics>,
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct GeneralMetrics {
    records: RecordMetrics,
    cigar: CigarMetrics,
    summary: Option<SummaryMetrics>,
}

impl AtomicGeneralMetrics {
    pub fn clone_nonatomic(&self) -> GeneralMetrics {
        GeneralMetrics {
            records: self.records.clone_nonatomic(),
            cigar: self.cigar.clone_nonatomic(),
            summary: self.summary.clone(),
        }
    }
}

#[derive(Debug, Default)]
pub struct GeneralMetricsFacet {
    pub metrics: AtomicGeneralMetrics,
    pub reference_sequences: HashMap<String, Sequence>,
}
