use noodles::sam::header;
use std::collections::HashSet;
use std::sync::Arc;
use tracing::warn;

use lazy_static::lazy_static;

// Strings used to index into the HashMaps used to store the Read Group ordering flags.
// Lazy statics are used to save memory.
lazy_static! {
    /// String used to index into the HashMaps used to store the "overall" ordering flags.
    pub static ref OVERALL: Arc<String> = Arc::new(String::from("overall"));

    /// String used to index into th HashMaps used to store the "unknown_read_group" ordering flags.
    pub static ref UNKNOWN_READ_GROUP: Arc<String> = Arc::new(String::from("unknown_read_group"));
}

/// Compares the read group tags found in the records
/// and the read groups found in the header.
/// Returns a vector of read group names that were found in the header
/// but not in the records.
pub fn validate_read_group_info(
    found_rgs: &HashSet<Arc<String>>,
    header: &header::Header,
) -> Vec<String> {
    let mut rgs_in_header_not_records = Vec::new();
    let mut rgs_in_records_not_header = Vec::new();

    for (rg_id, _) in header.read_groups() {
        if !found_rgs.contains(rg_id) {
            rgs_in_header_not_records.push(rg_id.to_string());
        }
    }
    if !rgs_in_header_not_records.is_empty() {
        warn!(
            "The following read groups were not found in the file: {:?}",
            rgs_in_header_not_records
        );
    }

    for rg_id in found_rgs {
        if !header.read_groups().contains_key(rg_id.as_str()) {
            rgs_in_records_not_header.push(rg_id.to_string());
        }
    }
    if !rgs_in_records_not_header.is_empty() {
        warn!(
            "The following read groups were not found in the header: {:?}",
            rgs_in_records_not_header
        );
    }

    rgs_in_header_not_records
}
