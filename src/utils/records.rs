//! Extensions to and utilities concerning [`noodles::sam::record`]s.

//====================================//
// Command line parsing utility types //
//====================================//

use tracing::debug;

/// Utility enum to designate whether we are reviewing all records in the file
/// or just some of them.
pub enum NumberOfRecords {
    /// Designates that we should review _all_ of the records in the file.
    All,

    /// Designates that we should review _some_ of the records in the file. The
    /// exact count of records is stored in the `usize`.
    Some(usize),
}

impl From<Option<usize>> for NumberOfRecords {
    fn from(num_records: Option<usize>) -> Self {
        match num_records {
            Some(n) => {
                debug!("Reading a maximum of {} records.", n);
                NumberOfRecords::Some(n)
            }
            None => {
                debug!("Reading all available records.");
                NumberOfRecords::All
            }
        }
    }
}
