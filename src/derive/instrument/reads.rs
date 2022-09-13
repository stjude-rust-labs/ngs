use std::process;
use tracing::error;

use crate::errors::NGSExitCodes;

#[derive(Debug)]
pub struct IlluminaReadName {
    pub instrument_name: String,
    pub run: Option<String>,
    pub flowcell: Option<String>,
    pub lane: String,
    pub tile: String,
    pub x: String,
    pub y: String,
}

impl IlluminaReadName {
    /// Parses the machine name from the standard Illumina read naming convention.
    /// The read naming convention is as follows:
    ///
    /// ```
    /// INSTRUMENT:RUN:FLOWCELL:LANE:TILE:X:Y
    /// ```
    ///
    pub fn from(name: String) -> Self {
        let parts: Vec<&str> = name.split(':').collect();
        let num_parts = parts.len();

        if num_parts == 5 {
            IlluminaReadName {
                instrument_name: parts[0].to_string(),
                run: None,
                flowcell: None,
                lane: parts[1].to_string(),
                tile: parts[2].to_string(),
                x: parts[3].to_string(),
                y: parts[4].to_string(),
            }
        } else if num_parts == 7 {
            IlluminaReadName {
                instrument_name: parts[0].to_string(),
                run: Some(parts[1].to_string()),
                flowcell: Some(parts[2].to_string()),
                lane: parts[3].to_string(),
                tile: parts[4].to_string(),
                x: parts[5].to_string(),
                y: parts[6].to_string(),
            }
        } else {
            error!(
                "Could not parse Illumina-formatted query names for read: {}.",
                name
            );
            error!(
                "Illumina-formatted reads are expected to be colon (:) delimited with \
             either five or seven fields. Please see the Wikipedia page on \
             FASTQ files (https://en.wikipedia.org/wiki/FASTQ_format#Illumina_sequence_identifiers) \
             for more details."
            );

            process::exit(NGSExitCodes::InvalidInputData as i32);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    pub fn test_parse_illumina_1_4_fmt() {
        let result = IlluminaReadName::from("MACHINE:0:1234:55555:66666".to_string());
        assert_eq!(result.instrument_name, "MACHINE");
        assert_eq!(result.run, None);
        assert_eq!(result.flowcell, None);
        assert_eq!(result.lane, "0".to_string());
        assert_eq!(result.tile, "1234".to_string());
        assert_eq!(result.x, "55555".to_string());
        assert_eq!(result.y, "66666".to_string());
    }

    #[test]
    pub fn test_parse_illumina_1_8_fmt() {
        let result = IlluminaReadName::from("MACHINE:0:A0A00AAAA:0:1234:55555:66666".to_string());
        assert_eq!(result.instrument_name, "MACHINE");
        assert_eq!(result.run, Some("0".to_string()));
        assert_eq!(result.flowcell, Some("A0A00AAAA".to_string()));
        assert_eq!(result.lane, "0".to_string());
        assert_eq!(result.tile, "1234".to_string());
        assert_eq!(result.x, "55555".to_string());
        assert_eq!(result.y, "66666".to_string());
    }
}
