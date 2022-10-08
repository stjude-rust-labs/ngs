//! Functionality related to parsing Illumina read names.
//!
//! Only Illumina 1.4 and 1.8 read names are supported. The expected convention
//! is the following pattern, where `[]` denotes optional sections of the name.
//!
//! `INSTRUMENT:[RUN:FLOWCELL:]LANE:TILE:X:Y`
//!
use std::str::FromStr;

/// An Illumina read name.
#[derive(Debug)]
pub struct IlluminaReadName {
    /// The name of the instrument that produced this read.
    pub instrument_name: String,

    /// The id of the run that produced this read, if Illumina 1.8 read names.
    pub run: Option<String>,

    /// The id of the flowcell that produced this read, if Illumina 1.8 read names.
    pub flowcell: Option<String>,

    /// The lane of the flowcell that produced this read.
    pub lane: String,

    /// The of tile of the flowcell lane that produced this read.
    pub tile: String,

    /// The of X position of the flowcell well that produced this read.
    pub x: String,

    /// The of Y position of the flowcell well that produced this read.
    pub y: String,
}

impl FromStr for IlluminaReadName {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let segments: Vec<&str> = s.split(':').collect();
        let num_segments = segments.len();

        match num_segments {
            5 => Ok(IlluminaReadName {
                instrument_name: segments[0].into(),
                run: None,
                flowcell: None,
                lane: segments[1].into(),
                tile: segments[2].into(),
                x: segments[3].into(),
                y: segments[4].into(),
            }),
            7 => Ok(IlluminaReadName {
                instrument_name: segments[0].into(),
                run: Some(segments[1].into()),
                flowcell: Some(segments[2].into()),
                lane: segments[3].into(),
                tile: segments[4].into(),
                x: segments[5].into(),
                y: segments[6].into(),
            }),
            _ => Err(String::from("invalid number of segments for read name")),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    pub fn test_parse_illumina_1_4_fmt() {
        let result = "MACHINE:0:1234:55555:66666".parse::<IlluminaReadName>();
        assert!(result.is_ok());

        let read_name = result.unwrap();
        assert_eq!(read_name.instrument_name, "MACHINE");
        assert_eq!(read_name.run, None);
        assert_eq!(read_name.flowcell, None);
        assert_eq!(read_name.lane, "0");
        assert_eq!(read_name.tile, "1234");
        assert_eq!(read_name.x, "55555");
        assert_eq!(read_name.y, "66666");
    }

    #[test]
    pub fn test_parse_illumina_1_8_fmt() {
        let result = "MACHINE:0:A0A00AAAA:0:1234:55555:66666".parse::<IlluminaReadName>();
        assert!(result.is_ok());

        let read_name = result.unwrap();
        assert_eq!(read_name.instrument_name, "MACHINE");
        assert_eq!(read_name.run, Some("0".into()));
        assert_eq!(read_name.flowcell, Some("A0A00AAAA".into()));
        assert_eq!(read_name.lane, "0");
        assert_eq!(read_name.tile, "1234");
        assert_eq!(read_name.x, "55555");
        assert_eq!(read_name.y, "66666");
    }

    #[test]
    pub fn test_parse_failure() {
        let result = "MACHINE:0:".parse::<IlluminaReadName>();
        assert!(result.is_err());
    }
}
