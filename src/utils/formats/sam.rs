//! Utilities related to opening and manipulating Sequence Alignment Map (SAM) files.

use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::path::Path;

use anyhow::bail;
use noodles::sam;
use noodles::sam::Header;
use regex::Captures;
use regex::Regex;

use super::BioinformaticsFileFormat;

//=================//
// Utility Methods //
//=================//

/// Corrects common header mistakes. See the inline comments for the things that
/// are automatically corrected.
pub fn correct_common_header_mistakes(header: String) -> String {
    // (1) Corrects any lowercase platform units in the read group to be all
    // uppercase. This is especially important for data that contains 'illumina'
    // instead of the correct 'ILLUMINA'.
    let pattern = Regex::new("(\tPL:)(.+)").unwrap();
    let replaced = pattern.replace_all(&header, |c: &Captures<'_>| {
        format!("{}{}", &c[1], c[2].to_uppercase())
    });

    replaced.to_string()
}

/// Parses a SAM/BAM/CRAM header from a string while also correcting common
/// header mistakes.
pub fn parse_header(header: String) -> sam::Header {
    correct_common_header_mistakes(header)
        .parse()
        .expect("Could not parse SAM/BAM/CRAM header.")
}

//====================================//
// Sequence Alignment Map (SAM) files //
//====================================//

/// Attempts to open a SAM file from a given source. Note that this file is private
/// because it should never be called by an external module (use [`open_and_parse`]
/// instead).
fn open<P>(src: P) -> anyhow::Result<sam::Reader<Box<dyn BufRead>>>
where
    P: AsRef<Path>,
{
    let path = src.as_ref();
    let file = File::open(path);

    match BioinformaticsFileFormat::try_detect(path) {
        Some(BioinformaticsFileFormat::SAM) => {
            let reader = file.map(BufReader::new)?;
            Ok(sam::Reader::new(Box::new(reader)))
        }
        Some(format) => bail!("incompatible formats: required SAM, found {}", format),
        None => {
            let ext = path
                .extension()
                .expect("file extension to exist")
                .to_str()
                .expect("extension to be convertible to &str");
            bail!("Not able to determine filetype for extension: {}", ext)
        }
    }
}

/// Contains both the opened SAM file and the parsed header from the SAM file.
pub struct ParsedSAMFile(sam::Reader<Box<dyn BufRead>>, Header);

/// Opens and subsequently parses a SAM file's header. This is useful when opening SAM
/// files when you want the corrections applied by [`correct_common_header_mistakes`] to
/// apply.
pub fn open_and_parse<P>(src: P) -> anyhow::Result<ParsedSAMFile>
where
    P: AsRef<Path>,
{
    let mut sam = open(src)?;
    let header_as_text = sam.read_header()?;
    let header = parse_header(header_as_text);
    Ok(ParsedSAMFile(sam, header))
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    pub fn test_illumina_lowercase_fix() {
        let data = "@RG\tID:rg0\tPL:illumina\n";
        let expected = "@RG\tID:rg0\tPL:ILLUMINA\n";

        assert_eq!(correct_common_header_mistakes(data.to_string()), expected);
    }
}
