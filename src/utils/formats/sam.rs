//! Utilities related to opening and manipulating SAM files.

use noodles::sam;
use regex::{Captures, Regex};

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
