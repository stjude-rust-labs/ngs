use std::{convert::TryFrom, error, fmt};

//========================//
// Phred Conversion Error //
//========================//

#[derive(Debug)]
struct PhredConversionError {
    reason: String,
}

impl PhredConversionError {
    pub fn new<S>(reason: S) -> Self
    where
        S: Into<String>,
    {
        PhredConversionError {
            reason: reason.into(),
        }
    }
}

impl fmt::Display for PhredConversionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "invalid conversion into phred score: {}", self.reason)
    }
}

impl error::Error for PhredConversionError {}

//============================//
// Illumina 1.8 Phred Scoring //
//============================//

#[derive(Debug, Eq, Ord, PartialEq, PartialOrd)]
struct Illumina1Point8Score(usize);

impl Illumina1Point8Score {
    pub const MIN: Self = Illumina1Point8Score(0);
    pub const MAX: Self = Illumina1Point8Score(41);
    pub const ASCII_OFFSET: usize = 33;

    pub fn score(&self) -> usize {
        self.0
    }

    #[allow(dead_code)]
    pub fn probability(&self) -> f64 {
        f64::powf(10.0, -(self.score() as f64) / 10.0)
    }
}

impl TryFrom<usize> for Illumina1Point8Score {
    type Error = PhredConversionError;

    fn try_from(value: usize) -> Result<Self, Self::Error> {
        if value < Self::MIN.score() {
            Err(PhredConversionError::new(format!(
                "score must be at least {}",
                Self::MIN.score()
            )))
        } else if value > Self::MAX.score() {
            Err(PhredConversionError::new(format!(
                "score must be at most {}",
                Self::MAX.score()
            )))
        } else {
            Ok(Illumina1Point8Score(value))
        }
    }
}

impl From<Illumina1Point8Score> for char {
    fn from(phred: Illumina1Point8Score) -> Self {
        // SAFETY: phred score is guarenteed a range of MIN to MAX, and both
        // 33 + MIN and 33 + MAX are guaranteed to be in the ASCII range.
        unsafe {
            char::from_u32_unchecked((Illumina1Point8Score::ASCII_OFFSET + phred.score()) as u32)
        }
    }
}

impl TryFrom<char> for Illumina1Point8Score {
    type Error = PhredConversionError;

    fn try_from(c: char) -> Result<Self, Self::Error> {
        Illumina1Point8Score::try_from(c as usize - Self::ASCII_OFFSET)
    }
}

#[cfg(test)]
mod tests {
    use super::{Illumina1Point8Score, PhredConversionError};

    #[test]
    fn it_is_valid_from_a_zero_score() {
        let result = Illumina1Point8Score::try_from(0);
        assert!(result.is_ok());

        let phred = result.unwrap();
        assert_eq!(phred.score(), 0);
    }

    #[test]
    fn it_is_not_valid_from_a_score_thats_too_high() {
        let result = Illumina1Point8Score::try_from(74);
        assert!(result.is_err());

        let PhredConversionError { reason } = result.unwrap_err();
        assert_eq!(reason, "score must be at most 41");
    }

    #[test]
    fn it_correctly_converts_to_chars() {
        let min_char = char::from(Illumina1Point8Score(0));
        assert_eq!(min_char, '!');

        let max_char: char = Illumina1Point8Score(41).into();
        assert_eq!(max_char, 'J');
    }

    #[test]
    fn it_correctly_converts_from_chars() {
        let result = Illumina1Point8Score::try_from('!');
        assert!(result.is_ok());

        let min_char = result.unwrap();
        assert_eq!(min_char.score(), 0);

        let result = Illumina1Point8Score::try_from('J');
        assert!(result.is_ok());

        let max_char = result.unwrap();
        assert_eq!(max_char.score(), 41);

        let result = Illumina1Point8Score::try_from('K');
        assert!(result.is_err());

        let PhredConversionError { reason } = result.unwrap_err();
        assert_eq!(reason, "score must be at most 41");
    }

    #[test]
    fn it_correctly_computes_probabilities() {
        assert!((Illumina1Point8Score(0).probability() - 1.0).abs() < f64::EPSILON);
        assert!((Illumina1Point8Score(10).probability() - 0.1).abs() < f64::EPSILON);
        assert!((Illumina1Point8Score(20).probability() - 0.01).abs() < f64::EPSILON);
        assert!((Illumina1Point8Score(30).probability() - 0.001).abs() < f64::EPSILON);
        assert!((Illumina1Point8Score(40).probability() - 0.0001).abs() < f64::EPSILON);
    }
}
