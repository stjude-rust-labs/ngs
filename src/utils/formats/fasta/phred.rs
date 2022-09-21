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

macro_rules! phred_encoding {
    ($name:ident,$min_score:literal,$max_score:literal,$offset:literal) => {
        #[derive(Debug, Eq, Ord, PartialEq, PartialOrd)]
        struct $name(usize);

        impl $name {
            pub const MIN: Self = $name($min_score);
            pub const MAX: Self = $name($max_score);
            pub const ASCII_OFFSET: usize = $offset;

            pub fn score(&self) -> usize {
                self.0
            }

            #[allow(dead_code)]
            pub fn probability(&self) -> f64 {
                f64::powf(10.0, -(self.score() as f64) / 10.0)
            }
        }

        impl TryFrom<usize> for $name {
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
                    Ok($name(value))
                }
            }
        }

        impl From<$name> for char {
            /// SAFETY: phred score is guarenteed a range of MIN to MAX, and both
            /// $offset + MIN and $offset + MAX are guaranteed to be in the ASCII range.
            fn from(phred: $name) -> Self {
                unsafe { char::from_u32_unchecked(($name::ASCII_OFFSET + phred.score()) as u32) }
            }
        }

        impl TryFrom<char> for $name {
            type Error = PhredConversionError;

            fn try_from(c: char) -> Result<Self, Self::Error> {
                $name::try_from(c as usize - Self::ASCII_OFFSET)
            }
        }
    };
}

phred_encoding!(Illumina1Point3PhredScore, 0, 40, 64);
phred_encoding!(Illumina1Point8PhredScore, 0, 41, 33);

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_is_valid_from_a_lowest_score() {
        //==============//
        // Illumina 1.3 //
        //==============//

        let result = Illumina1Point3PhredScore::try_from(0);
        assert!(result.is_ok());

        let phred = result.unwrap();
        assert_eq!(phred.score(), 0);

        //==============//
        // Illumina 1.8 //
        //==============//

        let result = Illumina1Point8PhredScore::try_from(0);
        assert!(result.is_ok());

        let phred = result.unwrap();
        assert_eq!(phred.score(), 0);
    }

    #[test]
    fn it_is_not_valid_from_a_score_thats_too_high() {
        //==============//
        // Illumina 1.3 //
        //==============//

        let result = Illumina1Point3PhredScore::try_from(41);
        assert!(result.is_err());

        let PhredConversionError { reason } = result.unwrap_err();
        assert_eq!(reason, "score must be at most 40");

        //==============//
        // Illumina 1.8 //
        //==============//

        let result = Illumina1Point8PhredScore::try_from(42);
        assert!(result.is_err());

        let PhredConversionError { reason } = result.unwrap_err();
        assert_eq!(reason, "score must be at most 41");
    }

    #[test]
    fn it_correctly_converts_to_chars() {
        //==============//
        // Illumina 1.3 //
        //==============//

        let min_char = char::from(Illumina1Point3PhredScore(0));
        assert_eq!(min_char, '@');

        let max_char: char = Illumina1Point3PhredScore(40).into();
        assert_eq!(max_char, 'h');

        //==============//
        // Illumina 1.8 //
        //==============//

        let min_char = char::from(Illumina1Point8PhredScore(0));
        assert_eq!(min_char, '!');

        let max_char: char = Illumina1Point8PhredScore(41).into();
        assert_eq!(max_char, 'J');
    }

    #[test]
    fn it_correctly_converts_from_chars() {
        //==============//
        // Illumina 1.3 //
        //==============//

        let result = Illumina1Point3PhredScore::try_from('@');
        assert!(result.is_ok());

        let min_char = result.unwrap();
        assert_eq!(min_char.score(), 0);

        let result = Illumina1Point3PhredScore::try_from('h');
        assert!(result.is_ok());

        let max_char = result.unwrap();
        assert_eq!(max_char.score(), 40);

        let result = Illumina1Point3PhredScore::try_from('i');
        assert!(result.is_err());

        let PhredConversionError { reason } = result.unwrap_err();
        assert_eq!(reason, "score must be at most 40");

        //==============//
        // Illumina 1.8 //
        //==============//

        let result = Illumina1Point8PhredScore::try_from('!');
        assert!(result.is_ok());

        let min_char = result.unwrap();
        assert_eq!(min_char.score(), 0);

        let result = Illumina1Point8PhredScore::try_from('J');
        assert!(result.is_ok());

        let max_char = result.unwrap();
        assert_eq!(max_char.score(), 41);

        let result = Illumina1Point8PhredScore::try_from('K');
        assert!(result.is_err());

        let PhredConversionError { reason } = result.unwrap_err();
        assert_eq!(reason, "score must be at most 41");
    }

    #[test]
    fn it_correctly_computes_probabilities() {
        //==============//
        // Illumina 1.3 //
        //==============//

        assert!((Illumina1Point3PhredScore(0).probability() - 1.0).abs() < f64::EPSILON);
        assert!((Illumina1Point3PhredScore(10).probability() - 0.1).abs() < f64::EPSILON);
        assert!((Illumina1Point3PhredScore(20).probability() - 0.01).abs() < f64::EPSILON);
        assert!((Illumina1Point3PhredScore(30).probability() - 0.001).abs() < f64::EPSILON);
        assert!((Illumina1Point3PhredScore(40).probability() - 0.0001).abs() < f64::EPSILON);

        //==============//
        // Illumina 1.8 //
        //==============//

        assert!((Illumina1Point8PhredScore(0).probability() - 1.0).abs() < f64::EPSILON);
        assert!((Illumina1Point8PhredScore(10).probability() - 0.1).abs() < f64::EPSILON);
        assert!((Illumina1Point8PhredScore(20).probability() - 0.01).abs() < f64::EPSILON);
        assert!((Illumina1Point8PhredScore(30).probability() - 0.001).abs() < f64::EPSILON);
        assert!((Illumina1Point8PhredScore(40).probability() - 0.0001).abs() < f64::EPSILON);
    }
}
