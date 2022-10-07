use noodles_fastq as fastq;
use rand_distr::Normal;

/// Simple utility struct used for passing around a sequence name and its
/// associated length. This is useful for, say, generating a random chromosome
/// and random coordinates within that chromosome to generate a read pair from.
#[derive(Debug)]
pub struct SeqLen(pub String, pub usize);

impl SeqLen {
    /// Gets the sequence name
    pub fn get_seq_name(&self) -> &String {
        &self.0
    }

    /// Gets the length of the sequence
    pub fn get_seq_len(&self) -> usize {
        self.1
    }
}

#[derive(Debug)]
pub struct NormalDistributionParams(f64, f64);

impl NormalDistributionParams {
    pub fn new(mu: f64, sigma: f64) -> Self {
        NormalDistributionParams(mu, sigma)
    }

    pub fn get_mu(&self) -> f64 {
        self.0
    }

    pub fn get_sigma(&self) -> f64 {
        self.1
    }
}

impl From<NormalDistributionParams> for Normal<f64> {
    fn from(params: NormalDistributionParams) -> Self {
        Normal::new(params.get_mu(), params.get_sigma()).unwrap_or_else(|_| {
            panic!(
                "Could not create normal distribution from parameters: {:?}",
                params
            )
        })
    }
}

/// Struct representing a pair of fastq::Records: one for the forward read and
/// one for the reverse read. This is useful for returning the read pair to a
/// caller of the `SequenceProvider.generate_read_pair` method.
#[derive(Debug)]
pub struct PairedRead(pub fastq::Record, pub fastq::Record);

impl PairedRead {
    /// Gets the forward read.
    pub fn get_forward_read(&self) -> &fastq::Record {
        &self.0
    }

    /// Gets the reverse read.
    pub fn get_reverse_read(&self) -> &fastq::Record {
        &self.1
    }
}

/// Utility method to compliment a bytes string. This method is considered to be
/// temporary until Michael implements a reverse compliment method within
/// noodles (see
/// https://github.com/zaeleus/noodles/issues/86#issuecomment-1124383295) for
/// more details.
///
/// # Arguments
///
/// * `seq`: the sequence of bytes to compliment.
///
/// # Examples
///
/// ```
/// use ngs::generate::utils;
/// let base = 'a' as u8;
/// let compliment = utils::compliment(&base).unwrap();
/// assert_eq!(compliment, 't' as u8)
/// ```
pub fn compliment(seq: &u8) -> Option<u8> {
    match seq {
        0x61 => Some(0x74), // 'a' => 't'
        0x63 => Some(0x67), // 'c' => 'g'
        0x67 => Some(0x63), // 'g' => 'c'
        0x74 => Some(0x61), // 't' => 'a'
        0x41 => Some(0x54), // 'A' => 'T'
        0x43 => Some(0x47), // 'C' => 'G'
        0x47 => Some(0x43), // 'G' => 'C'
        0x54 => Some(0x41), // 'T' => 'A'
        _ => None,
    }
}

/// Reverse compliments a byte string, failing if any of the characters in the
/// string fail to be complimented. This is wrapped in an Option accordingly.
///
/// # Arguments
///
/// * `seq`: the sequence of bytes to reverse compliment.
pub fn reverse_compliment(seq: &[u8]) -> Option<Vec<u8>> {
    let iter = seq.iter().map(compliment);

    if iter.clone().any(|x| x.is_none()) {
        return None;
    }

    let mut result: Vec<u8> = iter.map(|x| x.unwrap()).collect();
    result.reverse(); // then reverse

    Some(result)
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_compliment_valid() {
        let input = "ACTGactg".as_bytes();
        assert_eq!(reverse_compliment(input), Some(Vec::from("cagtCAGT")));
    }

    #[test]
    fn test_compliment_invalid() {
        let input: Vec<u8> = "n".as_bytes().to_vec();
        assert_eq!(super::reverse_compliment(&input), None);
    }
}
