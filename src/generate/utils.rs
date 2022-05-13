/// Characters:
/// 'A' => 0x41, 'a' => 0x61
/// 'C' => 0x43, 'c' => 0x63
/// 'G' => 0x47, 'g' => 0x67
/// 'T' => 0x54, 't' => 0x67
fn compliment(c: &u8) -> Option<u8> {
    match c {
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

pub fn reverse_compliment(seq: &[u8]) -> Option<Vec<u8>> {
    let iter = seq.iter().map(compliment);

    if iter.clone().any(|x| x.is_none()) {
        return None;
    }

    let mut result: Vec<u8> = iter.map(|x| x.unwrap()).collect();
    result.reverse(); // then reverse

    Some(result)
}

mod tests {

    #[test]
    fn test_compliment_valid() {
        let input: Vec<u8> = "ACTGactg".as_bytes().to_vec();
        assert_eq!(
            super::reverse_compliment(&input),
            Some(Vec::from("cagtCAGT"))
        );
    }

    #[test]
    fn test_compliment_invalid() {
        let input: Vec<u8> = "n".as_bytes().to_vec();
        assert_eq!(super::reverse_compliment(&input), None);
    }
}
