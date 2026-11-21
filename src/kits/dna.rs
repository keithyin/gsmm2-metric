//! gskits `dna` 模块的本地移植(仅保留本项目用到的 reverse_complement)。

pub static COMPLEMENT_TABLE: [u8; 256] = {
    let mut table = [0u8; 256];
    table[b'A' as usize] = b'T';
    table[b'T' as usize] = b'A';
    table[b'C' as usize] = b'G';
    table[b'G' as usize] = b'C';
    table[b'-' as usize] = b'-';
    table[b'*' as usize] = b'*';
    table[b'N' as usize] = b'N';

    table[b'a' as usize] = b't';
    table[b't' as usize] = b'a';
    table[b'c' as usize] = b'g';
    table[b'g' as usize] = b'c';
    table[b'-' as usize] = b'-';
    table[b'*' as usize] = b'*';
    table[b'n' as usize] = b'n';

    table
};

/// b"ACGTAA" -> b"TTACGT"
pub fn reverse_complement(dna: &[u8]) -> Vec<u8> {
    dna.iter()
        .rev()
        .map(|&base| COMPLEMENT_TABLE[base as usize])
        .collect()
}

#[cfg(test)]
mod test {
    use super::reverse_complement;

    #[test]
    fn test_reverse_complement() {
        let dna_sequence = "ATCGTAGC";
        let res = reverse_complement(dna_sequence.as_bytes());
        assert_eq!(res, b"GCTACGAT");
    }
}
