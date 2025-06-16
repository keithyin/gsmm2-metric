use mm2::minimap2::{Aligner, Mapping};

pub fn get_target_start_end(
    ori_start: usize,
    ori_end: usize,
    target_len: usize,
    is_rev: bool,
) -> (usize, usize) {
    if is_rev {
        (target_len - ori_end, target_len - ori_start)
    } else {
        (ori_start, ori_end)
    }
}

pub fn get_target_substr<'a>(
    ori_start: usize,
    ori_end: usize,
    is_rev: bool,
    target_fwd: &'a str,
    target_rev: &'a str,
) -> &'a str {
    let (start, end) = get_target_start_end(ori_start, ori_end, target_fwd.len(), is_rev);
    if is_rev {
        &target_rev[start..end]
    } else {
        &target_fwd[start..end]
    }
}

pub fn do_align_4_homo(ori_query_seq: &str, qstart: usize, qend: usize, target_substr: &str) -> Option<Mapping> {
    let mut aligner = Aligner::builder()
        .map_ont()
        .with_cigar() // cigar_str has bug in minimap2="0.1.20+minimap2.2.28"
        .with_sam_out()
        .with_index_threads(1);

    aligner.mapopt.best_n = 1;
    aligner.idxopt.k = 11;
    aligner.idxopt.w = 1;

    aligner.mapopt.a = 4;
    aligner.mapopt.b = 10;
    aligner.mapopt.q = 4;
    aligner.mapopt.q2 = 24;
    aligner.mapopt.e = 2;
    aligner.mapopt.e2 = 1;

    let aligner = aligner
        .with_seq_and_id(target_substr.as_bytes(), "target".as_bytes())
        .unwrap();
    let mut hits = aligner
        .map(
            ori_query_seq[qstart..qend].as_bytes(),
            false,
            false,
            None,
            Some(&[0x4000000, 0x40000000]),
            Some(b"query"),
        )
        .unwrap()
        .into_iter()
        .filter(|hit| hit.is_primary)
        .map(|v| Some(v))
        .collect::<Vec<_>>();
    
    if hits.is_empty() {
        None
    } else {
        hits[0].take()
    }
}
