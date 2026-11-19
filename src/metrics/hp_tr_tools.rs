use std::{collections::HashMap, sync::Arc};

use mm2::minimap2::{Aligner, Mapping};

/// ref 位点 -> 覆盖该位点的 motif 区间，即 GlobalDataKey::TargetRegion2Motif4HpTr。
pub type RefPos2Motifs = HashMap<usize, Vec<((usize, usize), Arc<String>)>>;

/// 把比对端点收缩到 HP/TR motif 区域的边界。
///
/// 起点落在某 motif 区域内部时，起点后移到该区域的终点；终点落在区域内部时，
/// 终点前移到该区域的起点，使目标窗口不含残缺的 motif 区域。
///
/// 返回 `None` 表示收缩后窗口为空。最常见的原因是比对的两个端点都落在同一个
/// motif 区域内（整条 read 位于一段 poly-A/T 内部）。此时无法确定该区域的
/// 完整重复数，调用方必须跳过该 hit，不能继续切片。
pub fn shrink_to_motif_boundary(
    target_start: usize,
    target_end: usize,
    ref_pos2motifs: &RefPos2Motifs,
) -> Option<(usize, usize)> {
    if target_end == 0 {
        return None;
    }

    let mut start = target_start;
    let mut end = target_end;

    if let Some(regions) = ref_pos2motifs.get(&start) {
        let new_start = regions.iter().max_by_key(|v| v.0.1).unwrap().0.1;
        tracing::info!("old_target_start: {start} -> {new_start}");
        start = new_start;
    }

    if let Some(regions) = ref_pos2motifs.get(&(end - 1)) {
        let new_end = regions.iter().min_by_key(|v| v.0.0).unwrap().0.0;
        tracing::info!("old_target_end: {end} -> {new_end}");
        end = new_end;
    }

    (start < end).then_some((start, end))
}

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
    let ori_seq = if is_rev { target_rev } else { target_fwd };

    // 端点顺序颠倒或越过序列末端时返回空串，由调用方按“重比对失败”跳过。
    // is_rev 分支做减法，越界会 panic，所以必须先检查。
    if ori_start > ori_end || ori_end > ori_seq.len() {
        tracing::warn!(
            "invalid target range: {ori_start}-{ori_end}, target len {}",
            ori_seq.len()
        );
        return &ori_seq[0..0];
    }

    let (mut start, mut end) = get_target_start_end(ori_start, ori_end, target_fwd.len(), is_rev);

    let ori_seq_bytes = ori_seq.as_bytes();

    // 如果 slice 边界落在了 homopolymer 内部，就持续收缩
    if end < ori_seq.len() {
        while end > start && end >= 2 {
            if ori_seq_bytes[end - 2] == ori_seq_bytes[end - 1] {
                end -= 1;
            } else {
                break;
            }
        }
    }

    while start < end {
        if ori_seq_bytes[start] == ori_seq_bytes[start + 1] {
            start += 1
        } else {
            break;
        }
    }

    &ori_seq[start..end]
}

pub fn do_align_4_homo(ori_query_seq: &str, target_substr: &str) -> Option<Mapping> {
    if ori_query_seq.is_empty() || target_substr.is_empty() {
        return None;
    }

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
            ori_query_seq.as_bytes(),
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
    } else if !hits[0].as_ref().unwrap().is_primary {
        None
    } else {
        hits[0].take()
    }
}

#[cfg(test)]
mod test {
    use super::*;

    /// 一个 115 bp 的 poly-A 区域：区间 [308, 423) 覆盖 308..423 的每个位点。
    fn poly_a_region() -> RefPos2Motifs {
        let mut map = RefPos2Motifs::new();
        for pos in 308..423 {
            map.insert(pos, vec![((308, 423), Arc::new("(A)115".to_string()))]);
        }
        map
    }

    /// 两个端点都落在同一个 motif 区域内时，收缩后的窗口为空，必须返回 None，
    /// 不能返回 start > end 的区间（旧代码在这里会 panic）。
    #[test]
    fn test_shrink_same_region_returns_none() {
        let map = poly_a_region();
        assert_eq!(shrink_to_motif_boundary(350, 400, &map), None);
        assert_eq!(shrink_to_motif_boundary(308, 423, &map), None);
    }

    /// 起点在区域内、终点在区域外：起点后移到区域终点，窗口仍有效。
    #[test]
    fn test_shrink_start_inside_region() {
        let map = poly_a_region();
        assert_eq!(shrink_to_motif_boundary(350, 500, &map), Some((423, 500)));
    }

    /// 端点都不在区域内时保持原区间。
    #[test]
    fn test_shrink_no_region() {
        let map = poly_a_region();
        assert_eq!(shrink_to_motif_boundary(100, 200, &map), Some((100, 200)));
    }

    /// 端点颠倒或越界时返回空串，不 panic。
    #[test]
    fn test_get_target_substr_invalid_range() {
        let fwd = "ACGTACGT";
        let rev = "ACGTACGT";
        assert_eq!(get_target_substr(6, 2, false, fwd, rev), "");
        assert_eq!(get_target_substr(6, 2, true, fwd, rev), "");
        assert_eq!(get_target_substr(0, 100, false, fwd, rev), "");
    }
}
