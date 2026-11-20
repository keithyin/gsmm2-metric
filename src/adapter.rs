//! bandmap-core 与被移除的 mm2 之间的适配层。
//!
//! 职责：
//! 1. 定义项目内自有 `Mapping` 结构（字段与 mm2 `minimap2::Mapping` 同名同类型，
//!    可构造可测试），并从 bandmap 的比对结果转换而来；
//! 2. 复刻 mm2 `MappingExt` 的 `identity()`（EQX cigar op7 计数 / cigar 总长，
//!    分母含全部 op、不含 clip）与 `aligned_2_str()`；
//! 3. 提供单条 query 对全部 target 的比对入口（bandmap 内部已过滤 unmapped，
//!    并按 minimap2 语义标注 primary/supplementary）。

use std::sync::Arc;

use bandmap::api::CigarOp;
use gskits::dna::reverse_complement;
use gskits::ds::ReadInfo;

pub use bandmap::Strand;

/// 项目内自有 Mapping：由 bandmap 的比对结果转换而来。
/// 字段语义与 mm2 `minimap2::Mapping` 一致：
/// - 坐标为 i32，query_start/query_end 为原始 read 正向坐标；
/// - `cigar` 为 EQX 模式的 cigar（op 编码与 minimap2 一致：1=I 2=D 7='=' 8='X'），
///   不含 clip；None 表示 unmapped（bandmap 已在管线内过滤，正常不会出现）。
#[derive(Debug, Clone)]
pub struct Mapping {
    pub query_name: Option<Arc<String>>,
    pub query_start: i32,
    pub query_end: i32,
    pub strand: Strand,
    pub target_name: Option<Arc<String>>,
    pub target_len: i32,
    pub target_start: i32,
    pub target_end: i32,
    pub mapq: u32,
    pub is_primary: bool,
    pub is_supplementary: bool,
    pub cigar: Option<Vec<(u32, u8)>>,
    pub alignment_score: Option<i32>,
}

impl Mapping {
    pub fn is_rev(&self) -> bool {
        matches!(self.strand, Strand::Reverse)
    }

    /// 复刻 mm2 `MappingExt::identity`：EQX cigar 中 op7 长度 / 全部 op 长度之和。
    /// 未对齐（cigar None）返回 0.0。
    pub fn identity(&self) -> f32 {
        match &self.cigar {
            Some(cigar) => {
                let mut matched = 0u32;
                let mut tot = 0u32;
                for &(cnt, op) in cigar {
                    assert!(op != 0, "eqx required");
                    tot += cnt;
                    if op == 7 {
                        matched += cnt;
                    }
                }
                matched as f32 / tot as f32
            }
            None => 0.0,
        }
    }

    /// 复刻 mm2 `MappingExt::aligned_2_str`：target/query 分别取自身区间切片，
    /// reverse 时 query 切片整体反向互补，然后按 ref-forward 顺序消费 cigar，
    /// I 在 target 侧补 '-'，D 在 query 侧补 '-'。
    pub fn aligned_2_str(&self, target: &[u8], query: &[u8]) -> (String, String) {
        let sbr_target = &target[self.target_start as usize..self.target_end as usize];
        let sbr_target = String::from_utf8(sbr_target.to_vec()).unwrap();
        let sbr_query = &query[self.query_start as usize..self.query_end as usize];
        let sbr_query = if !self.is_rev() {
            String::from_utf8(sbr_query.to_vec()).unwrap()
        } else {
            String::from_utf8(reverse_complement(sbr_query)).unwrap()
        };

        let mut aligned_target = String::new();
        let mut aligned_query = String::new();

        let mut target_start = 0usize;
        let mut query_start = 0usize;

        self.cigar.as_ref().unwrap().iter().for_each(|&(cnt, code)| {
            let cnt = cnt as usize;
            match code {
                0 | 7 | 8 => {
                    aligned_target.push_str(&sbr_target[target_start..target_start + cnt]);
                    aligned_query.push_str(&sbr_query[query_start..query_start + cnt]);
                    target_start += cnt;
                    query_start += cnt;
                }
                1 => {
                    let ins = vec!["-"; cnt].join("");
                    aligned_target.push_str(&ins);
                    aligned_query.push_str(&sbr_query[query_start..query_start + cnt]);
                    query_start += cnt;
                }
                2 => {
                    let del = vec!["-"; cnt].join("");
                    aligned_target.push_str(&sbr_target[target_start..target_start + cnt]);
                    aligned_query.push_str(&del);
                    target_start += cnt;
                }
                _ => {
                    panic!("invalid code. {}", code);
                }
            }
        });

        (aligned_target, aligned_query)
    }
}

/// bandmap `Mapping` -> 项目内 `Mapping`。
///
/// bandmap 的 `cigar_ops`（eqx=true 时含 op7/8）不含 clip，与 mm2 的
/// `Alignment.cigar` 语义一致，直接复制。query_name 由调用方补充
/// （bandmap map_seq 不回传 qname）。
pub fn from_bandmap(m: &bandmap::Mapping, query_name: Option<Arc<String>>) -> Mapping {
    let cigar: Option<Vec<(u32, u8)>> = m
        .cigar_ops
        .as_ref()
        .map(|ops: &Vec<CigarOp>| ops.iter().map(|op| (op.len, op.op)).collect());

    Mapping {
        query_name,
        query_start: m.query_start as i32,
        query_end: m.query_end as i32,
        strand: m.strand,
        target_name: Some(Arc::new(m.target_name.to_string())),
        target_len: m.target_len as i32,
        target_start: m.target_start as i32,
        target_end: m.target_end as i32,
        mapq: m.mapq as u32,
        is_primary: m.is_primary,
        is_supplementary: m.is_supplementary,
        cigar,
        alignment_score: Some(m.score),
    }
}

/// 复刻 mm2 `set_primary_supp_alignment_according_2_align_score` 的 hit 重排语义。
///
/// 两引擎标注体系的对应关系：
/// - bandmap `is_primary && !is_supplementary`（sam_pri，恒有一个）→ mm2 的引擎 primary，
///   参与按 score 的 primary 竞争（+1e9 加权）；
/// - bandmap `is_primary && is_supplementary`（2nd+ 非 secondary）→ 对应 mm2 里
///   minimap2 已判为 supplementary/split 的 hit（parent==id 且 sam_pri==0，is_primary=false），
///   **不参与** primary 竞争，只接受 supplementary 去重（与已收 supp/primary 的
///   query 重叠 > 0.5 则降级 secondary）；
/// - bandmap secondary → 保持 secondary。
///
/// ovlp_ratio 与 mm2 相同：重叠长度 / min(两侧长度)。
fn reorder_like_mm2(hits: &mut Vec<Mapping>) {
    if hits.is_empty() {
        return;
    }

    fn ovlp_ratio(s1: i32, e1: i32, s2: i32, e2: i32) -> f32 {
        if s2 >= e1 || s1 >= e2 {
            return 0.0;
        }
        let s = std::cmp::max(s1, s2);
        let e = std::cmp::min(e1, e2);
        let ovlp_len = e - s;
        let min_len = std::cmp::min(e1 - s1, e2 - s2);
        if min_len == 0 {
            return 0.0;
        }
        ovlp_len as f32 / min_len as f32
    }

    // 旗标归一：bandmap 的 2nd+ 非 secondary（primary=true && supp=true）对应
    // minimap2 的 split/supplementary（parent==id && sam_pri==0 → is_primary=false）。
    // 归一后 is_primary=true 只剩 sam_pri 型（恒一个），与 mm2 的引擎旗标语义一致。
    // 不做归一会导致第 4 步 supp 去重降级时留下 is_primary=true 的「假 primary」。
    for h in hits.iter_mut() {
        if h.is_primary && h.is_supplementary {
            h.is_primary = false;
        }
    }

    hits.sort_by_key(|v| {
        -v.alignment_score.unwrap()
            - if v.is_primary && !v.is_supplementary {
                1_000_000_000
            } else {
                0
            }
    });

    let primary_hit = &hits[0];
    let (primary_qstart, primary_qend) = (primary_hit.query_start, primary_hit.query_end);

    // 多 target 场景下可能存在的其余 sam_pri 型 engine-primary：按 ovlp 判定
    // supplementary 或降级 secondary（bandmap 单索引下通常为空，保留 mm2 逻辑）。
    for hit in hits.iter_mut().skip(1) {
        if hit.is_primary {
            let ratio = ovlp_ratio(primary_qstart, primary_qend, hit.query_start, hit.query_end);
            hit.is_primary = false;
            hit.is_supplementary = ratio <= 0.5;
        }
    }

    hits.sort_by_key(|v| {
        -v.alignment_score.unwrap()
            - if v.is_primary && !v.is_supplementary {
                1_000_000_000
            } else if v.is_supplementary {
                100_000_000
            } else {
                0
            }
    });

    let mut visited_ranges = vec![(primary_qstart, primary_qend)];
    for v in hits.iter_mut().skip(1) {
        if v.is_supplementary {
            let (start, end) = (v.query_start, v.query_end);
            let cnt = visited_ranges
                .iter()
                .map(|&(range_s, range_e)| ovlp_ratio(range_s, range_e, start, end))
                .filter(|&ovlp| ovlp > 0.5)
                .count();
            if cnt > 0 {
                v.is_supplementary = false;
            } else {
                visited_ranges.push((start, end));
            }
        }
    }
}

/// 单条 query 对全部 target 的比对，等价于 mm2 `align_single_query_to_targets`。
///
/// bandmap 单索引一次性完成多 target 比对，unmapped 已在管线内过滤；
/// 返回顺序为 mm2 语义：primary → supplementary（按分）→ secondary（按分）。
pub fn align_single_query_to_targets(aligner: &bandmap::Aligner, read: &ReadInfo) -> Vec<Mapping> {
    let qname = Some(Arc::new(read.name.clone()));
    let result = aligner.map_seq(&read.name, read.seq.as_bytes());
    let mut hits: Vec<Mapping> = result
        .mappings
        .iter()
        .map(|m| from_bandmap(m, qname.clone()))
        .collect();
    reorder_like_mm2(&mut hits);
    hits
}
