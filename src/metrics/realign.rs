//! 三个 metric mode 共用的「重对齐前缀」流水线与 motif 正则表。
//!
//! hp-tr-v2 / hp-tr-v3 / hp-v2 的 `compute_metric` 对每个 hit 都执行同一段
//! 前缀流程：取 target 序列（hp-tr 两个 mode 先把比对端点收缩到 motif 区域
//! 边界）→ 切出 target/read 子串 → refAnchored 时两端反向互补 → 本地重对齐
//! → 丢弃空结果/仍反向的结果。差异只剩两点：是否做边界收缩、用哪张 motif
//! 正则表。本模块把这段流程收敛为 `realign_hit`，把正则表收敛为 `MotifSpec`。

use std::collections::HashMap;

use hp_tr_finder::UnitAndRepeats;
use regex::Regex;

use crate::{
    adapter::{Mapping, Strand},
    aligned_pairs::TAlignedPairs,
    global_data::{GlobalData, GlobalDataKey, GlobalDataValue},
    kits::ReadInfo,
    kits::reverse_complement,
    metrics::hp_tr_tools::{do_align_4_homo, get_target_substr, shrink_to_motif_boundary},
};

/// mode 的 motif 发现参数，派生 `single_seq_hp_tr_finder` 所需的正则表。
pub struct MotifSpec {
    /// motif unit 长度列表。HP 只认 1（同聚物），HP-TR 认 1..=4。
    pub unit_lens: &'static [u8],
    /// 最小重复次数。
    pub min_repeats: u8,
}

impl MotifSpec {
    pub fn build_regs(&self) -> Vec<HashMap<String, Regex>> {
        self.unit_lens
            .iter()
            .map(|&unit_len| UnitAndRepeats::new(unit_len, self.min_repeats).build_finder_regrex())
            .collect()
    }
}

pub const HP_MOTIF_SPEC: MotifSpec = MotifSpec {
    unit_lens: &[1],
    min_repeats: 3,
};

pub const HP_TR_MOTIF_SPEC: MotifSpec = MotifSpec {
    unit_lens: &[1, 2, 3, 4],
    min_repeats: 3,
};

lazy_static::lazy_static! {
    /// hp-v2 的 motif 正则表（unit 长度 1）。
    pub static ref HP_REG: Vec<HashMap<String, Regex>> = HP_MOTIF_SPEC.build_regs();
    /// hp-tr-v2/v3 的 motif 正则表（unit 长度 1..=4），GlobalData 的
    /// TargetRegion2Motif4HpTr 索引也由它构建。
    pub static ref HP_TR_REG: Vec<HashMap<String, Regex>> = HP_TR_MOTIF_SPEC.build_regs();
}

/// 重对齐前缀的 mode 差异配置。
pub struct RealignConfig {
    /// hp-tr-v2/v3 = true：比对端点先收缩到 HP-TR motif 区域边界；
    /// hp-v2 = false：不触碰 TargetRegion2Motif4HpTr 索引（保持其不初始化的现状）。
    pub shrink_to_boundary: bool,
    /// refAnchored：reverse hit 的 read/target 子串均反向互补后再重对齐。
    pub reference_anchored: bool,
}

/// 重对齐成功后的产物，三个 mode 的专属统计都只消费这三个字段。
pub struct RealignedHit {
    /// 已按 refAnchored 归一（reverse + refAnchored 时为反向互补序列）。
    pub read_seq: String,
    /// 已做 homopolymer 边界收缩与 refAnchored 归一。
    pub target_substr: String,
    /// `do_align_4_homo` 的结果，保证 strand == Forward。
    pub align: Mapping,
}

/// 单个 hit 的重对齐前缀流水线。
///
/// 返回 `None` 表示该 hit 应跳过（target 端点为 0 / 收缩成空 / 子串切片非法 /
/// 重对齐无结果 / 仍反向），并已打过对应的 warn。
pub fn realign_hit(
    read_info: &ReadInfo,
    mapping: &Mapping,
    global_data: &GlobalData,
    cfg: &RealignConfig,
) -> Option<RealignedHit> {
    let target_name = mapping.target_name.as_ref().unwrap();

    let global_data_value = global_data.get(GlobalDataKey::TargetName2SeqAndRev);
    let [target_seq_fwd, target_seq_rev] = match global_data_value {
        GlobalDataValue::TargetName2SeqAndRev(v) => v.get().unwrap().get(target_name).unwrap(),
        _ => panic!("TargetName2SeqAndRev not found"),
    };

    // hp-tr 两个 mode：如果比对的 reference_start, reference_end 正好处于
    // hp-str 区域的中心位置，就从两头往中间缩。两个端点都落在同一个区域内时
    // 窗口会被缩成空，此时无法确定完整重复数，跳过该 hit。
    let shrunk = if cfg.shrink_to_boundary {
        let global_data_value = global_data.get(GlobalDataKey::TargetRegion2Motif4HpTr);
        let region2motif = match global_data_value {
            GlobalDataValue::TargetRegion2Motif4HpTr(v) => {
                v.get().unwrap().get(target_name).unwrap()
            }
            _ => panic!("TargetRegion2Motif4HpTr not found"),
        };

        let raw_target_start = mapping.target_start as usize;
        let raw_target_end = mapping.target_end as usize;
        if raw_target_end == 0 {
            tracing::warn!("aligned target end is zero.");
            return None;
        }
        let Some((start, end)) =
            shrink_to_motif_boundary(raw_target_start, raw_target_end, region2motif)
        else {
            tracing::warn!(
                "target window collapses after motif-boundary shrink. QueryName:{}. TargetStartEnd:{}-{}, strand:{:?}",
                read_info.name,
                raw_target_start,
                raw_target_end,
                mapping.strand
            );
            return None;
        };
        Some((start, end))
    } else {
        None
    };

    let (tstart, tend) =
        shrunk.unwrap_or((mapping.target_start as usize, mapping.target_end as usize));

    let target_substr = get_target_substr(
        tstart,
        tend,
        mapping.is_reverse(),
        target_seq_fwd,
        target_seq_rev,
    );

    let read_seq = &read_info.seq[mapping.query_start as usize..mapping.query_end as usize];

    let (read_seq, target_substr) = if mapping.is_reverse() && cfg.reference_anchored {
        let read_rev = reverse_complement(read_seq.as_bytes());
        let read_rev = String::from_utf8(read_rev).unwrap();
        let target_rev = reverse_complement(target_substr.as_bytes());
        let target_rev = String::from_utf8(target_rev).unwrap();
        (read_rev, target_rev)
    } else {
        (read_seq.to_string(), target_substr.to_string())
    };

    let align = do_align_4_homo(&read_seq, &target_substr);
    let Some(align) = align else {
        tracing::warn!(
            "no aligned result. QueryName:{}. QueryStartEnd:{}-{}, TargetStartEnd:{}-{}, strand:{:?}",
            read_info.name,
            mapping.query_start,
            mapping.query_end,
            mapping.target_start,
            mapping.target_end,
            mapping.strand
        );
        return None;
    };

    if !matches!(align.strand, Strand::Forward) {
        tracing::warn!(
            "still reverse aligment. QueryName:{}. QueryStartEnd:{}-{}",
            read_info.name,
            mapping.query_start,
            mapping.query_end,
        );
        return None;
    }

    Some(RealignedHit {
        read_seq,
        target_substr,
        align,
    })
}
