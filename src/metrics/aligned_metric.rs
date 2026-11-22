//! 全比对质量 metric（mm2 `gsmm2-aligned-metric` bin 的移植）。
//!
//! 与 hp-tr / hp 系 mode 的根本差异：直接消费 **raw mapping**（不做重对齐、
//! 不做边界收缩、不做 identity 过滤——它本身就是评估比对质量的工具），
//! 对一条 read 的全部 primary + supplementary hit 做多片段合并统计，
//! **每条输入 read 恰好输出一行**（unmapped / `--noMar` 过滤的 read 输出
//! 全 0 占位行，identity 列为 `NaN`，与 mm2 版逐字节对齐）。
//!
//! 移植时绕开了 mm2 的 `TseqAndRecord`/rust-htslib record 机制：metric 路径的
//! record cigar 恒不含 soft clip（query 跨度 == query_len），其本质只是
//! 「(reverse 时反序的) cigar 逐 op 展开」，由本模块的 `alignment_columns`
//! 直接完成，列序与 mm2 `TseqAndRecord::alignment_pair` 逐字节一致：
//! - query 侧恒为 read 的**原始方向**切片（per-base 桶按 read 原始方向碱基计数）；
//! - reverse hit 的 target 侧取反向互补切片（read 方向），cigar 反序后按
//!   read 方向消费。

use crate::adapter::Mapping;
use crate::kits::dna::SEQ_NT4_TABLE;
use crate::kits::regions::Regions;
use crate::metrics::hp_tr_tools::get_target_start_end;
use crate::metrics::{Metric, MetricContext, join_row};

#[derive(Debug, Default)]
pub struct AlignedMetric {
    qlen: usize,

    matched_cnt: [usize; 4],
    mismatched_cnt: [usize; 4],
    homodel_cnt: [usize; 4],
    non_homodel_cnt: [usize; 4],
    homoins_cnt: [usize; 4],
    non_homoins_cnt: [usize; 4],

    num_segs: usize,
    q_ovlp: usize,
    q_ovlp_ratio: f32,
    r_ovlp_ratio: f32,
    covlen: usize,
    primary_covlen: usize,
    ori_align_info: String,
    ori_q_gaps: String,
    merged_qry_span: String,

    long_indels: Vec<i32>,
}

impl AlignedMetric {
    pub fn matched(&self) -> usize {
        self.matched_cnt.iter().copied().sum::<usize>()
    }

    pub fn mismatched(&self) -> usize {
        self.mismatched_cnt.iter().copied().sum::<usize>()
    }

    pub fn homoins(&self) -> usize {
        self.homoins_cnt.iter().copied().sum::<usize>()
    }

    pub fn non_homoins(&self) -> usize {
        self.non_homoins_cnt.iter().copied().sum::<usize>()
    }

    pub fn homodel(&self) -> usize {
        self.homodel_cnt.iter().copied().sum::<usize>()
    }

    pub fn non_homodel(&self) -> usize {
        self.non_homodel_cnt.iter().copied().sum::<usize>()
    }

    pub fn aligned_span(&self) -> usize {
        self.matched()
            + self.mismatched()
            + self.homoins()
            + self.non_homoins()
            + self.homodel()
            + self.non_homodel()
    }

    /// 行构造器：列顺序/格式与 mm2 `align_metric::Metric::Display` 逐列一致。
    /// 占位行（全零状态）也走同一构造器：identity = 0/0 → `NaN`、
    /// queryCoverage = 0/qlen，与 mm2 对 unmapped read 的输出字节对齐。
    fn build_row(&self, qname: &str, tname: &str) -> String {
        let mut cells = Vec::with_capacity(Self::COLUMNS.len());
        cells.push(qname.to_string());
        cells.push(tname.to_string());
        cells.push(self.qlen.to_string());
        cells.push(self.num_segs.to_string());
        cells.push(self.covlen.to_string());
        cells.push(self.primary_covlen.to_string());
        cells.push(format!("{:.6}", self.covlen as f64 / self.qlen as f64));
        cells.push(format!(
            "{:.6}",
            self.matched() as f64 / self.aligned_span() as f64
        ));
        cells.push(self.ori_align_info.clone());
        cells.push(self.ori_q_gaps.clone());
        cells.push(self.q_ovlp.to_string());
        cells.push(format!("{:.6}", self.q_ovlp_ratio));
        cells.push(format!("{:.6}", self.r_ovlp_ratio));
        cells.push(self.merged_qry_span.clone());
        cells.push(self.matched().to_string());
        cells.push(self.mismatched().to_string());
        cells.push(self.non_homoins().to_string());
        cells.push(self.homoins().to_string());
        cells.push(self.non_homodel().to_string());
        cells.push(self.homodel().to_string());
        cells.push(
            self.long_indels
                .iter()
                .map(|v| v.to_string())
                .collect::<Vec<_>>()
                .join(","),
        );
        for idx in 0..4 {
            cells.push(self.matched_cnt[idx].to_string());
            cells.push(self.mismatched_cnt[idx].to_string());
            cells.push(self.non_homoins_cnt[idx].to_string());
            cells.push(self.homoins_cnt[idx].to_string());
            cells.push(self.non_homodel_cnt[idx].to_string());
            cells.push(self.homodel_cnt[idx].to_string());
        }
        debug_assert_eq!(cells.len(), Self::COLUMNS.len());
        join_row(&cells)
    }
}

/// 列名与 mm2 `align_metric::METRIC_CSV_HEADER` 完全一致（45 列）。
impl Metric for AlignedMetric {
    const COLUMNS: &'static [&'static str] = &[
        "qname",
        "rname",
        "qlen",
        "segs",
        "covlen",
        "primaryCovlen",
        "queryCoverage",
        "identity",
        "oriAlignInfo",
        "oriQGaps",
        "qOvlp",
        "qOvlpRatio",
        "rOvlpRatio",
        "mergedQrySpan",
        "match",
        "misMatch",
        "ins",
        "homoIns",
        "del",
        "homoDel",
        "longIndel",
        "match-A",
        "misMatch-A",
        "ins-A",
        "homoIns-A",
        "del-A",
        "homoDel-A",
        "match-C",
        "misMatch-C",
        "ins-C",
        "homoIns-C",
        "del-C",
        "homoDel-C",
        "match-G",
        "misMatch-G",
        "ins-G",
        "homoIns-G",
        "del-G",
        "homoDel-G",
        "match-T",
        "misMatch-T",
        "ins-T",
        "homoIns-T",
        "del-T",
        "homoDel-T",
    ];

    /// 质量评估类 metric：不做 identity 过滤（对齐 mm2 原版语义），
    /// 低 identity read 恰恰是它要度量的对象。
    const IDENTITY_THRESHOLD: Option<f32> = None;

    fn compute(&mut self, ctx: MetricContext) -> Vec<String> {
        self.qlen = ctx.read.seq.len();

        // unmapped / --noMar 过滤的 read：占位行（mm2 语义：每条输入 read 恰好一行）
        if ctx.mappings.is_empty() {
            return vec![self.build_row(&ctx.read.name, &ctx.target_name)];
        }

        let seqs = ctx.global_data.seqs_with_rev();
        let [target_fwd, target_rev] = seqs.get(&ctx.target_name).unwrap();
        let (target_fwd, target_rev) = (target_fwd.as_str(), target_rev.as_str());

        // primary_covlen：最后一个 primary 的 query 跨度（mm2 循环覆盖语义）
        for v in &ctx.mappings {
            if v.is_primary {
                self.primary_covlen = (v.query_end - v.query_start) as usize;
            }
        }

        // mm2：按 qstart 排序后做区间运算与去重叠截断
        let mut sorted: Vec<&Mapping> = ctx.mappings.iter().collect();
        sorted.sort_by_key(|v| v.query_start);

        let qstart_ends = sorted
            .iter()
            .map(|v| (v.query_start as usize, v.query_end as usize))
            .collect::<Vec<(usize, usize)>>();

        let rstart_ends = sorted
            .iter()
            .map(|v| {
                let s = v.target_start as usize;
                let e = v.target_end as usize;
                (s.min(e), e.max(s))
            })
            .collect::<Vec<(usize, usize)>>();

        let qstart_ends_region: Regions = (&qstart_ends).into();
        self.q_ovlp_ratio = qstart_ends_region.ovlp_ratio();

        self.ori_q_gaps = qstart_ends_region
            .gaps(Some(0), Some(self.qlen))
            .iter()
            .map(|v| v.to_string())
            .collect::<Vec<_>>()
            .join(",");

        let rstart_ends_region: Regions = (&rstart_ends).into();
        self.r_ovlp_ratio = rstart_ends_region.ovlp_ratio();

        // 去重叠截断（照抄 mm2 align_metric.rs：后者与前段重叠时裁短先到者）
        let mut truncated_qstarts_ends = vec![qstart_ends[0]];
        let mut ovlp_len = 0;
        qstart_ends
            .iter()
            .skip(1)
            .for_each(|&(mut cur_start, cur_end)| {
                let (pre_start, pre_end) = truncated_qstarts_ends.last_mut().unwrap();
                if cur_start >= *pre_end {
                    truncated_qstarts_ends.push((cur_start, cur_end));
                } else {
                    ovlp_len += *pre_end - cur_start;

                    if (cur_end - cur_start) > (*pre_end - *pre_start) {
                        *pre_end = cur_start;
                    } else {
                        cur_start = *pre_end;
                    }

                    truncated_qstarts_ends.push((cur_start, cur_end));
                }
            });

        self.num_segs = truncated_qstarts_ends.len();
        let coverlen = truncated_qstarts_ends
            .iter()
            .map(|&(s, e)| e - s)
            .sum::<usize>();
        self.covlen = coverlen;

        sorted
            .iter()
            .zip(truncated_qstarts_ends.iter())
            .for_each(|(mapping, se)| {
                let (aligned_ref, aligned_qry) =
                    alignment_columns(mapping, target_fwd, target_rev, &ctx.read.seq, se.0, se.1);
                fill_metric_from_align_pair(aligned_ref.as_bytes(), aligned_qry.as_bytes(), self);
            });

        // 消费 query 的列数（match/mismatch/ins）必须恰好等于截断后的总跨度
        let cov2 = self.matched_cnt.iter().copied().sum::<usize>()
            + self.mismatched_cnt.iter().copied().sum::<usize>()
            + self.homoins_cnt.iter().copied().sum::<usize>()
            + self.non_homoins_cnt.iter().copied().sum::<usize>();

        assert_eq!(cov2, coverlen);

        self.q_ovlp = ovlp_len;

        let mut ori_align_info = ctx
            .mappings
            .iter()
            .map(|v| (v.query_start, align_info_str(v)))
            .collect::<Vec<_>>();
        ori_align_info.sort_by_key(|v| v.0);
        self.ori_align_info = ori_align_info
            .iter()
            .map(|v| v.1.as_str())
            .collect::<Vec<_>>()
            .join(";");

        self.merged_qry_span = truncated_qstarts_ends
            .iter()
            .map(|v| format!("{}:{}", v.0, v.1))
            .collect::<Vec<_>>()
            .join(";");

        vec![self.build_row(&ctx.read.name, &ctx.target_name)]
    }
}

/// oriAlignInfo 列的单 hit 片段，与 mm2 `AlignInfo::Display` 一致：
/// `{qstart}:{qend}->{rstart}:{rend}{+/-}`，target 坐标用原始（ref-forward）
/// 值，不因 strand 换序。
fn align_info_str(v: &Mapping) -> String {
    format!(
        "{}:{}->{}:{}{}",
        v.query_start,
        v.query_end,
        v.target_start,
        v.target_end,
        if v.is_rev() { "-" } else { "+" }
    )
}

/// 为单个 hit 生成 [q_lo, q_hi)（绝对 query 坐标，read 原始方向）内的比对列串，
/// 返回 (aligned_ref, aligned_qry)，'-' 表示该侧缺口。
///
/// 等价于 mm2 `TseqAndRecord::alignment_pair`（经 `fwd_record` 构造 clip-0 record
/// 后的逐列重放），列序与之一致。走查终止/包含判定也照抄：
/// - del 列（q 缺口）不推进 q 游标，包含判定沿用当前游标值；
/// - `q 游标 >= q_hi` 即终止整条走查（`r 游标 >= rend` 在本构造下不可达，
///   保留以对齐 mm2 的终止条件）。
fn alignment_columns(
    mapping: &Mapping,
    target_fwd: &str,
    target_rev: &str,
    read_seq: &str,
    q_lo: usize,
    q_hi: usize,
) -> (String, String) {
    let ts = mapping.target_start as usize;
    let te = mapping.target_end as usize;
    let qs = mapping.query_start as usize;
    let qe = mapping.query_end as usize;

    // target 子串：reverse 时取全 target 反向互补序列的对应切片，
    // 等价于 RC(target[ts..te])（与 mm2 fwd_record 的 aligned_target_seq 一致）
    let (tseq, qseq) = if mapping.is_rev() {
        let (rs, re) = get_target_start_end(ts, te, target_fwd.len(), true);
        (&target_rev.as_bytes()[rs..re], &read_seq.as_bytes()[qs..qe])
    } else {
        (&target_fwd.as_bytes()[ts..te], &read_seq.as_bytes()[qs..qe])
    };

    let rend = te - ts; // record 局部 ref 跨度
    let qstart_rel = q_lo - qs; // record 局部 query 区间（截断后）
    let qend_rel = q_hi - qs;

    let cigar = mapping.cigar.as_ref().unwrap();
    let cigar_iter: Box<dyn Iterator<Item = &(u32, u8)>> = if mapping.is_rev() {
        Box::new(cigar.iter().rev())
    } else {
        Box::new(cigar.iter())
    };

    let mut q_cursor: Option<usize> = None;
    let mut r_cursor: Option<usize> = None;
    let mut aligned_ref = String::new();
    let mut aligned_qry = String::new();

    // 逐碱基游标（record 局部坐标：q 相对 hit 的 qstart，r 相对 hit 的 target start）
    let mut q_base = 0usize;
    let mut r_base = 0usize;

    'outer: for &(cnt, op) in cigar_iter {
        for _ in 0..cnt {
            let (qpos, rpos) = match op {
                // op 0（通用 M）与 7/8 同为双推进（EQX 管线中不会出现 0，
                // 与 adapter::aligned_2_str 的口径一致）
                0 | 7 | 8 => {
                    let q = q_base;
                    let r = r_base;
                    q_base += 1;
                    r_base += 1;
                    (Some(q), Some(r))
                }
                1 => {
                    let q = q_base;
                    q_base += 1;
                    (Some(q), None)
                }
                2 => {
                    let r = r_base;
                    r_base += 1;
                    (None, Some(r))
                }
                v => panic!("invalid cigar op: {}", v),
            };

            if rpos.is_some() {
                r_cursor = rpos;
            }
            if qpos.is_some() {
                q_cursor = qpos;
            }

            if r_cursor.unwrap_or(0) >= rend || q_cursor.unwrap_or(0) >= qend_rel {
                break 'outer;
            }

            if q_cursor.unwrap_or(0) >= qstart_rel {
                if let Some(qpos) = qpos {
                    aligned_qry.push(qseq[qpos] as char);
                } else {
                    aligned_qry.push('-');
                }

                if let Some(rpos) = rpos {
                    aligned_ref.push(tseq[rpos] as char);
                } else {
                    aligned_ref.push('-');
                }
            }
        }
    }

    (aligned_ref, aligned_qry)
}

/// 逐列扫描比对串，填充 match/mismatch/homo ins/del 计数
/// （mm2 `fill_metric_from_align_pair` 逐行移植）。
pub fn fill_metric_from_align_pair(
    aligned_ref: &[u8],
    aligned_qry: &[u8],
    metric: &mut AlignedMetric,
) {
    assert_eq!(aligned_ref.len(), aligned_qry.len());

    let mut ins_bases_idx = vec![];
    let mut del_bases_idx = vec![];

    for idx in 0..aligned_qry.len() {
        let qbase = aligned_qry[idx];
        let rbase = aligned_ref[idx];

        match (qbase as char, rbase as char) {
            ('-', '-') => panic!("both ref and qry are gap at column {}", idx),

            ('-', _rbase) => {
                // deletion
                del_bases_idx.push(idx);
                fill_insertion_info(&ins_bases_idx, aligned_ref, aligned_qry, metric);
                ins_bases_idx.clear();
            }

            (_qbase, '-') => {
                // insertion
                ins_bases_idx.push(idx);
                fill_deletion_info_v2(&del_bases_idx, aligned_ref, metric);
                del_bases_idx.clear();
            }

            (qbase, rbase) => {
                if rbase == qbase {
                    metric.matched_cnt[SEQ_NT4_TABLE[qbase as usize] as usize] += 1;
                } else {
                    metric.mismatched_cnt[SEQ_NT4_TABLE[qbase as usize] as usize] += 1;
                }

                fill_deletion_info_v2(&del_bases_idx, aligned_ref, metric);
                fill_insertion_info(&ins_bases_idx, aligned_ref, aligned_qry, metric);

                del_bases_idx.clear();
                ins_bases_idx.clear();
            }
        }
    }

    fill_deletion_info_v2(&del_bases_idx, aligned_ref, metric);
    fill_insertion_info(&ins_bases_idx, aligned_ref, aligned_qry, metric);
}

/// homo del 判定（mm2 v2 口径）：被删 ref 碱基与其左/右 ref 邻碱基任一一致即 homo。
/// 连续删除 >= 10 时额外记一条负值长 indel。
fn fill_deletion_info_v2(del_bases_idx: &[usize], aligned_ref: &[u8], metric: &mut AlignedMetric) {
    if del_bases_idx.len() >= 10 {
        metric.long_indels.push(-(del_bases_idx.len() as i32));
    }
    if !del_bases_idx.is_empty() {
        del_bases_idx.iter().for_each(|&base_idx| {
            let pre_ref_base = get_pre_rbase(aligned_ref, base_idx);
            let next_ref_base = get_next_rbase(aligned_ref, base_idx);
            let cur_base = aligned_ref[base_idx];
            let mut homo = false;
            if let Some(pre_ref_base) = pre_ref_base {
                homo |= cur_base == pre_ref_base;
            }

            if let Some(next_ref_base) = next_ref_base {
                homo |= cur_base == next_ref_base;
            }

            if homo {
                metric.homodel_cnt[SEQ_NT4_TABLE[cur_base as usize] as usize] += 1;
            } else {
                metric.non_homodel_cnt[SEQ_NT4_TABLE[cur_base as usize] as usize] += 1;
            }
        });
    }
}

/// homo ins 判定（mm2 口径）：连续插入碱基一致，且与左/右 ref 邻碱基其一一致。
/// 连续插入 >= 10 时额外记一条正值长 indel。
fn fill_insertion_info(
    ins_bases_idx: &[usize],
    aligned_ref: &[u8],
    aligned_qry: &[u8],
    metric: &mut AlignedMetric,
) {
    if ins_bases_idx.len() >= 10 {
        metric.long_indels.push(ins_bases_idx.len() as i32);
    }

    if !ins_bases_idx.is_empty() {
        let pre = get_pre_rbase(aligned_ref, ins_bases_idx.first().copied().unwrap());
        let next = get_next_rbase(aligned_ref, ins_bases_idx.last().copied().unwrap());
        let mut homo = false;
        let ins_bases = ins_bases_idx
            .iter()
            .map(|&idx| aligned_qry[idx])
            .collect::<Vec<_>>();
        if let Some(pre) = pre {
            homo |= ins_bases.iter().filter(|&&b| b == pre).count() == ins_bases.len();
        }
        if let Some(next) = next {
            homo |= ins_bases.iter().filter(|&&b| b == next).count() == ins_bases.len();
        }
        if homo {
            ins_bases
                .iter()
                .for_each(|&v| metric.homoins_cnt[SEQ_NT4_TABLE[v as usize] as usize] += 1);
        } else {
            ins_bases
                .iter()
                .for_each(|&v| metric.non_homoins_cnt[SEQ_NT4_TABLE[v as usize] as usize] += 1);
        }
    }
}

fn get_next_rbase(aligned_ref: &[u8], idx: usize) -> Option<u8> {
    let gap = b'-';

    let mut idx_cursor_for_next = idx;

    loop {
        if (idx_cursor_for_next + 1) >= aligned_ref.len() {
            break None;
        }

        idx_cursor_for_next += 1;
        if aligned_ref[idx_cursor_for_next] != gap {
            break Some(aligned_ref[idx_cursor_for_next]);
        }
    }
}

fn get_pre_rbase(aligned_ref: &[u8], idx: usize) -> Option<u8> {
    let gap = b'-';
    let mut idx_cursor_for_pre = idx;

    loop {
        if idx_cursor_for_pre == 0 {
            break None;
        }
        idx_cursor_for_pre -= 1;

        if aligned_ref[idx_cursor_for_pre] != gap {
            break Some(aligned_ref[idx_cursor_for_pre]);
        }
    }
}

#[cfg(test)]
mod test {
    use std::{collections::HashMap, sync::Arc};

    use super::{AlignedMetric, alignment_columns, fill_metric_from_align_pair};
    use crate::adapter::{self, Mapping, Strand};
    use crate::aligner::build_aligner_from_params;
    use crate::global_data::GlobalData;
    use crate::kits::{
        ReadInfo,
        dna::SEQ_NT4_TABLE,
        fastx::{FastaFileReader, read_fastx},
        reverse_complement,
    };
    use crate::metrics::{Metric, MetricContext};
    use crate::params::{AlignParams, IndexParams};

    /// mm2 test_fill_align_info_from_align_pair 的移植 + 精确断言。
    /// ref: C-GGG-G
    /// qry: GG---CA
    #[test]
    fn test_fill_metric_from_align_pair() {
        let aligned_ref = b"C-GGG-G";
        let aligned_qry = b"GG---CA";
        let mut metric = AlignedMetric::default();
        fill_metric_from_align_pair(aligned_ref, aligned_qry, &mut metric);

        assert_eq!(
            metric.matched(),
            0,
            "该比对串没有 match 列（C/G、A/G 均错配）"
        );
        assert_eq!(metric.mismatched(), 2, "idx0: G vs C；idx6: A vs G");
        assert_eq!(
            metric.mismatched_cnt[SEQ_NT4_TABLE[b'G' as usize] as usize],
            1
        );
        assert_eq!(
            metric.mismatched_cnt[SEQ_NT4_TABLE[b'A' as usize] as usize],
            1
        );

        // ins run [1]: q=G, pre-ref=C, next-ref=G → 与 next 一致 → homo
        assert_eq!(metric.homoins(), 1);
        assert_eq!(metric.homoins_cnt[SEQ_NT4_TABLE[b'G' as usize] as usize], 1);
        // ins run [5]: q=C, pre/next-ref=G → 不一致 → non-homo
        assert_eq!(metric.non_homoins(), 1);
        assert_eq!(
            metric.non_homoins_cnt[SEQ_NT4_TABLE[b'C' as usize] as usize],
            1
        );

        // del run [2..5]: 3 个 G，各自与左/右 ref 邻碱基（G）一致 → 全 homo
        assert_eq!(metric.homodel(), 3);
        assert_eq!(metric.homodel_cnt[SEQ_NT4_TABLE[b'G' as usize] as usize], 3);
        assert_eq!(metric.non_homodel(), 0);

        assert_eq!(
            metric.aligned_span(),
            7,
            "2 mismatch + 1 ins + 1 ins + 3 del"
        );
        assert!(metric.long_indels.is_empty());
    }

    /// 长删除（>=10）记录到 longIndel（del 记负值）。
    /// ref: C*11 + A*4，qry: 11 个删除 + 4 个 G（mismatch）。
    #[test]
    fn test_long_deletion() {
        let aligned_ref = b"CCCCCCCCCCCAAAA";
        let aligned_qry = b"-----------GGGG";
        let mut metric = AlignedMetric::default();
        fill_metric_from_align_pair(aligned_ref, aligned_qry, &mut metric);

        assert_eq!(metric.long_indels, vec![-11], "del run 11 → 记 -11");
        assert_eq!(metric.homodel(), 11, "11 个 C 与右邻 C 一致 → homo del");
        assert_eq!(metric.mismatched(), 4, "4 个 A vs G mismatch");
    }

    /// 长插入：ins run 10 → 记 +10。
    /// ref: A -×10 GGGG，qry: A C×10 GGGG。
    #[test]
    fn test_long_insertion() {
        let aligned_ref = b"A----------GGGG";
        let aligned_qry = b"ACCCCCCCCCCGGGG";
        let mut metric = AlignedMetric::default();
        fill_metric_from_align_pair(aligned_ref, aligned_qry, &mut metric);

        assert_eq!(metric.long_indels, vec![10], "ins run 10 → 记 +10");
        assert_eq!(
            metric.non_homoins(),
            10,
            "C 与 pre=A / next=G 均不一致 → non-homo"
        );
    }

    fn make_mapping(
        query_start: i32,
        query_end: i32,
        target_start: i32,
        target_end: i32,
        strand: Strand,
        cigar: Vec<(u32, u8)>,
    ) -> Mapping {
        Mapping {
            query_name: None,
            query_start,
            query_end,
            strand,
            target_name: Some(Arc::new("t".to_string())),
            target_len: 64,
            target_start,
            target_end,
            mapq: 60,
            is_primary: true,
            is_supplementary: false,
            cigar: Some(cigar),
            alignment_score: None,
        }
    }

    /// forward hit：列 = target-forward × read-original。
    /// target[2..10) = "AACCCCGG"，cigar = 4= 1I 4=，read = "AACC" + 插入 A + "CCGG"。
    #[test]
    fn test_alignment_columns_fwd() {
        let target_fwd = "AAAACCCCGGGGTTTTACGTAC";
        let target_rev = String::from_utf8(reverse_complement(target_fwd.as_bytes())).unwrap();
        let read = "AACCACCGG";
        let mapping = make_mapping(0, 9, 2, 10, Strand::Forward, vec![(4, 7), (1, 1), (4, 7)]);

        let (aligned_ref, aligned_qry) =
            alignment_columns(&mapping, target_fwd, &target_rev, read, 0, 9);

        assert_eq!(aligned_ref, "AACC-CCGG");
        assert_eq!(aligned_qry, "AACCACCGG");
        assert_eq!(aligned_ref.len(), aligned_qry.len());
    }

    /// reverse hit：target 侧取 RC 切片、cigar 反序、query 恒为原始方向——
    /// per-base 桶因此按 read 原始方向碱基计数。
    #[test]
    fn test_alignment_columns_rev() {
        // target_fwd[0..8) = "AAAACCCC"，其 RC = "GGGGTTTT"
        let target_fwd = "AAAACCCCGGGGTTTTACGTAC";
        let target_rev = String::from_utf8(reverse_complement(target_fwd.as_bytes())).unwrap();
        // read 即该片段的反向互补（与 RC(target[0..8)) 相同）→ 8 列全相等
        let read = "GGGGTTTT";
        let mapping = make_mapping(0, 8, 0, 8, Strand::Reverse, vec![(8, 7)]);

        let (aligned_ref, aligned_qry) =
            alignment_columns(&mapping, target_fwd, &target_rev, read, 0, 8);

        assert_eq!(
            aligned_ref, "GGGGTTTT",
            "target 侧应为 RC 后的 read 方向切片"
        );
        assert_eq!(aligned_qry, "GGGGTTTT", "query 侧恒为 read 原始方向");

        // 若方向处理错误（target 用正向切片），会得到全 mismatch
        let mut metric = AlignedMetric::default();
        fill_metric_from_align_pair(aligned_ref.as_bytes(), aligned_qry.as_bytes(), &mut metric);
        assert_eq!(metric.matched(), 8);
        assert_eq!(metric.mismatched(), 0);
        assert_eq!(metric.matched_cnt[SEQ_NT4_TABLE[b'G' as usize] as usize], 4);
        assert_eq!(metric.matched_cnt[SEQ_NT4_TABLE[b'T' as usize] as usize], 4);
    }

    /// reverse hit 截断：[q_lo, q_hi) 过滤（截断区间是 hit 区间的子集）。
    #[test]
    fn test_alignment_columns_rev_truncation() {
        let target_fwd = "AAAACCCCGGGGTTTTACGTAC";
        let target_rev = String::from_utf8(reverse_complement(target_fwd.as_bytes())).unwrap();
        let read = "GGGGTTTT";
        let mapping = make_mapping(0, 8, 0, 8, Strand::Reverse, vec![(8, 7)]);

        // 截取 query [2, 6) → 4 列
        let (aligned_ref, aligned_qry) =
            alignment_columns(&mapping, target_fwd, &target_rev, read, 2, 6);
        assert_eq!(aligned_ref, "GGTT");
        assert_eq!(aligned_qry, "GGTT");
    }

    fn default_metric_context(
        read: ReadInfo,
        target_name: Arc<String>,
        mappings: Vec<Mapping>,
        global_data: Arc<GlobalData>,
    ) -> MetricContext {
        MetricContext {
            read,
            target_name,
            mappings,
            global_data,
            reference_anchored: false,
        }
    }

    /// 占位行：unmapped（空 mappings）read 输出全 0 行，45 列、identity 为 NaN。
    #[test]
    fn test_placeholder_row() {
        let global_data = Arc::new(GlobalData::new(HashMap::new()));
        let read = ReadInfo::new_fa_record("q1".to_string(), "ACGTACGTAC".to_string());
        let mut metric = AlignedMetric::default();
        let rows = metric.compute(default_metric_context(
            read,
            Arc::new(String::new()),
            vec![],
            global_data,
        ));

        assert_eq!(rows.len(), 1, "每条输入 read 恰好一行");
        let cells: Vec<&str> = rows[0].split('\t').collect();
        assert_eq!(cells.len(), AlignedMetric::COLUMNS.len());
        assert_eq!(cells[0], "q1");
        assert_eq!(cells[1], "", "rname 为空");
        assert_eq!(cells[2], "10", "qlen");
        assert_eq!(cells[3], "0", "segs");
        assert_eq!(cells[6], "0.000000", "queryCoverage = 0/qlen");
        assert_eq!(cells[7], "NaN", "identity = 0/0");
        assert_eq!(cells[20], "", "longIndel 为空");
        assert_eq!(cells[14], "0", "match");
    }

    /// 去重叠截断逻辑（照抄 mm2）：[(0,10),(5,15),(20,30)] →
    /// [(0,10),(10,15),(20,30)]，ovlp=5，covlen=25。
    #[test]
    fn test_truncation() {
        // 通过公共路径验证：3 个不重叠-部分重叠的 hit 手工构造（跳过真实比对）。
        // 这里直接构造 mappings + 简单序列，走 AlignedMetric::compute。
        let target = "AAAACCCCGGGGTTTTACGTACGTAAAACCCCGGGGTTTTACGTACGT";
        let targetname2seq =
            HashMap::from([(Arc::new("t".to_string()), Arc::new(target.to_string()))]);
        let global_data = Arc::new(GlobalData::new(targetname2seq));

        // 三个 hit：q [0,10)→t[0,10)；q [5,15)→t[16,26)；q [20,30)→t[32,42)
        // （第二/三个的 query 与前序重叠，触发截断逻辑）
        let read = "AAAACCCCGGGGTTTTACGTAAAACCCCGGGGTTTTACGTACGTAAAACC";
        let mappings = vec![
            make_mapping(0, 10, 0, 10, Strand::Forward, vec![(10, 7)]),
            make_mapping(5, 15, 16, 26, Strand::Forward, vec![(10, 7)]),
            make_mapping(20, 30, 32, 42, Strand::Forward, vec![(10, 7)]),
        ];
        let read = ReadInfo::new_fa_record("q1".to_string(), read.to_string());

        let mut metric = AlignedMetric::default();
        let rows = metric.compute(default_metric_context(
            read,
            Arc::new("t".to_string()),
            mappings,
            global_data,
        ));

        assert_eq!(rows.len(), 1);
        let cells: Vec<&str> = rows[0].split('\t').collect();
        assert_eq!(cells.len(), AlignedMetric::COLUMNS.len());
        assert_eq!(cells[3], "3", "segs");
        assert_eq!(cells[4], "25", "covlen = 10+5+10");
        assert_eq!(cells[10], "5", "qOvlp");
        assert_eq!(cells[13], "0:10;10:15;20:30", "mergedQrySpan");
        assert_eq!(
            cells[8], "0:10->0:10+;5:15->16:26+;20:30->32:42+",
            "oriAlignInfo 按 qstart 排序"
        );
    }

    /// 端到端：真实 aligner + fasta 参考，fwd/rev 各产出一行 45 列。
    #[test]
    fn test_e2e_aligned_metric() {
        let ref_file = format!(
            "{}/test_data/ref_Saureus_ATCC25923.m.new.corrected.fasta",
            env!("CARGO_MANIFEST_DIR")
        );
        let fa_iter = FastaFileReader::new(ref_file);
        let targets = read_fastx(fa_iter);
        let targetname2seq = targets
            .iter()
            .map(|v| (Arc::new(v.name.clone()), Arc::new(v.seq.clone())))
            .collect::<HashMap<_, _>>();
        let (tname, tseq) = {
            let (k, v) = targetname2seq.iter().next().unwrap();
            (k.clone(), v.clone())
        };
        let global_data = Arc::new(GlobalData::new(targetname2seq));

        let index_params = IndexParams::default().set_kmer(11).set_wins(1);
        let align_params = AlignParams::default()
            .set_m_score(2)
            .set_mm_score(5)
            .set_gap_open_penalty("2,24".to_string())
            .set_gap_extension_penalty("1,0".to_string());
        let mut aligner = build_aligner_from_params(
            "map-ont",
            &index_params,
            &align_params,
            &[ReadInfo::new_fa_record(
                tname.as_ref().to_string(),
                tseq.as_ref().to_string(),
            )],
        );
        aligner.options_mut().filtering.best_n = 100000;
        aligner.options_mut().filtering.pri_ratio = 0.1;

        // read = 参考 [0..300)，正向
        let fwd_read =
            ReadInfo::new_fa_record("query".to_string(), tseq.as_ref()[..300].to_string());
        let hits = adapter::align_single_query_to_targets(&aligner, &fwd_read);
        assert!(!hits.is_empty(), "fwd query 应至少产生一个 hit");

        let mut metric = AlignedMetric::default();
        let rows = metric.compute(default_metric_context(
            fwd_read,
            tname.clone(),
            hits,
            global_data.clone(),
        ));
        assert_eq!(rows.len(), 1, "aligned mode 每条 read 恰好一行");
        let cells: Vec<&str> = rows[0].split('\t').collect();
        assert_eq!(cells.len(), AlignedMetric::COLUMNS.len());
        assert_eq!(cells[1], tname.as_ref(), "rname");
        assert_eq!(cells[2], "300", "qlen");
        let covlen: usize = cells[4].parse().unwrap();
        assert!(covlen > 0, "有 hit 时 covlen > 0");
        let identity: f64 = cells[7].parse().unwrap();
        assert!(
            identity > 0.9,
            "exact read 的 identity 应接近 1，got {identity}"
        );

        // read = 同一片段的反向互补 → 仍一行，oriAlignInfo 带 '-' 链型
        let rev_read = ReadInfo::new_fa_record(
            "query".to_string(),
            String::from_utf8(reverse_complement(&tseq.as_ref().as_bytes()[..300])).unwrap(),
        );
        let hits = adapter::align_single_query_to_targets(&aligner, &rev_read);
        assert!(!hits.is_empty(), "rev query 应至少产生一个 hit");

        let mut metric = AlignedMetric::default();
        let rows = metric.compute(default_metric_context(
            rev_read,
            tname.clone(),
            hits,
            global_data,
        ));
        assert_eq!(rows.len(), 1);
        let cells: Vec<&str> = rows[0].split('\t').collect();
        assert_eq!(cells.len(), AlignedMetric::COLUMNS.len());
        assert!(
            cells[8].contains('-'),
            "reverse hit 的 oriAlignInfo 带负链标记"
        );
        let identity: f64 = cells[7].parse().unwrap();
        assert!(
            identity > 0.9,
            "exact rev read 的 identity 应接近 1，got {identity}"
        );
    }
}
