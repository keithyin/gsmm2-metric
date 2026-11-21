use std::{collections::HashMap, sync::Arc};

use hp_tr_finder::single_seq_hp_tr_finder;

use crate::{
    adapter::Mapping,
    aligned_pairs::{AlignOp, TAlignedPairs},
    global_data::GlobalData,
    metrics::{
        TMetric,
        realign::{HP_REG, RealignConfig, realign_hit},
    },
};

#[derive(Debug, Default)]
struct CounterWithAlignSpan {
    counts: [HashMap<usize, usize>; 2], // [pure, mixed]
    align_span: usize,
    eq: usize,
}

impl CounterWithAlignSpan {
    pub fn update(&mut self, misc: bool, cnt: usize, eq: usize, align_span: usize) {
        let idx = if misc { 1 } else { 0 };
        *self.counts[idx].entry(cnt).or_default() += 1;
        self.eq += eq;
        self.align_span += align_span;
    }
}

/// 统计 HP 区域 `[start, end)` 投射到 query 上的窗口内碱基计数。
///
/// 窗口切片规则（起点/终点两侧对称，口径与 hp-tr-v3 一致：两侧边界插入均计入）：
/// - 起点侧：`start > ts` 时跳到 ref 位点 `start-1` 的柱并跳过它，位于 `start-1`
///   与 `start` 之间的 leading insertion 柱（rpos=None）随之被计入；
/// - 终点侧：迭代到首个 `rpos >= end` 的柱（match 或 deletion）即停止，位于
///   `end-1` 与 `end` 之间的 trailing insertion 柱在停止前被计入。
///   （修复 Bug 5：原 take_while 实现在 rpos==end-1 柱后立即截断，trailing
///   insertion 被丢弃，与起点侧不对称。）
///
/// 计数语义：
/// - `eq`：仅 match 柱（qpos/rpos 均 Some）中 read 碱基 == target_base 的次数；
/// - `cnt`/`is_misc`：所有 qpos 为 Some 的柱（match + insertion）；
/// - `align_span`：被消费的柱数（含 insertion 柱）；deletion 柱（qpos=None）
///   只占 span，不参与任何计数。
fn count_window(
    read_seq_bytes: &[u8],
    target_base: u8,
    pairs: impl Iterator<Item = (Option<i64>, Option<i64>, AlignOp)>,
    start: usize,
    end: usize,
    ts: usize,
) -> (usize, bool, usize, usize) {
    // 起点侧：丢弃首个 rpos == start-1 的柱之前的所有柱（含 leading deletion/insertion），
    // 并跳过该柱本身；等价于 skip_while(rpos != Some(start-1)).skip(1)。
    let mut skip_to_start = start > ts;

    let mut cnt = 0;
    let mut is_misc = false;
    let mut eq = 0;
    let mut align_span = 0;

    for (qpos, rpos, _align_op) in pairs {
        if skip_to_start {
            if rpos == Some(start as i64 - 1) {
                skip_to_start = false;
            }
            continue;
        }

        if let Some(rpos) = rpos
            && rpos >= end as i64
        {
            break;
        }

        if let (Some(qpos_), Some(_rpos_)) = (qpos, rpos)
            && read_seq_bytes[qpos_ as usize] == target_base
        {
            eq += 1;
        }

        if let Some(qpos) = qpos {
            cnt += (read_seq_bytes[qpos as usize] == target_base) as usize;
            is_misc |= read_seq_bytes[qpos as usize] != target_base;
        }
        align_span += 1;
    }

    (cnt, is_misc, eq, align_span)
}

#[derive(Debug, Default)]
pub struct HpMetricV2 {
    qname: Option<String>,
    tname: Option<Arc<String>>,
    global_data: Option<Arc<GlobalData>>,
    align_infos: Vec<Mapping>,
    //Arc<String> for motif string. [fwd_counter, rev_counter]
    metric_core: HashMap<Arc<String>, CounterWithAlignSpan>,
    metric_str: Option<String>,
}

impl TMetric for HpMetricV2 {
    fn csv_header() -> String {
        let csv_header = [
            "qname".to_string(),
            "rname".to_string(),
            "motif".to_string(),
            "tag".to_string(),
            "called".to_string(),
            "num".to_string(),
            "eq".to_string(),
            "alignspan".to_string(),
        ];
        csv_header.join("\t")
    }

    fn set_qname(&mut self, qname: String) {
        self.qname = Some(qname);
    }
    fn get_qname(&self) -> &str {
        self.qname.as_ref().unwrap()
    }

    fn set_global_data(&mut self, global_data: Arc<GlobalData>) {
        self.global_data = Some(global_data);
    }
    fn get_global_data(&self) -> &GlobalData {
        self.global_data.as_ref().unwrap()
    }
    fn compute_metric(&mut self, read_info: &crate::kits::ReadInfo, reference_anchored: bool) {
        let cfg = RealignConfig {
            shrink_to_boundary: false,
            reference_anchored,
        };
        for old_align_info in &self.align_infos {
            let Some(hit) = realign_hit(read_info, old_align_info, self.get_global_data(), &cfg)
            else {
                continue;
            };

            let align_info = &hit.align;
            let read_seq_bytes = hit.read_seq.as_bytes();

            let mut match_patterns: HashMap<String, Arc<String>> = HashMap::new();
            let region2motif =
                single_seq_hp_tr_finder(&HP_REG, &mut match_patterns, &hit.target_substr);

            for v in region2motif.iter() {
                let (start, end) = (v.0.0, v.0.1);
                if start < align_info.target_start as usize || end > align_info.target_end as usize
                {
                    continue;
                }

                let motif = v.1.clone();
                let target_base = motif.as_bytes()[1];
                let expected_cnt = motif[3..].parse::<i32>().unwrap();

                let (cnt, is_misc, eq, align_span) = count_window(
                    read_seq_bytes,
                    target_base,
                    align_info.aligned_pairs(),
                    start,
                    end,
                    align_info.target_start as usize,
                );

                if (cnt as i64 - expected_cnt as i64).abs() >= 10 {
                    let (aligned_target, aligned_query) =
                        align_info.aligned_2_str(hit.target_substr.as_bytes(), read_seq_bytes);
                    let aligned_len = align_info.query_end - align_info.query_start;
                    tracing::warn!(
                        "\n query_name:{}, motif:{}, baseCnt:{}, expectedCnt:{}, QueryStartEnd:{}-{}, TargetStartEnd:{}-{}\ntarget:{}\nquery :{}",
                        read_info.name,
                        motif,
                        cnt,
                        expected_cnt,
                        old_align_info.query_start + align_info.query_start,
                        old_align_info.query_start + align_info.query_start + aligned_len,
                        old_align_info.target_start,
                        old_align_info.target_end,
                        aligned_target,
                        aligned_query
                    );
                }

                self.metric_core
                    .entry(motif.clone())
                    .or_default()
                    .update(is_misc, cnt, eq, align_span);
            }
        }
    }

    fn build_metric_str(&mut self) -> String {
        let mut result_items = vec![];
        for (key, counter) in &self.metric_core {
            for (misc_idx, cnt) in counter.counts.iter().enumerate() {
                cnt.iter().for_each(|(&called, &num)| {
                    let mut innner_items = vec![];
                    innner_items.push(self.qname.as_ref().unwrap().clone());
                    innner_items.push(self.tname.as_ref().unwrap().as_ref().clone());
                    innner_items.push(key.as_ref().clone());

                    innner_items.push(if misc_idx == 0 {
                        "pure".to_string()
                    } else {
                        "mixed".to_string()
                    });

                    innner_items.push(format!("{}", called));
                    innner_items.push(format!("{}", num));
                    innner_items.push(format!("{}", counter.eq));
                    innner_items.push(format!("{}", counter.align_span));
                    result_items.push(innner_items.join("\t"));
                });
            }
        }
        result_items.join("\n")
    }
    fn set_target_name(&mut self, target_name: Arc<String>) {
        self.tname = Some(target_name);
    }
    fn get_target_name(&self) -> &str {
        self.tname.as_ref().unwrap()
    }

    fn get_mappings_mut(&mut self) -> &mut Vec<Mapping> {
        &mut self.align_infos
    }

    fn get_metric_str_mut(&mut self) -> &mut Option<String> {
        &mut self.metric_str
    }
    fn get_metric_str(&self) -> &Option<String> {
        &self.metric_str
    }
}

#[cfg(test)]
mod test {
    use std::{collections::HashMap, sync::Arc};

    use crate::kits::{
        ReadInfo,
        fastx::{FastaFileReader, read_fastx},
    };

    use crate::{
        adapter,
        aligned_pairs::AlignOp,
        aligner::build_aligner_from_params,
        global_data::GlobalData,
        metrics::{TMetric, hp_metric_v2::HpMetricV2},
        params::{AlignParams, IndexParams},
    };

    /// 行序来自 HashMap 迭代（每次运行可能不同），排序后比较行集合。
    /// 期望值以当前 Cargo.lock 锁定的 bandmap-core 基线为准；依赖升级导致
    /// 变化时按新基线更新。
    fn assert_metric_rows(metric_str: &str, expected_rows: &[&str]) {
        let mut actual: Vec<&str> = metric_str.lines().collect();
        actual.sort_unstable();
        let mut expected: Vec<&str> = expected_rows.to_vec();
        expected.sort_unstable();
        assert_eq!(actual, expected, "metric 行集合与基线不一致");
    }

    fn ref_file() -> String {
        format!(
            "{}/test_data/ref_Saureus_ATCC25923.m.new.corrected.fasta",
            env!("CARGO_MANIFEST_DIR")
        )
    }

    #[test]
    fn test_hp_metric() {
        let fa_iter = FastaFileReader::new(ref_file());
        let targets = read_fastx(fa_iter);

        let targetname2seq = targets
            .iter()
            .map(|v| (Arc::new(v.name.clone()), Arc::new(v.seq.clone())))
            .collect::<HashMap<_, _>>();

        let _target_seq = targetname2seq.values().take(1).collect::<Vec<_>>()[0].clone();

        let global_data = Arc::new(GlobalData::new(targetname2seq));

        let index_params = IndexParams {
            kmer: Some(11),
            wins: Some(1),
        };

        let align_params = AlignParams::default();
        let mut aligner =
            build_aligner_from_params("map-ont", &index_params, &align_params, &targets);

        aligner.options_mut().filtering.best_n = 100000;
        aligner.options_mut().filtering.pri_ratio = 0.1;
        let seq = b"CCCTCCCAATGATGTATAAACAATTATATGTTATGTTCATTATCCTACAAATCTCCAACATTGATGATTGGGCACAACAATTTTACCTGTTTAATAAGGTGAACAAAAAACAAACGAAAAAGGTGATAACAATGAACAACATTTACATTAGGAAATCCAAATTAACTGTAACTCATGTCCATGAAGTGAAAGCCGGTATTAAACCGACATCGGTTGTCGACAGTGTTCAATAGAAATCAAGAAATGATTATGGAAAAAAGATGCACTGTGGAAATGCGCTGGCGAAAAATTATATATTAATATTTTCAGACCAAAAAAGATGGCATTTCCCCTGTAAGTTATGTCTGCAGATACTTAAACGGTAAAGATAATAAAGTAAATCACATTGGGTGCCCTTTGCCAACATTAGGCGCATTCCCGACATCAGTTTTACACCTGAAAGAATCACCAGACCCAGATTTTCGGTGCCAATGATTATGTTGGTAAGTTAAAGTTGCAATTACGGCGGTAGTGACAAATCAAAGGCGTCTTTATCTCCATGTCAAAAAGAGCAAGCGGAAGATTATTACAAAGTGATGGTGGGCAGCAAATCAGTCATGAGTAATGGAAACATCGGACAAATGGCGTTTCTTATCTTGCGGTGACTCAATGGTGGGTCGCAATCATTAAATCCACCACATTAAAGCAATGATCCTTGGGAAGGCTTAAATGATATGTATGTATAGAGAAGTAAGGCTTTCACGGAGTATAACGGAATACTGGCTTTTAATCGTTTCTGGACTCAAGTATTTTGCGAGATGACAGATAATCCCAAATATCAAGATTTTGATTCAACCAACAAGAACACCCTTTGTCGATGATTTTGAAACAGCGCAAGTGGAGCCATTATCACAATTAGAAAAAAACACCTCTACTAAAACATGTGCTAGTTGGTCTACACAGTTGCACAACCGGGCTCTTTTGATGGTTTAAAACAAGCAGCATCAGAAGAAAAAGGCTATATGGTGCCATGGCGTAAAGAGTGGGAAATATCTCTCTCAAAACAAACAACGGAGGAGAGGAAAAAAGAGAGAGATAACTTTTTCCTCTTTACTCCCATGCACATATAGCCATTTTTCTTCTGATGCTGCTTGTTTAATCCTTCAAAAGAGCACGTTGTGAAAAAAAAAAAAAAGGAAAAAAAAAAAAAACCTTGTGTAGCCAACTAGCACATGTTAGTAAGAGTGTTTAATTTTGTGATATGGCACCTTGACGCTGTGTCCAAAAATATCGAACCAAAGGGTGTTCTTGTTGTGCCTGAATCAAATCTTCGAATTGGTTACGTCCACTTCCAAAAAAAATACTTGAGTCCCAAAACGATAAACCAGTATCCGTTTATACTCCTGAAAGGCTACTTCCTCTATACATATCATTTAAGCCTGCCCAAGAATCAATTGCTTTTTAATGTGGTGGATTTAAATGATGCGAACCCACCCATTGAGTCACGCAAAGATAAGAAACGCCATTTTCCGATGTTTCCATTACTCCATGACTTGATTGCTGTCACCAATCACTTCCGTAATGATACTTCGCTTCTCTTTTTGACCATGGAGATAAGACGCCTTTTGAATTTGTCACTACCGCGTAATGCACCTTTAACTACAACATAATCATTTGGCACCATATCCTGGTCGGTGATTTCTTCAGGTGAAAACTAGTGTCGGAAAATGGCGCCTATTGGCCAAAGGGCCATATTTGTGATTTTAAGCTTATTATCTTTACCGTAAGTATCTGCAGAACATAACTACAGGGAATTGCCATCTTTTATTTGGGTCTGAAAATATTAATATATAATTTTCGCCATCGCGCAATTTCCACAGTGGACATCTTTTTCGCATAATCATTTCTGATTTCCATATTGAACCTGTCGACAACGATGTGTTAATCCGGTTTACTTCCATTTGACATGAGTTACAAGTTAATTTGGATTCCTAGTAAATGTTGGTTCATTGTTAAACACCCTTTTCTGTTTGTTTGTATCACCATATTAAACAGGTAAACCATTGTGTGCCAATCATCAATGTTGGATTTGTAGATAATGAACATAACATTTATATTGTTATAATCATTGGGGAGACTTGAATACAAATGACTATATCTCTCTTCAAACAACAAACACCACCACAGACAACAGCAGAAAGGCAA";
        let fwd_part = &seq[..100];
        let rev_part = &seq[1074..1074 + 100];
        let fwd_query_record = ReadInfo::new_fa_record(
            "query".to_string(),
            String::from_utf8(fwd_part.to_vec()).unwrap(),
        );

        let hits = adapter::align_single_query_to_targets(&aligner, &fwd_query_record);
        assert!(!hits.is_empty(), "fwd query 应至少产生一个 hit");

        let mut metric = HpMetricV2::default();
        metric.set_qname("query".to_string());
        metric.set_target_name(Arc::new("target".to_string()));
        metric.set_global_data(global_data.clone());
        metric.set_mappings(hits);
        metric.compute_metric(&fwd_query_record, false);
        metric.set_metric_str();
        println!("metric:\n{}", metric.get_metric_str().as_ref().unwrap());
        assert_metric_rows(
            metric.get_metric_str().as_ref().unwrap(),
            &[
                "query\ttarget\t(T)3\tpure\t3\t2\t6\t6",
                "query\ttarget\t(A)3\tpure\t3\t2\t7\t9",
                "query\ttarget\t(A)3\tpure\t1\t1\t7\t9",
            ],
        );

        let rev_query_record = ReadInfo::new_fa_record(
            "query".to_string(),
            String::from_utf8(rev_part.to_vec()).unwrap(),
        );

        let hits = adapter::align_single_query_to_targets(&aligner, &rev_query_record);
        assert!(!hits.is_empty(), "rev query 应至少产生一个 hit");

        let mut metric = HpMetricV2::default();
        metric.set_qname("query".to_string());
        metric.set_target_name(Arc::new("target".to_string()));
        metric.set_global_data(global_data.clone());
        metric.set_mappings(hits);
        metric.compute_metric(&rev_query_record, false);
        metric.set_metric_str();
        println!("metric:\n{}", metric.get_metric_str().as_ref().unwrap());
        assert_metric_rows(
            metric.get_metric_str().as_ref().unwrap(),
            &[
                "query\ttarget\t(A)3\tpure\t2\t1\t2\t3",
                "query\ttarget\t(A)4\tpure\t4\t1\t4\t4",
                "query\ttarget\t(C)3\tpure\t1\t1\t1\t3",
                "query\ttarget\t(T)3\tpure\t3\t2\t9\t11",
                "query\ttarget\t(T)3\tpure\t5\t1\t9\t11",
                "query\ttarget\t(T)5\tpure\t5\t1\t5\t5",
            ],
        );
    }

    // ---- count_window 回归测试（手造 aligned_pairs，Bug 5）----

    fn pair(q: i64, r: i64) -> (Option<i64>, Option<i64>, AlignOp) {
        (Some(q), Some(r), AlignOp::Equal(1))
    }
    fn ins(q: i64) -> (Option<i64>, Option<i64>, AlignOp) {
        (Some(q), None, AlignOp::Ins(1))
    }
    fn del(r: i64) -> (Option<i64>, Option<i64>, AlignOp) {
        (None, Some(r), AlignOp::Del(1))
    }

    /// 区域 [2,6)，ref 位点 5（end-1）与 6（end）之间有 1 个 trailing insertion。
    /// 修复前 take_while 在 ref5 柱后截断，insertion 被丢弃（cnt=4）；修复后计入（cnt=5）。
    #[test]
    fn test_count_window_trailing_insertion_counted() {
        let read_seq = b"AAAAAAAA"; // qpos 0..7 全 A
        let pairs = vec![
            pair(0, 0),
            pair(1, 1),
            pair(2, 2),
            pair(3, 3),
            pair(4, 4),
            pair(5, 5),
            ins(6), // ref5 与 ref6 之间的插入
            pair(7, 6),
        ];

        let (cnt, is_misc, eq, align_span) =
            super::count_window(read_seq, b'A', pairs.into_iter(), 2, 6, 0);

        assert_eq!(cnt, 5, "trailing insertion 碱基应计入 called");
        assert!(!is_misc);
        assert_eq!(eq, 4, "eq 仅统计 match 柱（ref2..ref5）");
        assert_eq!(align_span, 5, "4 个 match 柱 + 1 个 insertion 柱");
    }

    /// trailing insertion 碱基与 target_base 不同时应置 is_misc（修复前被丢弃则不影响 misc）。
    #[test]
    fn test_count_window_trailing_insertion_misc() {
        let read_seq = b"AAAAAACA"; // qpos 6 是 C
        let pairs = vec![
            pair(0, 0),
            pair(1, 1),
            pair(2, 2),
            pair(3, 3),
            pair(4, 4),
            pair(5, 5),
            ins(6), // C 插在 ref5 与 ref6 之间
            pair(7, 6),
        ];

        let (cnt, is_misc, _, _) = super::count_window(read_seq, b'A', pairs.into_iter(), 2, 6, 0);

        assert_eq!(cnt, 4, "C 碱基不计入 called，但参与 is_misc 判定");
        assert!(is_misc, "边界插入的非 motif 碱基应将区域标为 mixed");
    }

    /// 起点侧 leading insertion 计入行为与修复前保持一致（本就计入，非本次改动目标）。
    #[test]
    fn test_count_window_leading_insertion_unchanged() {
        let read_seq = b"AAAAAAA"; // qpos 0..6 全 A
        let pairs = vec![
            pair(0, 0),
            pair(1, 1),
            ins(2), // ref1 与 ref2 之间的插入（leading，相对窗口起点 start=2）
            pair(3, 2),
            pair(4, 3),
            pair(5, 4),
            pair(6, 5),
        ];

        let (cnt, is_misc, eq, align_span) =
            super::count_window(read_seq, b'A', pairs.into_iter(), 2, 6, 0);

        assert_eq!(cnt, 5, "1 个 leading insertion + 4 个 match 碱基");
        assert!(!is_misc);
        assert_eq!(eq, 4, "eq 仅统计 match 柱");
        assert_eq!(align_span, 5);
    }

    /// 无边界插入时与旧 take_while 逻辑完全一致：窗口 [2,6) 恰好消费 ref2..ref5，
    /// ref6 处的 deletion 柱触发停止且不计入 span。
    #[test]
    fn test_count_window_no_boundary_indel() {
        // 纯 A 区域
        let read_seq = b"AAAAAA";
        let pairs = vec![
            pair(0, 0),
            pair(1, 1),
            pair(2, 2),
            pair(3, 3),
            pair(4, 4),
            pair(5, 5),
            del(6), // ref6 柱应触发 break
        ];

        let (cnt, is_misc, eq, align_span) =
            super::count_window(read_seq, b'A', pairs.clone().into_iter(), 2, 6, 0);

        assert_eq!(cnt, 4);
        assert!(!is_misc);
        assert_eq!(eq, 4);
        assert_eq!(align_span, 4, "deletion 柱被 break 前不消费");

        // 含 1 个替换碱基：qpos3 为 C
        let read_seq_mixed = b"AAACAA";
        let (cnt, is_misc, eq, _) =
            super::count_window(read_seq_mixed, b'A', pairs.into_iter(), 2, 6, 0);

        assert_eq!(cnt, 3);
        assert!(is_misc);
        assert_eq!(eq, 3);
    }

    /// 起点侧跳过逻辑：窗口起点之前的 deletion 柱与 match 柱被正确丢弃。
    #[test]
    fn test_count_window_skip_before_start() {
        let read_seq = b"AAAAAA";
        let pairs = vec![
            del(0), // start-1 之前的 deletion，应被跳过
            pair(0, 0),
            del(1), // start-1 柱本身，应被跳过
            ins(1), // start-1 与 start 之间的 leading insertion，应被计入
            pair(1, 2),
            pair(2, 3),
            pair(3, 4),
            pair(4, 5),
            del(6),
        ];

        let (cnt, is_misc, _, align_span) =
            super::count_window(read_seq, b'A', pairs.into_iter(), 2, 6, 0);

        assert_eq!(cnt, 5, "1 个 leading insertion + 4 个 match 碱基");
        assert!(!is_misc);
        assert_eq!(align_span, 5);
    }
}
