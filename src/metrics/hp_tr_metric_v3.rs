use std::{collections::HashMap, sync::Arc};

use gskits::dna::reverse_complement;
use hp_tr_finder::single_seq_hp_tr_finder;
use mm2::minimap2::Mapping;

use crate::{
    aligned_pairs::TAlignedPairs,
    global_data::{GlobalData, GlobalDataKey, GlobalDataValue},
    metrics::{
        TMetric,
        hp_tr_metric_v2::HP_TR_REG,
        hp_tr_tools::{do_align_4_homo, get_target_substr},
    },
};

#[derive(Debug, Default)]
struct Counter {
    ref_repeats: usize,
    counts: [HashMap<usize, usize>; 2], // [pure, mixed]: called -> num
    align_span: usize,
    eq: usize,
}

impl Counter {
    pub fn update(&mut self, ref_repeats: usize, misc: bool, cnt: usize, eq: usize, align_span: usize) {
        let idx = if misc { 1 } else { 0 };
        *self.counts[idx].entry(cnt).or_default() += 1;
        self.ref_repeats = ref_repeats;
        self.eq += eq;
        self.align_span += align_span;
    }
}

/// 从 motif 字符串中解析出 unit 序列、unit 长度和 reference 的 repeats 数。
/// 例如 "(T)8" -> ("T", 1, 8), "(CA)5" -> ("CA", 2, 5)。
fn parse_motif(motif: &str) -> (&str, usize, usize) {
    let unit = &motif[motif.find('(').unwrap() + 1..motif.rfind(')').unwrap()];
    let repeats = motif[motif.rfind(')').unwrap() + 1..]
        .parse::<usize>()
        .unwrap();
    (unit, unit.len(), repeats)
}

/// 统计 query 子串中 motif 的 repeats 数（完整 unit 计数，贪婪左起非重叠滑动）。
///
/// 从最左位置开始，找到第一个完整匹配 motif 的 unit，计数后跳到该 unit 之后，
/// 重复"找下一个 match"。前导垃圾/中间错配不会因固定相位而吞掉后续 unit。
///
/// - `called`：找到的完整匹配 unit 个数（非重叠）；
/// - `is_misc`：窗口内存在未被任何匹配 unit 覆盖的碱基（前导/尾随/夹缝）→ mixed；
/// - `eq`：被计为完整 unit 覆盖的碱基数（= unit_len × called）。
fn count_repeats(motif_unit: &[u8], query_substr: &[u8]) -> (usize, bool, usize) {
    let unit_len = motif_unit.len();
    let mut called = 0;
    let mut i = 0;
    while i + unit_len <= query_substr.len() {
        if &query_substr[i..i + unit_len] == motif_unit {
            called += 1;
            i += unit_len; // 命中后跳到该 unit 之后，重复"找下一个 match"
        } else {
            i += 1; // 未命中则滑动一格，继续找第一个 match
        }
    }
    let eq = called * unit_len;
    let is_misc = eq != query_substr.len();
    (called, is_misc, eq)
}

/// 基于 aligned_pairs() 构建 ref_position -> query_position 的映射。
/// 数组下标即 ref 局部坐标（index i <-> 局部 ref 位点 target_start + i），
/// 长度固定为对齐覆盖的 ref 跨度 (target_end - target_start)。
/// indel 处理规则：
///   - deletion (ref 有, query 无)：ref_pos 指向下一个未 deletion 的 query 位点（query 游标不变）
///   - insertion (query 有, ref 无)：ref_pos 指向当前 ref 游标 rcur（上一个 match 的下一个 query 位点）
///
/// query 游标初值为 align_info.query_start()：aligned_pairs() 的 query 坐标从该位置开始
/// 递增（见 aligned_pairs::TAlignedPairs::aligned_pairs），所以前导 deletion 应映射到
/// 第一个 match 的 query 位点 query_start，而非 0。
fn build_ref2qpos(align_info: &Mapping) -> Vec<usize> {
    let ts = align_info.target_start() as usize;
    let te = align_info.target_end() as usize;
    let mut ref2qpos = vec![0usize; te - ts];
    let mut rcur = ts;
    let mut qcur = align_info.query_start() as usize;
    for (qpos, rpos, _op) in align_info.aligned_pairs() {
        match (qpos, rpos) {
            (Some(q), Some(r)) => {
                ref2qpos[r as usize - ts] = q as usize;
                rcur = r as usize + 1;
                qcur = q as usize + 1;
            }
            (Some(q), None) => {
                // insertion 位于 ref 游标处；若后面跟 match，会被该 match 覆盖
                if rcur < te {
                    ref2qpos[rcur - ts] = q as usize;
                }
                qcur = q as usize + 1;
            }
            (None, Some(r)) => {
                ref2qpos[r as usize - ts] = qcur;
                rcur = r as usize + 1;
            }
            (None, None) => unreachable!(),
        }
    }
    ref2qpos
}

#[derive(Debug, Default)]
pub struct HpTrMetricV3 {
    qname: Option<String>,
    tname: Option<Arc<String>>,
    global_data: Option<Arc<GlobalData>>,
    align_infos: Vec<Mapping>,
    //Arc<String> for motif string.
    metric_core: HashMap<Arc<String>, Counter>,
    metric_str: Option<String>,
}

impl TMetric for HpTrMetricV3 {
    fn csv_header() -> String {
        let csv_header = vec![
            "qname".to_string(),
            "rname".to_string(),
            "motif".to_string(),
            "tag".to_string(),
            "ref".to_string(),
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

    fn compute_metric(&mut self, read_info: &mm2::gskits::ds::ReadInfo, reference_anchored: bool) {
        if self.align_infos.is_empty() {
            return;
        }

        let target_name = self.align_infos[0].target_name.as_ref().unwrap();
        let global_data_value = self
            .global_data
            .as_ref()
            .unwrap()
            .get(GlobalDataKey::TargetName2SeqAndRev);
        let [target_seq_fwd, target_seq_rev] = match global_data_value {
            GlobalDataValue::TargetName2SeqAndRev(v) => v.get().unwrap().get(target_name).unwrap(),
            _ => panic!(""),
        };

        let global_data_value = self
            .global_data
            .as_ref()
            .unwrap()
            .get(GlobalDataKey::TargetRegion2Motif4HpTr);
        let target_region_2_motif_all = match global_data_value {
            GlobalDataValue::TargetRegion2Motif4HpTr(v) => {
                v.get().unwrap().get(target_name).unwrap()
            }
            _ => panic!("TargetRegion2Motif4HpTr not found"),
        };

        for old_align_info in &self.align_infos {
            let mut old_target_start = old_align_info.target_start as usize;
            let mut old_target_end = old_align_info.target_end as usize;
            // 如果比对的 reference_start, reference_end 正好处于 hp-str 区域的中心位置，那就从两头往中间缩
            if let Some(regions) = target_region_2_motif_all.get(&old_target_start) {
                let new_start = regions.iter().max_by_key(|v| v.0.1).unwrap().0.1;
                tracing::info!("old_target_start: {old_target_start} -> {new_start}");
                old_target_start = new_start;
            }
            if old_target_end == 0 {
                tracing::warn!("aligned target end is zero.");
                continue;
            }
            if let Some(regions) = target_region_2_motif_all.get(&(old_target_end - 1)) {
                let new_end = regions.iter().min_by_key(|v| v.0.0).unwrap().0.0;
                tracing::info!("old_target_end: {old_target_end} -> {new_end}");
                old_target_end = new_end;
            }

            let target_substr = get_target_substr(
                old_target_start,
                old_target_end as usize,
                old_align_info.is_reverse(),
                target_seq_fwd,
                target_seq_rev,
            );

            let read_seq = &read_info.seq
                [old_align_info.query_start as usize..old_align_info.query_end as usize];

            let (read_seq, target_substr) = if old_align_info.is_reverse() && reference_anchored {
                let read_rev = reverse_complement(read_seq.as_bytes());
                let read_rev = String::from_utf8(read_rev).unwrap();
                let target_rev = reverse_complement(target_substr.as_bytes());
                let target_rev = String::from_utf8(target_rev).unwrap();
                (read_rev, target_rev)
            } else {
                (read_seq.to_string(), target_substr.to_string())
            };

            let (read_seq, target_substr) = (&read_seq, &target_substr);
            let read_seq_bytes = read_seq.as_bytes();

            let align_info = do_align_4_homo(read_seq, target_substr);
            if align_info.is_none() {
                tracing::warn!(
                    "no aligned result. QueryName:{}. QueryStartEnd:{}-{}, TargetStartEnd:{}-{}, strand:{:?}",
                    read_info.name,
                    old_align_info.query_start,
                    old_align_info.query_end,
                    old_align_info.target_start,
                    old_align_info.target_end,
                    old_align_info.strand
                );
                continue;
            }
            let align_info = align_info.unwrap();

            if !matches!(align_info.strand, mm2::minimap2::Strand::Forward) {
                tracing::warn!(
                    "still reverse aligment. QueryName:{}. QueryStartEnd:{}-{}",
                    read_info.name,
                    old_align_info.query_start,
                    old_align_info.query_end,
                );
                continue;
            }

            let mut match_patterns: HashMap<String, Arc<String>> = HashMap::new();
            let region2motif =
                single_seq_hp_tr_finder(&HP_TR_REG, &mut match_patterns, target_substr);

            // 基于 aligned_pairs() 构建 ref_position -> query_position 的映射，
            // 数组下标即 ref 局部坐标（index i <-> 局部 ref 位点 target_start + i），
            // indel 规则见 build_ref2qpos。
            let ts = align_info.target_start as usize;
            let te = align_info.target_end as usize;
            let ref2qpos = build_ref2qpos(&align_info);

            for v in region2motif.iter() {
                let (start, end) = (v.0.0, v.0.1);
                if start < ts || end > te {
                    continue;
                }

                let motif = v.1.clone();
                let (motif_unit, _unit_len, ref_repeats) = parse_motif(motif.as_ref());

                // motif 窗口用两侧 flanking 公式：前导/尾随插入（在 motif 边界两侧、
                // ref 侧无对齐柱）被 [ref2qpos[start-1-ts]+1, ref2qpos[start-ts]] 和
                // [ref2qpos[end-1-ts], ref2qpos[end-ts]-1] 的区间覆盖；
                // 当 motif 顶到对齐末端（end == te）时，尾随插入只能靠回落 query_end-1 捕获。
                let q_start = if start > ts {
                    ref2qpos[start - 1 - ts] + 1
                } else {
                    align_info.query_start as usize
                };
                let q_end = if end < te {
                    ref2qpos[end - ts].saturating_sub(1)
                } else {
                    (align_info.query_end - 1) as usize
                };
                // 该 motif 区域在 query 侧没有碱基（整个被删掉），跳过
                if q_start > q_end {
                    continue;
                }
                let query_substr = &read_seq_bytes[q_start..=q_end];

                let (called, is_misc, eq) = count_repeats(motif_unit.as_bytes(), query_substr);

                self.metric_core.entry(motif.clone()).or_default().update(
                    ref_repeats,
                    is_misc,
                    called,
                    eq,
                    query_substr.len(),
                );
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

                    innner_items.push(counter.ref_repeats.to_string());
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

    fn get_mappings_mut(&mut self) -> &mut Vec<mm2::minimap2::Mapping> {
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

    use gskits::{
        ds::ReadInfo,
        fastx_reader::{fasta_reader::FastaFileReader, read_fastx},
    };
    use mm2::{
        align_single_query_to_targets, build_aligner,
        params::{AlignParams, IndexParams, MapParams, OupParams},
    };

    use crate::{
        global_data::GlobalData,
        metrics::{
            TMetric, hp_tr_metric_v3::HpTrMetricV3, hp_tr_metric_v3::build_ref2qpos,
            hp_tr_metric_v3::count_repeats,
        },
    };

    /// 前导 deletion 时，ref2qpos[0] 应映射到第一个 match 的 query 位点 query_start，
    /// 而不是写死的 0。
    #[test]
    fn test_build_ref2qpos() {
        use mm2::minimap2::{Alignment, Mapping, Strand};

        fn make_mapping(query_start: i32) -> Mapping {
            Mapping {
                query_name: None,
                query_len: None,
                query_start,
                query_end: query_start + 1,
                strand: Strand::Forward,
                target_name: None,
                target_len: 22,
                target_start: 0,
                target_end: 2,
                match_len: 0,
                block_len: 0,
                mapq: 0,
                is_primary: true,
                is_supplementary: false,
                alignment: Some(Alignment {
                    nm: 0,
                    cigar: Some(vec![(1, 2), (1, 7)]), // 1xDel + 1xEqual
                    cigar_str: None,
                    md: None,
                    cs: None,
                    alignment_score: None,
                }),
            }
        }

        // query_start=5：aligned_pairs() 产出 (None, Some(0)) 然后 (Some(5), Some(1))
        let mapping = make_mapping(5);
        assert_eq!(build_ref2qpos(&mapping), vec![5, 5]);

        // 对照：query_start=0 同 cigar，回归旧行为
        let mapping = make_mapping(0);
        assert_eq!(build_ref2qpos(&mapping), vec![0, 0]);
    }

    #[test]
    fn test_hp_tr_metric() {
        let ref_file = "test_data/ref_Saureus_ATCC25923.m.new.corrected.fasta";
        let fa_iter = FastaFileReader::new(ref_file.to_string());
        let targets = read_fastx(fa_iter);

        let targetname2seq = targets
            .iter()
            .map(|v| (Arc::new(v.name.clone()), Arc::new(v.seq.clone())))
            .collect::<HashMap<_, _>>();

        let global_data = Arc::new(GlobalData::new(targetname2seq));

        let mut index_params = IndexParams::default();
        index_params.kmer = Some(11);
        index_params.wins = Some(1);

        let align_params = AlignParams::default()
            .set_m_score(2)
            .set_mm_score(5)
            .set_gap_open_penalty("2,24".to_string())
            .set_gap_extension_penalty("1,0".to_string());
        let mut aligners = build_aligner(
            "map-ont",
            &index_params,
            &MapParams::default(),
            &align_params,
            &OupParams::default(),
            &targets,
            10,
        );

        aligners[0].mapopt.best_n = 100000;
        aligners[0].mapopt.pri_ratio = 0.1;
        let seq = b"CCCTCCCAATGATGTATAAACAATTATATGTTATGTTCATTATCCTACAAATCTCCAACATTGATGATTGGGCACAACAATTTTACCTGTTTAATAAGGTGAACAAAAAACAAACGAAAAAGGTGATAACAATGAACAACATTTACATTAGGAAATCCAAATTAACTGTAACTCATGTCCATGAAGTGAAAGCCGGTATTAAACCGACATCGGTTGTCGACAGTGTTCAATAGAAATCAAGAAATGATTATGGAAAAAAGATGCACTGTGGAAATGCGCTGGCGAAAAATTATATATTAATATTTTCAGACCAAAAAAGATGGCATTTCCCCTGTAAGTTATGTCTGCAGATACTTAAACGGTAAAGATAATAAAGTAAATCACATTGGGTGCCCTTTGCCAACATTAGGCGCATTCCCGACATCAGTTTTACACCTGAAAGAATCACCAGACCCAGATTTTCGGTGCCAATGATTATGTTGGTAAGTTAAAGTTGCAATTACGGCGGTAGTGACAAATCAAAGGCGTCTTTATCTCCATGTCAAAAAGAGCAAGCGGAAGATTATTACAAAGTGATGGTGGGCAGCAAATCAGTCATGAGTAATGGAAACATCGGACAAATGGCGTTTCTTATCTTGCGGTGACTCAATGGTGGGTCGCAATCATTAAATCCACCACATTAAAGCAATGATCCTTGGGAAGGCTTAAATGATATGTATGTATAGAGAAGTAAGGCTTTCACGGAGTATAACGGAATACTGGCTTTTAATCGTTTCTGGACTCAAGTATTTTGCGAGATGACAGATAATCCCAAATATCAAGATTTTGATTCAACCAACAAGAACACCCTTTGTCGATGATTTTGAAACAGCGCAAGTGGAGCCATTATCACAATTAGAAAAAAACACCTCTACTAAAACATGTGCTAGTTGGTCTACACAGTTGCACAACCGGGCTCTTTTGATGGTTTAAAACAAGCAGCATCAGAAGAAAAAGGCTATATGGTGCCATGGCGTAAAGAGTGGGAAATATCTCTCTCAAAACAAACAACGGAGGAGAGGAAAAAAGAGAGAGATAACTTTTTCCTCTTTACTCCCATGCACATATAGCCATTTTTCTTCTGATGCTGCTTGTTTAATCCTTCAAAAGAGCACGTTGTGAAAAAAAAAAAAAAGGAAAAAAAAAAAAAACCTTGTGTAGCCAACTAGCACATGTTAGTAAGAGTGTTTAATTTTGTGATATGGCACCTTGACGCTGTGTCCAAAAATATCGAACCAAAGGGTGTTCTTGTTGTGCCTGAATCAAATCTTCGAATTGGTTACGTCCACTTCCAAAAAAAATACTTGAGTCCCAAAACGATAAACCAGTATCCGTTTATACTCCTGAAAGGCTACTTCCTCTATACATATCATTTAAGCCTGCCCAAGAATCAATTGCTTTTTAATGTGGTGGATTTAAATGATGCGAACCCACCCATTGAGTCACGCAAAGATAAGAAACGCCATTTTCCGATGTTTCCATTACTCCATGACTTGATTGCTGTCACCAATCACTTCCGTAATGATACTTCGCTTCTCTTTTTGACCATGGAGATAAGACGCCTTTTGAATTTGTCACTACCGCGTAATGCACCTTTAACTACAACATAATCATTTGGCACCATATCCTGGTCGGTGATTTCTTCAGGTGAAAACTAGTGTCGGAAAATGGCGCCTATTGGCCAAAGGGCCATATTTGTGATTTTAAGCTTATTATCTTTACCGTAAGTATCTGCAGAACATAACTACAGGGAATTGCCATCTTTTATTTGGGTCTGAAAATATTAATATATAATTTTCGCCATCGCGCAATTTCCACAGTGGACATCTTTTTCGCATAATCATTTCTGATTTCCATATTGAACCTGTCGACAACGATGTGTTAATCCGGTTTACTTCCATTTGACATGAGTTACAAGTTAATTTGGATTCCTAGTAAATGTTGGTTCATTGTTAAACACCCTTTTCTGTTTGTTTGTATCACCATATTAAACAGGTAAACCATTGTGTGCCAATCATCAATGTTGGATTTGTAGATAATGAACATAACATTTATATTGTTATAATCATTGGGGAGACTTGAATACAAATGACTATATCTCTCTTCAAACAACAAACACCACCACAGACAACAGCAGAAAGGCAA";
        let fwd_part = &seq[..100];
        let rev_part = &seq[1074..1074 + 100];
        let fwd_query_record = ReadInfo::new_fa_record(
            "query".to_string(),
            String::from_utf8(fwd_part.to_vec()).unwrap(),
        );

        let hits = align_single_query_to_targets(&fwd_query_record, &aligners);

        let mut metric = HpTrMetricV3::default();
        metric.set_qname("query".to_string());
        metric.set_target_name(Arc::new("target".to_string()));
        metric.set_global_data(global_data.clone());
        metric.set_mappings(hits);
        metric.compute_metric(&fwd_query_record, false);
        metric.set_metric_str();
        println!("metric:\n{}", metric.get_metric_str().as_ref().unwrap());

        let rev_query_record = ReadInfo::new_fa_record(
            "query".to_string(),
            String::from_utf8(rev_part.to_vec()).unwrap(),
        );

        let hits = align_single_query_to_targets(&rev_query_record, &aligners);

        let mut metric = HpTrMetricV3::default();
        metric.set_qname("query".to_string());
        metric.set_target_name(Arc::new("target".to_string()));
        metric.set_global_data(global_data.clone());
        metric.set_mappings(hits);
        metric.compute_metric(&rev_query_record, false);
        metric.set_metric_str();
        println!("metric:\n{}", metric.get_metric_str().as_ref().unwrap());
    }

    /// 在 (AC)4 @174 位点构造 read，跑 HpTrMetricV3 并断言 (ref, called, tag)。
    /// 参考序列从 fasta 读取（pos 174..182 = "ACACACAC"）。
    #[test]
    fn test_str_called_counting() {
        let ref_file = "test_data/ref_Saureus_ATCC25923.m.new.corrected.fasta";
        let fa_iter = FastaFileReader::new(ref_file.to_string());
        let targets = read_fastx(fa_iter);

        let targetname2seq = targets
            .iter()
            .map(|v| (Arc::new(v.name.clone()), Arc::new(v.seq.clone())))
            .collect::<HashMap<_, _>>();
        let ref_seq = targetname2seq.values().next().unwrap().clone();
        let ref_seq = ref_seq.as_bytes();
        // 验证 (AC)4 位点
        assert_eq!(&ref_seq[174..182], b"ACACACAC");

        let global_data = Arc::new(GlobalData::new(targetname2seq));

        let mut index_params = IndexParams::default();
        index_params.kmer = Some(11);
        index_params.wins = Some(1);

        let align_params = AlignParams::default()
            .set_m_score(2)
            .set_mm_score(5)
            .set_gap_open_penalty("2,24".to_string())
            .set_gap_extension_penalty("1,0".to_string());
        let mut aligners = build_aligner(
            "map-ont",
            &index_params,
            &MapParams::default(),
            &align_params,
            &OupParams::default(),
            &targets,
            10,
        );
        aligners[0].mapopt.best_n = 100000;
        aligners[0].mapopt.pri_ratio = 0.1;

        // read = ref[0..300]，但把 [174,182) 替换为指定序列
        let build_read = |replacement: &[u8]| -> ReadInfo {
            let mut read = ref_seq[0..174].to_vec();
            read.extend_from_slice(replacement);
            read.extend_from_slice(&ref_seq[182..300]);
            ReadInfo::new_fa_record(
                "query".to_string(),
                String::from_utf8(read).unwrap(),
            )
        };

        // 3 个用例：(AC)6 扩张, (AC)3 收缩, (AC)4 + 1 个 substitution
        let cases: Vec<(&str, &[u8], usize, bool)> = vec![
            ("read_ac6_exp", b"ACACACACACAC", 6, false),
            ("read_ac3_contr", b"ACACAC", 3, false),
            ("read_ac4_subst", b"ACACATAC", 3, true),
        ];

        for (qname, replacement, exp_called, exp_misc) in cases {
            let query_record = build_read(replacement);
            let hits = align_single_query_to_targets(&query_record, &aligners);

            let mut metric = HpTrMetricV3::default();
            metric.set_qname(qname.to_string());
            metric.set_target_name(Arc::new("target".to_string()));
            metric.set_global_data(global_data.clone());
            metric.set_mappings(hits);
            metric.compute_metric(&query_record, false);
            metric.set_metric_str();

            let metric_str = metric.get_metric_str().as_ref().unwrap();
            println!("{qname} metric:\n{metric_str}\n");

            // 找到 (AC)4 的行
            let row = metric_str
                .lines()
                .find(|line| line.contains("\t(AC)4\t"))
                .unwrap_or_else(|| panic!("{qname}: no (AC)4 row. metric:\n{metric_str}"));
            let cols: Vec<&str> = row.split('\t').collect();
            // cols: qname rname motif tag ref called num eq alignspan
            assert_eq!(cols[4], "4", "{qname}: ref");
            assert_eq!(cols[5], exp_called.to_string(), "{qname}: called");
            assert_eq!(
                cols[3],
                if exp_misc { "mixed" } else { "pure" },
                "{qname}: tag"
            );
        }
    }

    /// 直接测 count_repeats 的贪婪滑动计数逻辑。
    #[test]
    fn test_count_repeats() {
        // (AC)6 纯扩张
        assert_eq!(count_repeats(b"AC", b"ACACACACACAC"), (6, false, 12));
        // (AC)4 + 1 错配（AT）-> 3 个完整 unit，eq 只算完整 unit 覆盖的碱基数
        assert_eq!(count_repeats(b"AC", b"ACACATAC"), (3, true, 6));
        // 纯 HP：unit=1，T8
        assert_eq!(count_repeats(b"T", b"TTTTTTTT"), (8, false, 8));
        // HP 带 1 个杂散碱基 -> mixed
        assert_eq!(count_repeats(b"T", b"TTTCTT"), (5, true, 5));
        // 前导 1bp 垃圾：相位偏移不影响后续 unit 匹配
        assert_eq!(count_repeats(b"AC", b"TACACACACAC"), (5, true, 10));
        // 中间夹缝 1bp：跳过夹缝后仍能匹配后续 unit
        assert_eq!(count_repeats(b"AC", b"ACACACXAC"), (4, true, 8));
        // unit=4 非重叠
        assert_eq!(count_repeats(b"TATA", b"TATATATA"), (2, false, 8));
    }
}
