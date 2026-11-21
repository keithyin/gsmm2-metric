use std::{collections::HashMap, panic, sync::Arc};

use hp_tr_finder::single_seq_hp_tr_finder;

use crate::{
    adapter::Mapping,
    aligned_pairs::{AlignOp, TAlignedPairs},
    global_data::GlobalData,
    metrics::{
        TMetric,
        realign::{HP_TR_REG, RealignConfig, realign_hit},
    },
};

#[derive(Debug, Default)]
pub struct AlignCounter {
    pub eq: usize,
    pub diff: usize,
    pub ins: usize,
    pub del: usize,
}

impl AlignCounter {
    pub fn update(&mut self, align_op: AlignOp) {
        match align_op {
            AlignOp::Equal(n) => self.eq += n as usize,
            AlignOp::Diff(n) => self.diff += n as usize,
            AlignOp::Ins(n) => self.ins += n as usize,
            AlignOp::Del(n) => self.del += n as usize,
            _ => panic!("not a valid align op: {:?}", align_op),
        }
    }
}

#[derive(Debug, Default)]
pub struct HpTrMetricV2 {
    qname: Option<String>,
    tname: Option<Arc<String>>,
    global_data: Option<Arc<GlobalData>>,
    align_infos: Vec<Mapping>,
    //Arc<String> for motif string.
    metric_core: HashMap<Arc<String>, AlignCounter>,
    metric_str: Option<String>,
}

impl TMetric for HpTrMetricV2 {
    fn csv_header() -> String {
        let csv_header = [
            "qname".to_string(),
            "rname".to_string(),
            "motif".to_string(),
            "eq".to_string(),
            "diff".to_string(),
            "ins".to_string(),
            "del".to_string(),
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
            shrink_to_boundary: true,
            reference_anchored,
        };
        for old_align_info in &self.align_infos {
            let Some(hit) = realign_hit(read_info, old_align_info, self.get_global_data(), &cfg)
            else {
                continue;
            };

            let mut match_patterns: HashMap<String, Arc<String>> = HashMap::new();
            let region2motif =
                single_seq_hp_tr_finder(&HP_TR_REG, &mut match_patterns, &hit.target_substr)
                    .flatten();

            let align_info = &hit.align;
            let mut pre_ins_cnt = 0;
            let target_start = align_info.target_start() as usize;
            let target_end = align_info.target_end() as usize;

            for (_qpos, rpos, align_op) in align_info.aligned_pairs() {
                if let Some(rpos) = rpos
                    && let Some(motifs) = region2motif.get(&(rpos as usize))
                {
                    motifs
                        .iter()
                        .filter(|motif| motif.0.0 >= target_start && motif.0.1 <= target_end)
                        .for_each(|motif| {
                            self.metric_core
                                .entry(motif.1.clone())
                                .or_default()
                                .update(align_op);

                            self.metric_core
                                .entry(motif.1.clone())
                                .or_default()
                                .update(AlignOp::Ins(pre_ins_cnt));
                        });
                }

                if rpos.is_none() {
                    pre_ins_cnt += 1;
                } else {
                    pre_ins_cnt = 0;
                }
            }
        }
    }

    fn build_metric_str(&mut self) -> String {
        let mut result_items = vec![];
        for (key, counter) in &self.metric_core {
            let mut innner_items = vec![
                self.qname.as_ref().unwrap().clone(),
                self.tname.as_ref().unwrap().as_ref().clone(),
                key.as_ref().clone(),
            ];

            innner_items.push(counter.eq.to_string());
            innner_items.push(counter.diff.to_string());
            innner_items.push(counter.ins.to_string());
            innner_items.push(counter.del.to_string());

            result_items.push(innner_items.join("\t"));
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
        aligner::build_aligner_from_params,
        global_data::GlobalData,
        metrics::{TMetric, hp_tr_metric_v2::HpTrMetricV2},
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
    fn test_hp_tr_metric() {
        let fa_iter = FastaFileReader::new(ref_file());
        let targets = read_fastx(fa_iter);

        let targetname2seq = targets
            .iter()
            .map(|v| (Arc::new(v.name.clone()), Arc::new(v.seq.clone())))
            .collect::<HashMap<_, _>>();

        let target_seq = targetname2seq.values().take(1).collect::<Vec<_>>()[0].clone();

        let global_data = Arc::new(GlobalData::new(targetname2seq));

        let index_params = IndexParams::default().set_kmer(11).set_wins(1);

        let align_params = AlignParams::default()
            .set_m_score(2)
            .set_mm_score(5)
            .set_gap_open_penalty("2,24".to_string())
            .set_gap_extension_penalty("1,0".to_string());
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

        for hit in &hits {
            let (aligned_target, aligned_query) =
                hit.aligned_2_str(target_seq.as_bytes(), fwd_part);
            println!(
                "query_range:{}-{}, target_range:{}-{}, strand:{:?}",
                hit.query_start, hit.query_end, hit.target_start, hit.target_end, hit.strand
            );
            println!("target:{}\nquery :{}", aligned_target, aligned_query);
            println!();
        }

        let mut metric = HpTrMetricV2::default();
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
                "query\ttarget\t(T)3\t6\t0\t0\t0",
                "query\ttarget\t(C)4\t3\t0\t0\t1",
                "query\ttarget\t(A)3\t7\t0\t0\t2",
                "query\ttarget\t(CA)3\t6\t0\t1\t0",
            ],
        );

        let rev_query_record = ReadInfo::new_fa_record(
            "query".to_string(),
            String::from_utf8(rev_part.to_vec()).unwrap(),
        );

        let hits = adapter::align_single_query_to_targets(&aligner, &rev_query_record);
        assert!(!hits.is_empty(), "rev query 应至少产生一个 hit");

        for hit in &hits {
            let (aligned_target, aligned_query) =
                hit.aligned_2_str(target_seq.as_bytes(), rev_part);
            println!(
                "query_range:{}-{}, target_range:{}-{}, strand:{:?}",
                hit.query_start, hit.query_end, hit.target_start, hit.target_end, hit.strand
            );
            println!("target:{}\nquery :{}", aligned_target, aligned_query);
            println!();
        }

        let mut metric = HpTrMetricV2::default();
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
                "query\ttarget\t(C)3\t1\t0\t0\t2",
                "query\ttarget\t(A)3\t2\t0\t0\t1",
                "query\ttarget\t(A)4\t4\t0\t0\t0",
                "query\ttarget\t(T)3\t9\t0\t2\t0",
                "query\ttarget\t(T)5\t5\t0\t0\t0",
            ],
        );
    }
}
