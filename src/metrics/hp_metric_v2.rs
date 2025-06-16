use std::{collections::HashMap, i64, sync::Arc};

use hp_tr_finder::{UnitAndRepeats, single_seq_hp_tr_finder};
use mm2::{mapping_ext::MappingExt, minimap2::Mapping};
use regex::Regex;

use crate::{
    aligned_pairs::{AlignOp, TAlignedPairs},
    global_data::{GlobalData, GlobalDataKey, GlobalDataValue},
    metrics::{
        TMetric,
        hp_tr_tools::{do_align_4_homo, get_target_substr},
    },
};

#[derive(Debug, Default)]
struct Counter {
    counts: [HashMap<usize, usize>; 2], // [pure, mixed]
}

impl Counter {
    pub fn update(&mut self, misc: bool, cnt: usize) {
        let idx = if misc { 1 } else { 0 };
        *self.counts[idx].entry(cnt).or_default() += 1;
    }
}

lazy_static::lazy_static! {
    static ref HP_REG: Vec<HashMap<String, Regex>> = {
        vec![UnitAndRepeats::new(1, 3).build_finder_regrex()]
    };
}

#[derive(Debug, Default)]
pub struct HpMetricV2 {
    qname: Option<String>,
    tname: Option<Arc<String>>,
    global_data: Option<Arc<GlobalData>>,
    align_infos: Vec<Mapping>,
    //Arc<String> for motif string. [fwd_counter, rev_counter]
    metric_core: HashMap<Arc<String>, Counter>,
    metric_str: Option<String>,
}

impl TMetric for HpMetricV2 {
    fn csv_header() -> String {
        let csv_header = vec![
            "qname".to_string(),
            "rname".to_string(),
            "motif".to_string(),
            "tag".to_string(),
            "called".to_string(),
            "num".to_string(),
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
    fn compute_metric(&mut self, read_info: &mm2::gskits::ds::ReadInfo) {
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

        // let read_seq = read_info.seq.as_bytes();

        for old_align_info in &self.align_infos {
            let read_seq = &read_info.seq.as_bytes()
                [old_align_info.query_start as usize..old_align_info.query_end as usize];

            let target_substr = get_target_substr(
                old_align_info.target_start as usize,
                old_align_info.target_end as usize,
                old_align_info.is_reverse(),
                target_seq_fwd,
                target_seq_rev,
            );

            let align_info = do_align_4_homo(
                &read_info.seq,
                old_align_info.query_start as usize,
                old_align_info.query_end as usize,
                target_substr,
            );
            if align_info.is_none() {
                continue;
            }
            let align_info = align_info.unwrap();

            assert!(matches!(align_info.strand, mm2::minimap2::Strand::Forward));
            let mut match_patterns: HashMap<String, Arc<String>> = HashMap::new();
            let region2motif = single_seq_hp_tr_finder(&HP_REG, &mut match_patterns, target_substr);

            for v in region2motif.iter() {
                let (start, end) = (v.0.0, v.0.1);
                if start < align_info.target_start as usize || end > align_info.target_end as usize
                {
                    continue;
                }

                let motif = v.1.clone();
                let target_base = motif.as_bytes()[1];
                let expected_cnt = motif[3..].parse::<i32>().unwrap();
                let mut pre_rpos = -1;
                let mut aligned_pair_iter = Box::new(align_info.aligned_pairs().take_while(
                    |(_qpos, rpos, _align_op)| {
                        let pre_rpos_ = pre_rpos;
                        if let Some(rpos) = rpos {
                            pre_rpos = *rpos;
                        }
                        (pre_rpos_ + 1) < end as i64
                    },
                ))
                    as Box<dyn Iterator<Item = (Option<i64>, Option<i64>, AlignOp)>>;
                if start > align_info.target_start as usize {
                    aligned_pair_iter = Box::new(
                        aligned_pair_iter
                            .skip_while(|(_qpos, rpos, _align_op)| rpos != &Some(start as i64 - 1))
                            .skip(1),
                    )
                        as Box<dyn Iterator<Item = (Option<i64>, Option<i64>, AlignOp)>>;
                }

                let mut cnt = 0;
                let mut is_misc = false;
                
                for (qpos, _rpos, _align_op) in aligned_pair_iter {
                    if let Some(qpos) = qpos {
                        cnt += if read_seq[qpos as usize] == target_base {
                            1
                        } else {
                            0
                        };
                        is_misc |= read_seq[qpos as usize] != target_base;
                    }
                }

                if (cnt - expected_cnt).abs() >= 10 {
                    let (aligned_target, aligned_query) =
                        MappingExt(&align_info).aligned_2_str(target_substr.as_bytes(), read_seq);
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
                    .update(is_misc, cnt as usize);
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
        mapping_ext::MappingExt,
        params::{AlignParams, IndexParams, MapParams, OupParams},
    };

    use crate::{
        global_data::GlobalData,
        metrics::{TMetric, hp_metric_v2::HpMetricV2},
    };

    #[test]
    fn test_hp_metric() {
        let ref_file = "test_data/ref_Saureus_ATCC25923.m.new.corrected.fasta";
        let fa_iter = FastaFileReader::new(ref_file.to_string());
        let targets = read_fastx(fa_iter);

        let targetname2seq = targets
            .iter()
            .map(|v| (Arc::new(v.name.clone()), Arc::new(v.seq.clone())))
            .collect::<HashMap<_, _>>();

        let target_seq = targetname2seq
            .values()
            .take(1)
            .map(|v| v)
            .collect::<Vec<_>>()[0]
            .clone();

        let global_data = Arc::new(GlobalData::new(targetname2seq));

        let mut index_params = IndexParams::default();
        index_params.kmer = Some(11);
        index_params.wins = Some(1);

        let align_params = AlignParams::default();
        let mut aligners = build_aligner(
            "map-ont",
            &index_params,
            &MapParams::default(),
            &align_params,
            &OupParams::default(),
            &targets,
            10,
        );

        // // aligner.mapopt.min_cnt = 2;
        // aligner.mapopt.best_n = 100000;
        // // aligner.mapopt.min_dp_max = 10; // min dp score
        // aligner.mapopt.pri_ratio = 0.5; // min dp score
        // aligner.idxopt;
        println!("{:?}", aligners[0].mapopt);
        println!("{:?}", aligners[0].idxopt);

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

        // for hit in &hits {
        //     let mapping_ext = MappingExt(hit);
        //     let (aligned_target, aligned_query) =
        //         mapping_ext.aligned_2_str(target_seq.as_bytes(), fwd_part);
        //     println!(
        //         "query_range:{}-{}, target_range:{}-{}, strand:{:?}",
        //         hit.query_start, hit.query_end, hit.target_start, hit.target_end, hit.strand
        //     );
        //     println!("target:{}\nquery :{}", aligned_target, aligned_query);
        //     println!("");
        // }

        let mut metric = HpMetricV2::default();
        metric.set_qname("query".to_string());
        metric.set_target_name(Arc::new("target".to_string()));
        metric.set_global_data(global_data.clone());
        metric.set_mappings(hits);
        metric.compute_metric(&fwd_query_record);
        metric.set_metric_str();
        println!("metric:\n{}", metric.get_metric_str().as_ref().unwrap());

        let rev_query_record = ReadInfo::new_fa_record(
            "query".to_string(),
            String::from_utf8(rev_part.to_vec()).unwrap(),
        );

        let hits = align_single_query_to_targets(&rev_query_record, &aligners);

        // for hit in &hits {
        //     let mapping_ext = MappingExt(hit);
        //     let (aligned_target, aligned_query) =
        //         mapping_ext.aligned_2_str(target_seq.as_bytes(), rev_part);
        //     println!(
        //         "query_range:{}-{}, target_range:{}-{}, strand:{:?}",
        //         hit.query_start, hit.query_end, hit.target_start, hit.target_end, hit.strand
        //     );
        //     println!("target:{}\nquery :{}", aligned_target, aligned_query);
        //     println!("");
        // }

        let mut metric = HpMetricV2::default();
        metric.set_qname("query".to_string());
        metric.set_target_name(Arc::new("target".to_string()));
        metric.set_global_data(global_data.clone());
        metric.set_mappings(hits);
        metric.compute_metric(&rev_query_record);
        metric.set_metric_str();
        println!("metric:\n{}", metric.get_metric_str().as_ref().unwrap());
    }

    #[test]
    fn test_hp_tr_metric_case() {
        let ref_file = "test_data/ref_Saureus_ATCC25923.m.new.corrected.fasta";
        let fa_iter = FastaFileReader::new(ref_file.to_string());
        let targets = read_fastx(fa_iter);

        let targetname2seq = targets
            .iter()
            .map(|v| (Arc::new(v.name.clone()), Arc::new(v.seq.clone())))
            .collect::<HashMap<_, _>>();

        let target_seq = targetname2seq
            .values()
            .take(1)
            .map(|v| v)
            .collect::<Vec<_>>()[0]
            .clone();

        let global_data = Arc::new(GlobalData::new(targetname2seq));

        let mut index_params = IndexParams::default();
        index_params.kmer = Some(11);
        index_params.wins = Some(1);

        let align_params = AlignParams::default();
        let mut aligners = build_aligner(
            "map-ont",
            &index_params,
            &MapParams::default(),
            &align_params,
            &OupParams::default(),
            &targets,
            10,
        );

        // // aligner.mapopt.min_cnt = 2;
        // aligner.mapopt.best_n = 100000;
        // // aligner.mapopt.min_dp_max = 10; // min dp score
        // aligner.mapopt.pri_ratio = 0.5; // min dp score
        // aligner.idxopt;
        println!("{:?}", aligners[0].mapopt);
        println!("{:?}", aligners[0].idxopt);

        aligners[0].mapopt.best_n = 100000;
        aligners[0].mapopt.pri_ratio = 0.1;
        let seq = b"AAGGGGGGGGGGGGGGGAAAAAGGGGGGGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAATAAAAAAAAAAAAAAAAAAAAAAGGAAAAAAAGGAAAAAAAAAAAAAAAAAAAAGGGGAAAAAAAAAAGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGGGGGGGGGGGAAAAAAAAAAAAATTTTCCCCCCCAAAAAAAAAAAAACCCCCCCCCCCCCCTTGGGGGGGCGGGGGGGTTTTTTGGGGCATTGGGGGGGGGGGTGGGGGGGGGGGTTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATAAAAAATAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCAAAAAAAAAAAAAATAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTATTGGGGGGGGGTTTTTTTTTTTTTTTTTTTGGGGGGGGGGGGGCCCCCCCCCAAAATCCCCCCCCCCCTTTATCCCCCCCCCAAAAAAAAAAAAAATTTCCCCCCCCCCCCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTTTTTTTCCCCCCCCCCCCCCCCCCCCCCCCAAAAAAACAACCCCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATATAAAAAAAAAAAAAAAAAAATTTTTTTTTGTGGGGGGGGGGGGGGGGGGGGGGGTGGGGGTTTTTCCCCCCCCCCCCCCAAAAATAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGGCCCCCCCCCCCTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATGTTTTTTGTTTAGGGGGGGGGGGGCCGCCCCGTTTTATTTTTTTTTTTTTAACACCAATTTTTTTTTATTTTTGGTTGGGTTTTCACATGAAAGCGGGCCGGGGGCAAAATTAAATTCGGGGGGGGGGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAATTTTTTTTTTTTTTAAAAAAAAAAAAAAACAAAACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGGGGGGGGGGGCTTTATTTTTTTTTTTTATAAAAAGGGGCTTTCTGTGAAATGAAGAACAAACAAAAATCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACACTTCTAAAAAAAAAACCCCCCCCCCCCAAAAAGGCCGGGGGGAAAAAAAAAAAAACCTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGGGCCGGAAAAAAAAAAAAAAAAAAAAAGGGGGCCGAATTCCCCCCTGCCGAAAAAAACAAAACGAAAAAGCAAAAAGAAAAAAACAAAAAAAACAACAAAAAACAAAAAAAAGCCCCTGCCCCCCCCCCCCCCTTTAAAACCACAAACAATTAGTTTTTAAACGAGGGGGTTTTTTTTTTTTTTTTTTTTTTTCGTAAACGCGTTATATTTTTTTTTTTTTACAACCTGTTAACCATTTTTTCTTAACGTCCTACAACATTTTTTTTCCCTCTTTTTTTAAAAAATAAAAAAAAAAAAAAAAAAAGAGGTCAAAAAAAAACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGCTGCTTTCATTTTTTGGGGCAAAAAAAAAAAGAAAACAAAAAAAAAAAAAAAACAAAAAAAAAACAAAAAAAAAAAAACAACCAACACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACAAAAAAAACCAACACCCCGACAAAAAAACCCCCCCCAAAAAAAAAAAATTTTTTGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGGGCCCCCAAAAAAAAAAAAATTAGGGGGGGGGGGGGGGTGGTTTTTGGGGGGGTGGTGGGTTGGTTTTTTTGGGGGGGGGGGGGCGGTTTTTTTTTTTTTCCCCGGGGGCCCCCCCCCTTTTTTTTTTTTTTTTTTTTTTTTATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTAAAAAAAAAAAAAAAAAAATTTTTTTTTTTCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTAATATTTTTTTTCCCCCCCCCCCCCCCCCCCCGGGGGGAAAGGGGGGCGGTTAAAAAATTTTTTTTTTTTTTTTTTTTTTTTATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTAATTTGGGGGGGGGGGGGGGGGGGAGGGGGGTTTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAAAAGGGGGGGGGTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTTTTTTTTTTTATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGGGGTTTATTTTTTTTTTTTTAATTTTTTTTTTTTTTTTTTTTTTTTTTATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCCGGGGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTTTTTTTTTTAAAAAAAAAAAATTCCCCCCCCCCCCGGGGGGGGGGGGGGAAAATTTTTTTTTTTTTTTTTTTTTATTTTTTTTTTTTTTTTTTTTAAAAAAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGAAAAAAAAAATTTTTTTTTTTTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTTTTTTTTTTTTTTTTTTTTAAAAAAAAAAAAAGGGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTTTTATTTTTTTTTTTTTTTTTTTTTTTTTTTTTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGCCGGCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTGGGGGGGGGGGGGGCCTTTTTTTTTTTTTTTTTTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCAAAAAAAAAAAATTTTTTTTTTGTTGTGGTGGGGGGGGGGGGGGGGGGGGGGGAAAGGGGGGGGGGGGGGGGGGGGGGAAAAAAAAAAAAAAAAAAAAGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGGGGGGGGGTTTTTTATTTTTTTTAAAAAAAAAAAAAAAAAAAAAAAAAAAAATAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGGGGGAAAAAAAACCCCCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTATTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCGCCCCCCCCCCATAGGAGGGGGGGCCCCCGTTTATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTAAAAAAAAAAAAAAAAAAAAAAAGTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGCCCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATGCCCCCAAAAAAAAAAAAAAAATTAAAAAAAAAAAAAAAAAAAAAAAAAAGGGTTTTAAAAAAAAAAAAAAAAAAAAAAAATAAAAAAAAAAAAAAAAAAAAAAACAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGTCCCCCCCCCCCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTAAAAAAAAAAAAAAAAAAAAAGGGGGGGGGCCCCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACAAAAAAAAAAAAAAAAGGGGGGGGGCCCCGCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGGGGGGGCCCAAAAAAAAGGGGGGCCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGTTTTTTTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGGGGGGGGGGGGCCCCAAAAAAAAAAAAAAAAATTTTTTTTTTTTTTTTTTTTTTTTATTTTTTTGGGGGGGGGAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCACCCCCCCAAAAAAAAAAAATTTTTTTCCCAAAAAAAAAAAAAAAAGGGGGGGGAAATCCCCCCCCCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAGGGGGGGGGAAAAAATTCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGGAAAAAAAAAAAAAAAAAAAAAAAAAAGGGGAAAAAAAAAAAAAAAAAAAGGGAAGGGGGGGGGAAATTTCCCCCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        let fwd_part = &seq[..];
        let rev_part = &seq[1074..1074 + 100];
        let fwd_query_record = ReadInfo::new_fa_record(
            "query".to_string(),
            String::from_utf8(fwd_part.to_vec()).unwrap(),
        );

        let hits = align_single_query_to_targets(&fwd_query_record, &aligners);

        for hit in &hits {
            let mapping_ext = MappingExt(hit);
            let (aligned_target, aligned_query) =
                mapping_ext.aligned_2_str(target_seq.as_bytes(), fwd_part);
            println!(
                "query_range:{}-{}, target_range:{}-{}, strand:{:?}, identity:{}",
                hit.query_start,
                hit.query_end,
                hit.target_start,
                hit.target_end,
                hit.strand,
                mapping_ext.identity()
            );
            println!("target:{}\nquery :{}", aligned_target, aligned_query);
            println!("");
        }

        let mut metric = HpMetricV2::default();
        metric.set_qname("query".to_string());
        metric.set_target_name(Arc::new("target".to_string()));
        metric.set_global_data(global_data.clone());
        metric.set_mappings(hits);
        metric.compute_metric(&fwd_query_record);
        metric.set_metric_str();
        println!("metric:\n{}", metric.get_metric_str().as_ref().unwrap());

        // let rev_query_record = ReadInfo::new_fa_record(
        //     "query".to_string(),
        //     String::from_utf8(rev_part.to_vec()).unwrap(),
        // );

        // let hits = align_single_query_to_targets(&rev_query_record, &aligners);

        // for hit in &hits {
        //     let mapping_ext = MappingExt(hit);
        //     let (aligned_target, aligned_query) =
        //         mapping_ext.aligned_2_str(target_seq.as_bytes(), rev_part);
        //     println!(
        //         "query_range:{}-{}, target_range:{}-{}, strand:{:?}",
        //         hit.query_start, hit.query_end, hit.target_start, hit.target_end, hit.strand
        //     );
        //     println!("target:{}\nquery :{}", aligned_target, aligned_query);
        //     println!("");
        // }

        // let mut metric = HpMetric::default();
        // metric.set_qname("query".to_string());
        // metric.set_target_name(Arc::new("target".to_string()));
        // metric.set_global_data(global_data.clone());
        // metric.set_mappings(hits);
        // metric.compute_metric(&rev_query_record);
        // metric.set_metric_str();
        // println!("metric:\n{}", metric.get_metric_str().as_ref().unwrap());
    }
}
