use std::{collections::HashMap, i64, sync::Arc};

use gskits::dna::reverse_complement;
use mm2::minimap2::Mapping;

use crate::{
    aligned_pairs::{AlignOp, TAlignedPairs},
    global_data::{GlobalData, GlobalDataKey, GlobalDataValue},
    metrics::{BASE_COMP, TMetric},
};

#[derive(Debug, Default)]
struct Counter {
    counts: [HashMap<usize, usize>; 2], // [non-misc, misc]
}

impl Counter {
    pub fn update(&mut self, misc: bool, cnt: usize) {
        let idx = if misc { 1 } else { 0 };
        *self.counts[idx].entry(cnt).or_default() += 1;
    }
}

#[derive(Debug, Default)]
pub struct HpTrMetric {
    qname: Option<String>,
    tname: Option<Arc<String>>,
    global_data: Option<Arc<GlobalData>>,
    align_infos: Vec<Mapping>,
    //Arc<String> for motif string. [fwd_counter, rev_counter]
    metric_core: HashMap<Arc<String>, [Counter; 2]>,
    metric_str: Option<String>,
}

impl TMetric for HpTrMetric {
    fn csv_header() -> String {
        let csv_header = vec![
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
    fn compute_metric(&mut self, read_info: &mm2::gskits::ds::ReadInfo) {
        if self.align_infos.is_empty() {
            return;
        }

        let read_seq = read_info.seq.as_bytes();

        let region2motif = self
            .global_data
            .as_ref()
            .unwrap()
            .get(GlobalDataKey::TargetRegion2Motif4Hp);
        let region2motif = match region2motif {
            GlobalDataValue::TargetRegion2Motif4Hp(value) => value,
            _ => panic!(""),
        };
        let tname = self.align_infos[0].target_name.as_ref().unwrap().clone();
        let region2motif = region2motif.get().unwrap().get(&tname).unwrap();

        for align_info in &self.align_infos {
            let strand_idx = if align_info.is_reverse() { 1 } else { 0 };

            for v in
                region2motif.query(align_info.target_start as usize..align_info.target_end as usize)
            {
                let (start, end) = (v.range.start, v.range.end);

                if start < align_info.target_start as usize || end > align_info.query_end as usize {
                    continue;
                }

                let motif = v.value.clone();
                let target_base = if align_info.is_reverse() {
                    BASE_COMP[&motif.as_bytes()[1]]
                } else {
                    motif.as_bytes()[1]
                };

                let mut aligned_pair_iter = Box::new(
                    align_info
                        .aligned_pairs()
                        .take_while(|(_qpos, rpos, _align_op)| {
                            rpos.unwrap_or(i64::MIN) < end as i64
                        }),
                )
                    as Box<dyn Iterator<Item = (Option<i64>, Option<i64>, AlignOp)>>;
                if start > align_info.target_start as usize {
                    aligned_pair_iter = Box::new(
                        aligned_pair_iter
                            .skip_while(|(_qpos, rpos, _align_op)| rpos != &Some(start as i64 - 1)),
                    )
                        as Box<dyn Iterator<Item = (Option<i64>, Option<i64>, AlignOp)>>;
                }

                let mut cnt = 0;
                let mut is_misc = false;
                for (qpos, _rpos, align_op) in aligned_pair_iter {
                    if let Some(qpos) = qpos {
                        cnt += if read_seq[qpos as usize] == target_base {
                            1
                        } else {
                            1
                        };
                        is_misc |= read_seq[qpos as usize] != target_base;
                    }
                }

                self.metric_core.entry(motif.clone()).or_default()[strand_idx].update(is_misc, cnt);
            }
        }
    }

    fn build_metric_str(&mut self) -> String {
        let mut result_items = vec![];
        for (key, fwd_rev_counters) in &self.metric_core {
            for (idx, counter) in fwd_rev_counters.iter().enumerate() {
                for (misc_idx, cnt) in counter.counts.iter().enumerate() {
                    let mut innner_items = vec![];
                    innner_items.push(self.qname.as_ref().unwrap().clone());
                    innner_items.push(self.tname.as_ref().unwrap().as_ref().clone());
                    innner_items.push(if idx == 0 {
                        key.as_ref().clone()
                    } else {
                        reverse_complement_motif(key)
                    });

                    innner_items.push(if misc_idx == 0 {
                        "non-misc".to_string()
                    } else {
                        "misc".to_string()
                    });

                    innner_items.push(counter.eq.to_string());
                    innner_items.push(counter.diff.to_string());
                    innner_items.push(counter.ins.to_string());
                    innner_items.push(counter.del.to_string());

                    result_items.push(innner_items.join("\t"));
                }
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

fn reverse_complement_motif(s: &str) -> String {
    // 用正则提取括号里的序列
    let re = regex::Regex::new(r"\((?P<seq>[ACGT]+)\)").unwrap();

    re.replace_all(s, |caps: &regex::Captures| {
        let seq = &caps["seq"];
        let rc = reverse_complement(seq.as_bytes());
        format!("({})", String::from_utf8(rc).unwrap())
    })
    .into_owned()
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
        metrics::{TMetric, hp_tr_metric::HpTrMetric},
    };

    #[test]
    fn test_hp_tr_metric() {
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

        for hit in &hits {
            let mapping_ext = MappingExt(hit);
            let (aligned_target, aligned_query) =
                mapping_ext.aligned_2_str(target_seq.as_bytes(), fwd_part);
            println!(
                "query_range:{}-{}, target_range:{}-{}, strand:{:?}",
                hit.query_start, hit.query_end, hit.target_start, hit.target_end, hit.strand
            );
            println!("target:{}\nquery :{}", aligned_target, aligned_query);
            println!("");
        }

        let mut metric = HpTrMetric::default();
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

        for hit in &hits {
            let mapping_ext = MappingExt(hit);
            let (aligned_target, aligned_query) =
                mapping_ext.aligned_2_str(target_seq.as_bytes(), rev_part);
            println!(
                "query_range:{}-{}, target_range:{}-{}, strand:{:?}",
                hit.query_start, hit.query_end, hit.target_start, hit.target_end, hit.strand
            );
            println!("target:{}\nquery :{}", aligned_target, aligned_query);
            println!("");
        }

        let mut metric = HpTrMetric::default();
        metric.set_qname("query".to_string());
        metric.set_target_name(Arc::new("target".to_string()));
        metric.set_global_data(global_data.clone());
        metric.set_mappings(hits);
        metric.compute_metric(&rev_query_record);
        metric.set_metric_str();
        println!("metric:\n{}", metric.get_metric_str().as_ref().unwrap());
    }
}
