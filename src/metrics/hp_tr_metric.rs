use std::{collections::HashMap, sync::Arc};

use gskits::dna::reverse_complement;
use mm2::minimap2::Mapping;

use crate::{
    aligned_pairs::{AlignOp, TAlignedPairs},
    global_data::{GlobalData, GlobalDataKey, GlobalDataValue},
    metrics::TMetric,
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
pub struct HpTrMetric {
    qname: Option<String>,
    tname: Option<Arc<String>>,
    global_data: Option<Arc<GlobalData>>,
    align_infos: Vec<Mapping>,
    //Arc<String> for motif string. [fwd_counter, rev_counter]
    metric_core: HashMap<Arc<String>, [AlignCounter; 2]>,
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
    fn compute_metric(&mut self, _read_info: &mm2::gskits::ds::ReadInfo) {
        if self.align_infos.is_empty() {
            return;
        }

        let region2motif = self
            .global_data
            .as_ref()
            .unwrap()
            .get(GlobalDataKey::TargetRegion2Motif);
        let region2motif = match region2motif {
            GlobalDataValue::TargetRegion2Motif(value) => value,
            _ => panic!(""),
        };
        let tname = self.align_infos[0].target_name.as_ref().unwrap().clone();
        let region2motif = region2motif.get().unwrap().get(&tname).unwrap();

        for align_info in &self.align_infos {
            let mut pre_ins_cnt = 0;
            let target_start = align_info.target_start() as usize;
            let target_end = align_info.target_end() as usize;
            let strand_idx = if align_info.is_reverse() { 1 } else { 0 };

            for (_qpos, rpos, align_op) in align_info.aligned_pairs() {
                if let Some(rpos) = rpos {
                    if let Some(motifs) = region2motif.get(&(rpos as usize)) {
                        motifs
                            .iter()
                            .filter(|motif| motif.0.0 >= target_start && motif.0.1 <= target_end)
                            .for_each(|motif| {
                                self.metric_core.entry(motif.1.clone()).or_default()[strand_idx]
                                    .update(align_op);

                                self.metric_core.entry(motif.1.clone()).or_default()[strand_idx]
                                    .update(AlignOp::Ins(pre_ins_cnt));
                            });
                    }
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
        for (key, fwd_rev_counters) in &self.metric_core {
            for (idx, counter) in fwd_rev_counters.iter().enumerate() {
                let mut innner_items = vec![];
                innner_items.push(self.qname.as_ref().unwrap().clone());
                innner_items.push(self.tname.as_ref().unwrap().as_ref().clone());
                innner_items.push(if idx == 0 {
                    key.as_ref().clone()
                } else {
                    reverse_complement_motif(key)
                });

                innner_items.push(counter.eq.to_string());
                innner_items.push(counter.diff.to_string());
                innner_items.push(counter.ins.to_string());
                innner_items.push(counter.del.to_string());

                result_items.push(innner_items.join("\t"));
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
