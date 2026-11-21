use crate::kits::reverse_complement;
use hp_tr_finder::{UnitAndRepeats, all_seq_hp_tr_finder};
use once_cell::sync;
use std::{collections::HashMap, sync::Arc};

use crate::metrics::hp_tr_metric_v2::HP_TR_REG;

/// ref 位点 -> 覆盖该位点的 motif 区间（HP-TR 用 flatten 后的 HashMap 索引）。
pub type TargetRegion2Motif4HpTrMap =
    HashMap<Arc<String>, Arc<HashMap<usize, Vec<((usize, usize), Arc<String>)>>>>;
/// HP 用 IntervalTree 索引。
pub type TargetRegion2Motif4HpMap =
    HashMap<Arc<String>, Arc<hp_tr_finder::intervaltree::IntervalTree<usize, Arc<String>>>>;

#[derive(Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum GlobalDataKey {
    TargetName2Seq,
    TargetName2SeqAndRev,

    TargetRegion2Motif4HpTr,
    TargetRegion2Motif4Hp,
}

#[derive(Debug)]
pub enum GlobalDataValue {
    // 第一个Arc<String> 表示的是 target name
    // 第二个 HashMap 的 usize 表示的是 target 的 位置
    //      对应的 value Vec<((usize, usize), Arc<String>)> 表示覆盖该位置的 motifs。为什么是 motifs. 比如说 ACAACAAC[A]AAA. 框起来的位置既是 (ACA)3 覆盖的位置，也是 (A)4 覆盖的位置
    TargetName2Seq(sync::OnceCell<HashMap<Arc<String>, Arc<String>>>),
    TargetName2SeqAndRev(sync::OnceCell<HashMap<Arc<String>, [Arc<String>; 2]>>),

    TargetRegion2Motif4HpTr(sync::OnceCell<TargetRegion2Motif4HpTrMap>),

    TargetRegion2Motif4Hp(sync::OnceCell<TargetRegion2Motif4HpMap>),
}

#[derive(Debug)]
pub struct GlobalData {
    target_name2seq: HashMap<Arc<String>, Arc<String>>,

    key_values: HashMap<GlobalDataKey, GlobalDataValue>,
}

impl GlobalData {
    pub fn new(target_name2seq: HashMap<Arc<String>, Arc<String>>) -> Self {
        let mut key_values = HashMap::new();
        key_values.insert(
            GlobalDataKey::TargetName2Seq,
            GlobalDataValue::TargetName2Seq(sync::OnceCell::new()),
        );

        key_values.insert(
            GlobalDataKey::TargetName2SeqAndRev,
            GlobalDataValue::TargetName2SeqAndRev(sync::OnceCell::new()),
        );

        key_values.insert(
            GlobalDataKey::TargetRegion2Motif4HpTr,
            GlobalDataValue::TargetRegion2Motif4HpTr(sync::OnceCell::new()),
        );

        key_values.insert(
            GlobalDataKey::TargetRegion2Motif4Hp,
            GlobalDataValue::TargetRegion2Motif4Hp(sync::OnceCell::new()),
        );

        Self {
            target_name2seq,
            key_values,
        }
    }

    pub fn get(&self, key: GlobalDataKey) -> &GlobalDataValue {
        match self.key_values.get(&key).unwrap() {
            GlobalDataValue::TargetName2Seq(target_name2seq) => {
                target_name2seq.get_or_init(|| self.target_name2seq.clone());
            }

            GlobalDataValue::TargetName2SeqAndRev(target_name2seq) => {
                target_name2seq.get_or_init(|| {
                    self.target_name2seq
                        .iter()
                        .map(|(name, seq)| {
                            (
                                name.clone(),
                                [
                                    seq.clone(),
                                    Arc::new(
                                        String::from_utf8(reverse_complement(seq.as_bytes()))
                                            .unwrap(),
                                    ),
                                ],
                            )
                        })
                        .collect()
                });
            }

            GlobalDataValue::TargetRegion2Motif4HpTr(target_region2motif) => {
                target_region2motif.get_or_init(|| {
                    let region2motif: TargetRegion2Motif4HpTrMap =
                        all_seq_hp_tr_finder(&HP_TR_REG, &self.target_name2seq)
                            .into_iter()
                            .map(|(seqname, region2motif)| {
                                (seqname, Arc::new(region2motif.flatten()))
                            })
                            .collect::<HashMap<_, _>>();
                    region2motif
                });
            }

            GlobalDataValue::TargetRegion2Motif4Hp(target_region2motif) => {
                target_region2motif.get_or_init(|| {
                    let all_regs = vec![UnitAndRepeats::new(1, 3).build_finder_regrex()];
                    let region2motif: TargetRegion2Motif4HpMap =
                        all_seq_hp_tr_finder(&all_regs, &self.target_name2seq)
                            .into_iter()
                            .map(|(seqname, region2motif)| {
                                (seqname, Arc::new(region2motif.to_interval_search_tree()))
                            })
                            .collect::<HashMap<_, _>>();
                    region2motif
                });
            }
        };

        self.key_values.get(&key).unwrap()
    }
}
