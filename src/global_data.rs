use hp_tr_finder::{UnitAndRepeats, all_seq_hp_tr_finder};
use once_cell::sync;
use std::{collections::HashMap, sync::Arc};

#[derive(Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum GlobalDataKey {
    TargetName2Seq,
    TargetRegion2Motif,
}

#[derive(Debug)]
pub enum GlobalDataValue {
    // 第一个Arc<String> 表示的是 target name
    // 第二个 HashMap 的 usize 表示的是 target 的 位置
    //      对应的 value Vec<((usize, usize), Arc<String>)> 表示覆盖该位置的 motifs。为什么是 motifs. 比如说 ACAACAAC[A]AAA. 框起来的位置既是 (ACA)3 覆盖的位置，也是 (A)4 覆盖的位置
    TargetName2Seq(sync::OnceCell<HashMap<Arc<String>, Arc<String>>>),

    TargetRegion2Motif(
        sync::OnceCell<
            HashMap<Arc<String>, Arc<HashMap<usize, Vec<((usize, usize), Arc<String>)>>>>,
        >,
    ),
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
            GlobalDataKey::TargetRegion2Motif,
            GlobalDataValue::TargetRegion2Motif(sync::OnceCell::new()),
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

            GlobalDataValue::TargetRegion2Motif(target_region2motif) => {
                target_region2motif.get_or_init(|| {
                    let all_regs = vec![
                        UnitAndRepeats::new(1, 3).build_finder_regrex(),
                        UnitAndRepeats::new(2, 3).build_finder_regrex(),
                        UnitAndRepeats::new(3, 3).build_finder_regrex(),
                        UnitAndRepeats::new(4, 3).build_finder_regrex(),
                    ];

                    // let targetname2seq = self.target_name2seq.iter().map(|(k, v)| (k.as_str(), v.as_str())).collect();

                    let region2motif: HashMap<
                        Arc<String>,
                        Arc<HashMap<usize, Vec<((usize, usize), Arc<String>)>>>,
                    > = all_seq_hp_tr_finder(&all_regs, &self.target_name2seq)
                        .into_iter()
                        .map(|(seqname, region2motif)| (seqname, Arc::new(region2motif.flatten())))
                        .collect::<HashMap<_, _>>();
                    region2motif
                });
            }
        };

        self.key_values.get(&key).unwrap()
    }
}
