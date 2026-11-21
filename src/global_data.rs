use std::{collections::HashMap, sync::Arc};

use hp_tr_finder::all_seq_hp_tr_finder;
use once_cell::sync::OnceCell;

use crate::{kits::reverse_complement, metrics::realign::HP_TR_REG};

/// ref 位点 -> 覆盖该位点的 motif 区间（HP-TR 用 flatten 后的 HashMap 索引）。
pub type TargetRegion2Motif4HpTrMap =
    HashMap<Arc<String>, Arc<HashMap<usize, Vec<((usize, usize), Arc<String>)>>>>;

/// 按 target 惰性初始化的全局共享数据，metric worker 线程通过 `Arc` 共享。
/// 类型即文档：新增数据源加一个 `OnceCell<T>` 字段与对应 getter 即可，
/// 不再需要 enum 分派。
pub struct GlobalData {
    target_name2seq: HashMap<Arc<String>, Arc<String>>,
    /// target 名 -> [正向序列, 反向互补序列]。
    seqs_with_rev: OnceCell<HashMap<Arc<String>, [Arc<String>; 2]>>,
    /// target 名 -> 该 target 的 HP-TR motif 索引（由 HP_TR_REG 全序列扫描构建）。
    region2motif_hptr: OnceCell<TargetRegion2Motif4HpTrMap>,
}

impl GlobalData {
    pub fn new(target_name2seq: HashMap<Arc<String>, Arc<String>>) -> Self {
        Self {
            target_name2seq,
            seqs_with_rev: OnceCell::new(),
            region2motif_hptr: OnceCell::new(),
        }
    }

    pub fn seqs_with_rev(&self) -> &HashMap<Arc<String>, [Arc<String>; 2]> {
        self.seqs_with_rev.get_or_init(|| {
            self.target_name2seq
                .iter()
                .map(|(name, seq)| {
                    (
                        name.clone(),
                        [
                            seq.clone(),
                            Arc::new(
                                String::from_utf8(reverse_complement(seq.as_bytes())).unwrap(),
                            ),
                        ],
                    )
                })
                .collect()
        })
    }

    pub fn region2motif_hptr(&self) -> &TargetRegion2Motif4HpTrMap {
        self.region2motif_hptr.get_or_init(|| {
            all_seq_hp_tr_finder(&HP_TR_REG, &self.target_name2seq)
                .into_iter()
                .map(|(seqname, region2motif)| (seqname, Arc::new(region2motif.flatten())))
                .collect::<HashMap<_, _>>()
        })
    }
}
