use std::sync::Arc;

use crate::adapter::Mapping;
use crate::global_data::GlobalData;
use crate::kits::ReadInfo;

pub mod hp_metric_v2;
pub mod hp_tr_metric_v2;
pub mod hp_tr_metric_v3;
pub mod hp_tr_tools;
pub mod realign;

/// 一次性携带 `Metric::compute` 所需的全部上下文，取代原先散落在 trait
/// setter 里的可变状态（qname/tname/global_data/mappings）。
pub struct MetricContext {
    /// 当前 read（所有权随 ctx 移交）。
    pub read: ReadInfo,
    /// 目标名（来自第一个 primary hit），输出行的 rname 列取此值。
    pub target_name: Arc<String>,
    /// 该 read 的 hits，已过 main 层的 primary/supp、same-target、identity 过滤。
    pub mappings: Vec<Mapping>,
    pub global_data: Arc<GlobalData>,
    pub reference_anchored: bool,
}

/// `Metric::compute` 的显式结果，取代原先用「空 metric_str」当丢行哨兵的
/// 隐式协议。
pub enum MetricOutcome {
    /// 该 read 的输出行（可为空 = 该 read 没有 metric 行）。
    Rows(Vec<String>),
    /// 该 read 被过滤（如 --noMar 命中多 chain），不产出任何行。
    Skipped,
}

/// 一个 metric mode 的抽象：声明输出列 + 按 ctx 计算。
///
/// 与旧 `TMetric` 的差异：不再要求实现方维护 qname/tname/global_data/mappings/
/// metric_str 五个可变字段与 9 个 accessor；「跳过该 read」由 `Skipped` 显式
/// 表达；列名在 `COLUMNS` 单点声明，`compute` 产出的行字段顺序必须与之对应。
pub trait Metric: Send + 'static + Default {
    /// 输出列名，顺序与 `compute` 产出的行字段一一对应。
    const COLUMNS: &'static [&'static str];

    fn compute(&mut self, ctx: MetricContext) -> MetricOutcome;
}

/// 把行字段按 `COLUMNS` 顺序拼成 TSV 行（列集合的运行期一致性由此保证：
/// 字段按列名顺序 push，缺失字段会立刻在 golden 对比中暴露）。
pub fn join_row(cells: &[String]) -> String {
    cells.join("\t")
}
