use std::sync::Arc;

use crate::adapter::Mapping;
use crate::global_data::GlobalData;
use crate::kits::ReadInfo;

pub mod aligned_metric;
pub mod hp_metric_v2;
pub mod hp_tr_metric_v2;
pub mod hp_tr_metric_v3;
pub mod hp_tr_tools;
pub mod realign;

/// 参与 metric 统计的 hit 最低 identity（same-target 过滤后应用）。
/// `Metric::IDENTITY_THRESHOLD` 的默认值来源。
pub const DEFAULT_MIN_IDENTITY: f32 = 0.75;

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

/// 一个 metric mode 的抽象：声明输出列 + 按 ctx 计算。
///
/// 与旧 `TMetric` 的差异：不再要求实现方维护 qname/tname/global_data/mappings/
/// metric_str 五个可变字段与 9 个 accessor；列名在 `COLUMNS` 单点声明，`compute`
/// 产出的行字段顺序必须与之对应。
///
/// `compute` 返回该 read 的输出行（可为空 = 该 read 没有 metric 行）。
/// 同名 hit 过滤、`--noMar` 命中等「该 read 被流水线过滤」的场景也走
/// `compute`（空 mappings），由各 mode 自行决定是否产出占位行
/// （如 `aligned` 的「每条输入 read 恰好一行」语义）。
pub trait Metric: Send + 'static + Default {
    /// 输出列名，顺序与 `compute` 产出的行字段一一对应。
    const COLUMNS: &'static [&'static str];

    /// 参与 metric 统计的 hit 最低 identity（same-target 过滤后应用）；
    /// `None` 表示不做 identity 过滤（默认 `Some(DEFAULT_MIN_IDENTITY)`，
    /// 与迁移前行为一致；`aligned` 等质量评估类 mode 覆盖为 `None`）。
    const IDENTITY_THRESHOLD: Option<f32> = Some(DEFAULT_MIN_IDENTITY);

    fn compute(&mut self, ctx: MetricContext) -> Vec<String>;
}

/// 把行字段按 `COLUMNS` 顺序拼成 TSV 行（列集合的运行期一致性由此保证：
/// 字段按列名顺序 push，缺失字段会立刻在 golden 对比中暴露）。
pub fn join_row(cells: &[String]) -> String {
    cells.join("\t")
}
