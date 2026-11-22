# mm2 metric 功能独立化：overlap 差异与迁移缺口报告

> 对照基准：mm2 @ `339ab46`（main，2026-10-15，v1.2.1） vs gsmm2-metric @ `a37faff`（v0.10.3）。
> 方法：逐文件通读 mm2 的 5 个 bin 与 `align_processor/*`、`lib.rs`、`params.rs`、`aligned_pairs.rs`，
> 与本项目 `src/main.rs`、`src/cli.rs`、`src/params.rs`、`src/aligner.rs`、`src/metrics/*`、`src/kits/*` 对照。
> 结论先行：**overlap 的 hp-tr 算子谱 metric 已迁移并有意演化（golden 锁定），aligned metric 与 time-err 两个 bin 完全未迁移**；
> overlap 部分还发现 1 个疑似回归（fastq `--rq-range` 失效）和若干口径差异需要决策。

---

## 一、功能迁移状态总览

| mm2 侧（bin + 核心模块） | 功能 | gsmm2-metric 对应物 | 状态 |
|---|---|---|---|
| `gsmm2-hp-tr-metric.rs`（646 行）+ `align_processor/hp_tr_metric.rs`（156 行） | HP-TR motif 比对算子谱：`qname rname motif eq diff ins del` | mode `hp-tr-v2` | **已迁移并大幅演化**，详见 §三 |
| `gsmm2-aligned-metric.rs`（616 行）+ `align_processor/align_metric.rs`（498 行） | 全比对质量 metric：segs/covlen/identity/homo-ins-del/long-indel/per-base 分解（21+24 列） | mode `aligned` | **✅ 已迁移**（2026-11-22，golden 与 mm2 原版二进制逐字节一致），详见 §三 |
| `gsmm2-time-err.rs`（610 行）+ `align_processor/time_err_metric.rs`（200 行） | 按测序时间分桶（600s）的错误构成 metric，依赖 BAM 的 `ar`/`dw` tag | — | **❌ 未迁移**，详见 §4.2 |
| `gsmm2.rs`（407 行） | aligner 主 bin（比对→BAM，`--hpTrShrinkAlnRegion` / `--polyNGapLeftAlign` / `--query-forward`） | — | 非 metric 功能，不在独立范围（§五） |
| `bam-rw.rs`（92 行） | BAM 拷贝小工具 | — | 与 metric 无关 |
| `visualization.rs`（211 行） | blast 风格比对可视化打印（调试辅助） | — | 未迁移（可选，§4.3） |
| （无对应） | — | mode `hp-tr-v3` / `hp-v2` / `bam2-msa` | 本项目新增扩展 |

本项目相对 mm2 metric bin 的**结构性增益**已经完成且被 golden 锁定：比对引擎从 mm2 FFI 换为 bandmap-core（见
[mm2_bandmap_migration.md](mm2_bandmap_migration.md)，真实数据行级差异 ≤0.0025%）；三 mode 共用重对齐前缀
（[realign.rs](../src/metrics/realign.rs)）；rust-htslib/gskits → noodles + `src/kits`。

---

## 二、overlap 差异明细：hp-tr 算子谱 metric（mm2 `gsmm2-hp-tr-metric` ↔ `hp-tr-v2`）

### 2.1 管线差异（有意演化，golden 已锁定——需要与下游确认「权威口径」）

| 环节 | mm2 `gsmm2-hp-tr-metric` | 本项目 `hp-tr-v2` |
|---|---|---|
| 比对结果来源 | 直接消费 mm2 raw mapping 的 aligned pairs（[gsmm2-hp-tr-metric.rs:452-509](../../../mm2/src/bin/gsmm2-hp-tr-metric.rs)） | 每个 hit 先做边界收缩 + target 子串 + `do_align_4_homo` 精重对齐（[realign.rs](../src/metrics/realign.rs)） |
| motif 索引 | `all_seq_hp_tr_finder` 对**全 target** 预计算一次（bin 启动时） | 在重对齐后的 target 子串上现算 `single_seq_hp_tr_finder`；全 target 索引仅用于边界收缩（[global_data.rs](../src/global_data.rs)） |
| strand 处理 | 每个 motif 输出 `[fwd, rev]` **两行**，rev 行 motif 名反向互补；未被触碰的计数器也输出全 0 行（[hp_tr_metric.rs:69-93](../../../mm2/src/align_processor/hp_tr_metric.rs)） | 每个 motif **一行**，反向 hit 的统计在 read 方向的 RC 子串上完成，motif 名即 read 方向（`--refAnchored` 可切换为 reference 方向） |
| hit 过滤 | primary/supp 保留 + same-target，**无 identity 过滤** | 同前两项 + `identity > 0.75` 硬过滤（[main.rs:52](../src/main.rs)），0.10.2 起 `--noMar` 还可整条过滤多 chain read |
| 边界收缩 / refAnchored | 无 | 有（`shrink_to_motif_boundary` + `get_target_substr`，[hp_tr_tools.rs](../src/metrics/hp_tr_tools.rs)） |
| 输入限制 | `assert` 全部 query 必须是 BAM（[gsmm2-hp-tr-metric.rs:351-360](../../../mm2/src/bin/gsmm2-hp-tr-metric.rs)） | fa/fasta/fna/fq/fastq/bam 均可 |

这些差异是本项目演化方向（重对齐管线、golden 验收）的一部分，**不是回归**；但若期望本 bin 是 mm2
`gsmm2-hp-tr-metric` 的「直接平替」，需与下游确认：

1. **rev strand 的全 0 行不再输出**——mm2 每 motif 恒 2 行（其中 rev 行常为全 0），本项目 1 行。下游若有按
   strand 拆列的逻辑会受影响。
2. **unmapped / 被过滤 read 的占位行不再输出**——mm2 对 unmapped（或 `--noMar` 命中）read 仍写一行
   `qname\t\t0\t0\t0\t0\t0` 占位（复用 `TimeErrMetric` 的空兜底 Display，[gsmm2-hp-tr-metric.rs:463-470](../../../mm2/src/bin/gsmm2-hp-tr-metric.rs)）；
   映上了但无 motif 的 read 写**空行**。本项目 `MetricOutcome::Rows(vec![])` / `Skipped` 均不产出任何行。
   下游按行数/占位行对账时需注意。
3. **默认输出文件名不同**：mm2 `{q}.gsmm2-hp-tr-fact.csv` vs 本项目 `{q}.HpTrV2-fact.csv`。

### 2.2 参数与依赖口径差异

| 项 | mm2 | 本项目 | 处理建议 |
|---|---|---|---|
| `chaining.max_gap` | hp-tr / time-err bin 设 `mapopt.max_gap = 100`（[gsmm2-hp-tr-metric.rs:345-349](../../../mm2/src/bin/gsmm2-hp-tr-metric.rs)）；aligned-metric bin 不设 | 未设置（bandmap 默认 500，[aligner.rs:70-71](../src/aligner.rs) 只锁 best_n/pri_ratio） | bandmap 暴露 `chaining.max_gap`。若要与 mm2 bin 链合行为严格一致（大 gap 链被拆分/丢弃），补挂 `max_gap = 100`；当前 golden 不含该参数，加挂会改变行为、需重建 golden |
| `hp_tr_finder` 版本 | **0.2.3**（Cargo.lock） | **0.5.0** | 0.5.0 有 rotation-dedup 等新行为，motif 发现结果与 0.2.3 不同（迁移文档 §四已踩过版本漂移坑）。要做两版本对照验收必须先统一版本；否则应把「hp_tr_finder 版本差异」作为已知口径差异记录 |
| best_n / pri_ratio | 10000 / 0.2 | 10000 / 0.2（[aligner.rs:70-71](../src/aligner.rs)） | 一致 |
| `discard_secondary` 引擎层 flag | `OupParams::modify_aligner` 会翻转引擎 flag `0x1000000000`（[mm2 params.rs:331-335](../../../mm2/src/params.rs)）；secondary 的实际排除仍在 metric 层完成 | 不复刻引擎 flag；secondary 在 main 层无条件丢弃（`--noSeco` 仅为兼容旗标） | flag 对 map-ont 类长读 preset 大概率无效果，且两侧最终都排除 secondary、golden 验收通过；如需引擎级逐位对齐再核实该 flag 语义 |

### 2.3 CLI 语义差异

| 项 | mm2 三个 metric bin | 本项目 |
|---|---|---|
| `--short-aln` | 三个 bin 都定义了 `Cli::post_process_param` 但 **main 从不调用**——flag 全部无效（死代码） | 0.10.2 已修复接线（[main.rs:55-56](../src/main.rs)）。这是本项目领先 mm2 的修复，无需处理 |
| `--noMar` | 判定条件写错：`hits.is_empty() \|\| (hits.len() > 0 && discard_multi_align_reads)`——传入后**所有**有 hit 的 read 都被扔掉、只写占位行（[gsmm2-hp-tr-metric.rs:463](../../../mm2/src/bin/gsmm2-hp-tr-metric.rs)、[gsmm2-aligned-metric.rs:415](../../../mm2/src/bin/gsmm2-aligned-metric.rs)、[gsmm2-time-err.rs:454](../../../mm2/src/bin/gsmm2-time-err.rs)） | 按 mm2 `align_worker` 的 `records.len() > 1` 语义实现（多 chain 整条过滤，[main.rs:245-247](../src/main.rs)），并有单测。两个「noMar」语义不可互替，对齐时以本项目为准即可 |
| `--noSeco` | 不存在该 flag（`to_oup_params` 硬编码 discard_secondary=true） | 作为兼容旗标接受（0.10.3） |
| `--noSupp` | 一致（保留 primary，丢 supplementary） | 一致 |
| `--refAnchored` | 不存在 | 本项目新增 |
| 输出/日志 | tracing 打到 stdout；不写 log 文件 | 写 `{out}.log`（tracing non-blocking），hp-v2 的 `abs(cnt-expected)>=10` 调试对齐串也进 log |

### 2.4 输入过滤差异（⚠️ 发现一个疑似回归）

mm2 `query_seq_sender`（[mm2 lib.rs:245-304](../../../mm2/src/lib.rs)）与本项目重实现版（[main.rs:161-214](../src/main.rs)）逐分支对比：

| 分支 | mm2 | 本项目 | 结论 |
|---|---|---|---|
| fa | 无过滤 | 无过滤 | 一致 |
| **fq / fastq** | 计算 rq 后 **`input_filter_params.rq_valid(rq)` 过滤**（`--rq-range` 对 fastq 生效），再 send | 计算 rq 后**直接 send，无 rq 过滤** | **疑似回归**：`--rq-range` 对 fastq 输入静默失效。建议在 fq 分支补 `if input_filter_params.valid(&record) || 专门的 rq_valid` 语义（见 §六 P1）。注意 mm2 fq 分支同时有「后缀追加到 seq」的 bug，本项目修复时把 rq_valid 一并丢了 |
| bam | `valid(&record)`：np-range、rq-range（tag）、qname_suffix | `valid(&read_info)`：np-range、rq-range（tag）、qname_suffix | np/rq 一致。qname_suffix 时序不同：mm2 在**追加后缀之前**检查后缀（多 BAM 输入时后缀过滤会把全部 read 拒掉，mm2 自身的怪癖）；本项目在 valid 之后追加、等于不过滤。边缘场景，建议按本项目行为（不过滤）或干脆移除该过滤分支 |
| fq qname 后缀 | 追加到 `record.seq`（bug） | 追加到 `record.name`（已修） | 本项目正确 |

### 2.5 比对引擎差异（已知、已归档）

bandmap-core 与 minimap2 在 chimeric read 拆分链合上的引擎差异（影响 0.0046% adapter-dimer reads）已在
[mm2_bandmap_migration.md](mm2_bandmap_migration.md) §3.3 完整归因，可用 `--noMar` / `--noSupp` 收敛，此处不重复。

---

## 三、aligned metric（`gsmm2-aligned-metric` → mode `aligned`）——✅ 已迁移

> **迁移完成（2026-11-22）**：实现于本项目 `src/metrics/aligned_metric.rs`。绕开 mm2 的
> `TseqAndRecord`/rust-htslib record 机制，用 `alignment_columns` 直接按 mm2
> `fwd_record` + `alignment_pair` 的语义展开 cigar 列（reverse hit：target 侧取 RC
> 切片、cigar 反序、query 恒为原始方向）。口径决策：不做 identity 过滤（per-mode
> `Metric::IDENTITY_THRESHOLD = None`）；保留 mm2 占位行语义（unmapped / `--noMar`
> read 输出全 0 行，identity 列 `NaN`，每条输入 read 恰好一行）。验收：合成 32 条
> reads + BAM fixture 与 mm2 原版二进制输出**逐字节一致**（表头 + 排序数据行），
> `compare_golden.sh` 7/7、`compare_bam_golden.sh` 11/11 PASS。以下为迁移前的分析存档。

### 3.1 功能与输出

对每条 read 的**全部 primary + supplementary hit（raw mapping，无重对齐）**做多片段合并统计，输出
`qname rname qlen segs covlen primaryCovlen queryCoverage identity oriAlignInfo oriQGaps qOvlp qOvlpRatio rOvlpRatio mergedQrySpan match misMatch ins homoIns del homoDel longIndel` + 24 列 per-base（`match-A/C/G/T`、`misMatch-*`、`ins-*`、`homoIns-*`、`del-*`、`homoDel-*`）。

核心逻辑（[mm2 align_metric.rs](../../../mm2/src/align_processor/align_metric.rs)）：

1. `AlignInfo::fwd_record` 把每个 hit（含 reverse）转成 forward 的「target 子串 + 截取后 record」
   （reverse hit：target 取 RC、cigar 反序、query 切原始方向子串，[align_processor/mod.rs:150-198](../../../mm2/src/align_processor/mod.rs)）。
2. 多 hit 按 query 区间做**去重叠截断**（后 hit 与前一 hit 重叠时裁短先到者，[align_metric.rs:140-167](../../../mm2/src/align_processor/align_metric.rs)），据此得 `num_segs` / `covlen` / `qOvlp` / `mergedQrySpan` / `oriQGaps`（Regions 运算）。
3. 对每个截断后片段重放 `alignment_pair` 得逐列比对串，`fill_metric_from_align_pair` 扫描：
   - match / mismatch 按 **query 碱基** 计入 per-base 桶；
   - **homo ins**：连续插入碱基一致，且与左/右 ref 邻碱基一致（[align_metric.rs:415-449](../../../mm2/src/align_processor/align_metric.rs)）；
   - **homo del（v2 判定）**：被删 ref 碱基与**自身**左右 ref 邻碱基任一一致即 homo（逐碱基判定，[align_metric.rs:385-410](../../../mm2/src/align_processor/align_metric.rs)）；
   - **longIndel**：连续 indel ≥ 10 记录到 `longIndel` 列（ins 正值 / del 负值）。
4. `identity` 列 = `matched / aligned_span`（不是 cigar identity）；`queryCoverage` = `covlen/qlen`。

### 3.2 迁移要点（预估 1~2 天）

- **纯算法部分可近乎照搬**：`fill_metric_from_align_pair`、`fill_deletion_info_v2`、`fill_insertion_info`、
  `get_pre/next_rbase`、多片段截断逻辑都是 Vec/字符串纯函数，不依赖 mm2 API。
- **record 机制可绕开或复用**：mm2 用 rust-htslib `BamRecord` + `aligned_pairs_full` 做列重放；本项目
  `kits/bam_ext.rs` 已有等价的 `aligned_pairs_full`（bam2-msa 在用）。更简单的做法是直接用
  `Mapping::aligned_2_str` / `aligned_pairs()`（adapter 已有），不必重建 record。
- **过滤策略需决策**：mm2 原版无 identity 过滤（它本身就是比对质量评估工具）；迁移后是否套用 main 层的
  `identity > 0.75` / `--noMar`？建议 aligned metric 走「不过 identity 滤」的独立路径，否则会把低质量 read
  （恰恰是它要度量的对象）静默丢掉。
- **占位行语义需决策**：mm2 对 unmapped read 输出全 0 行（identity 列为 `NaN`）；本项目惯例是不产出行。
- 该 metric 消费 raw mapping，reverse hit 的坐标口径直接复用 `adapter::Mapping` 即可（bandmap 坐标已验证与
  mm2 一致，见迁移文档 §二）。

---

## 四、未迁移功能二：time-err metric（`gsmm2-time-err`）

### 4.1 功能与输出

按**测序时间**把一条 read 的比对算子分桶（默认 600s/桶），输出 `qname rname minutes eq diff ins del` 多行。
时间轴来自 BAM 的 `ar`（advance）与 `dw`（dwell）tag：`base_cycles[i] = cumsum(ar[i] + dw[i])`，
`base_secs = cycles × 0.002`（[gsmm2-time-err.rs:281-299, 480-484](../../../mm2/src/bin/gsmm2-time-err.rs)）。
每个 cigar op 落到其覆盖的 query 碱基的时间桶：op∈{1(I),7(=),8(X)} 逐碱基推进 qcursor，op 2(D) 记在当前
qcursor 上（[time_err_metric.rs:111-153](../../../mm2/src/align_processor/time_err_metric.rs)）。

### 4.2 迁移要点（预估 1 天 + 前置依赖）

1. **恢复 `ar`/`dw` tag 提取设备**。0.10.3 删除了 `pass_through_tags` 与 `BamRecordExt` 的整套 tag 提取
   （`ReadInfo` 只保留 name/seq/np/rq/qual）。迁移需在 `kits/bam_ext.rs` / `kits/read_info.rs` 按需恢复
   `ar`/`dw` 的读取（不需要完整 pass_through_tags 设备，做两个具名 getter 即可）。
2. **reverse hit 的 qcursor 口径疑似 mm2 既有 bug**：mm2 的 `TAlignedPairs` 对 reverse hit 用
   `query_cursor = qlen - query_end` 起步（[mm2 aligned_pairs.rs:59-63](../../../mm2/src/aligned_pairs.rs)），
   而 `TimeErrMetric` 对 reverse hit 直接 `qcursor = align_info.qstart` 且 cigar 反序迭代
   （[time_err_metric.rs:115-122](../../../mm2/src/align_processor/time_err_metric.rs)）。两者不可能同时正确
   （bandmap 与 mm2 的 query 坐标都是原始正向口径，迁移文档 §二已确认）。**迁移时应按 `qlen - query_end`
   起步实现，不要照抄 mm2 的 qstart 写法**；如需与 mm2 旧输出对账，先确认旧输出的 reverse read 分桶是否有偏差。
3. cigar op 编码（1/2/7/8）与 EQX 语义 bandmap 一致，cigar 不含 clip，逐 op 分桶逻辑可照搬。
4. mm2 bin 强制 BAM 输入（需要 ar/dw tag）——迁移后保持 BAM-only assert。
5. 空输出兜底行（`qname\t\t0\t0\t0\t0`）是否保留需与下游确认（本项目惯例：不产出行）。
6. 分桶参数 600s 在 mm2 中是硬编码（CLI 不可调）；迁移时可考虑提为 `--interval`，但为对齐口径默认 600。

---

## 五、明确不迁移（当前判断）

| mm2 内容 | 理由 |
|---|---|
| `gsmm2.rs` align 主 bin 及其 `--hpTrShrinkAlnRegion` / `--polyNGapLeftAlign` / `--query-forward`（`hp_tr_shrink.rs`、`cigar_adjust.rs`、`bam_writer.rs`、`mapping_ext.rs` 的对齐路径部分） | 属比对/BAM 输出管线，不是 metric。若未来想在本项目内做「比对+评估」闭环，可作为 phase 2 再评估 |
| `bam-rw.rs` | BAM 拷贝工具，与 metric 无关 |
| `visualization.rs` | 调试辅助；如需要，可基于 `Mapping::aligned_2_str` 十几行复刻，优先级最低 |

---

## 六、建议处理项（按优先级）

| 优先级 | 事项 | 类型 |
|---|---|---|
| **P1** | 修复 fastq 分支 `--rq-range` 失效：`query_seq_sender` 的 fq 分支补 rq 过滤（对齐 mm2 `rq_valid` 行为），补一条回归测试 | 回归修复 |
| **P1** | ~~决策并迁移 **aligned metric**~~ **✅ 已完成（2026-11-22）**：mode `aligned`，与 mm2 原版二进制逐字节一致，见 §三 | 功能迁移 |
| **P2** | 迁移 **time-err**：先恢复 `ar`/`dw` 具名 tag 读取；reverse hit qcursor 按 `qlen - qend` 实现（不照抄 mm2 疑似 bug 的 qstart 写法）；保留 600s 硬编码默认 | 功能迁移 |
| **P2** | 与下游确认 hp-tr-v2 的三个口径差异是否需要提供兼容开关：① per-strand 双行（含全 0 行）；② unmapped/无 motif read 的占位行/空行；③ 默认输出文件名 | 口径决策 |
| **P3** | 如需与 mm2 bin 链合行为严格对齐，评估补挂 `chaining.max_gap = 100`（bandmap 已暴露该参数；改动会变更行为，需重建 golden） | 参数对齐 |
| **P3** | 把 `hp_tr_finder 0.2.3 vs 0.5.0` 的差异作为已知口径记录（或做一次两版本 motif 发现 diff 报告），避免未来对照验收时再次踩版本漂移坑 | 文档/验收 |
| **P3** | `--qname-suffix` + 多 BAM 输入的边缘行为：mm2 会拒掉全部 read（后缀追加前检查）；本项目不过滤。保持现状并记录，或干脆从 valid() 中移除 qname_suffix 检查 | 边缘口径 |

---

## 七、证据索引（关键 file:line）

**mm2 侧**
- `src/bin/gsmm2-hp-tr-metric.rs:345-349`（max_gap=100）、`:351-360`（BAM-only assert）、`:463`（noMar 判定 bug）、`:452-509`（build_hp_tr_metric：raw mapping 直进 metric）
- `src/align_processor/hp_tr_metric.rs:69-93`（[fwd,rev] 双行 + RC motif 名）、`:95-131`（compute_metric）
- `src/bin/gsmm2-aligned-metric.rs:415`（noMar 判定 bug）、`src/align_processor/align_metric.rs`（metric 全部算法）
- `src/bin/gsmm2-time-err.rs:281-299`（ar/dw→base_secs）、`:454`（noMar 判定 bug）、`src/align_processor/time_err_metric.rs:111-153`（分桶，reverse qcursor 口径）
- `src/lib.rs:245-304`（query_seq_sender：fq rq_valid / seq 后缀 bug）、`:372-402`（metric 路径不做事后处理）
- `src/aligned_pairs.rs:59-63`（reverse 的 qlen-qend 口径）、`src/params.rs:63-70`（rq_valid）、`:331-335`（discard_secondary 引擎 flag）
- `Cargo.toml` / `Cargo.lock`：minimap2 =0.1.23、hp_tr_finder 0.2.3

**本项目侧**
- `src/main.rs:44-53`（常量）、`:161-214`（query_seq_sender，fq 分支缺 rq 过滤）、`:245-247`（has_multiple_chains）、`:249-302`（compute_metric 过滤链）
- `src/aligner.rs:70-71`（best_n/pri_ratio，无 max_gap）
- `src/metrics/realign.rs`（共用重对齐前缀）、`src/metrics/hp_tr_tools.rs`（shrink/子串/do_align_4_homo）
- `src/metrics/hp_tr_metric_v2.rs`、`hp_tr_metric_v3.rs`、`hp_metric_v2.rs`、`src/msa.rs`（bam2-msa）
- `src/params.rs:100-143`（bandmap 参数映射）、`src/kits/read_info.rs`（ReadInfo 仅 name/seq/np/rq/qual）、`src/kits/bam_ext.rs`（aligned_pairs_full，可复用于 aligned metric）
- `docs/mm2_bandmap_migration.md`（引擎替换验收）、`docs/CHANGELOG.md`（0.8.x-0.10.3 演化记录）
- `Cargo.toml` / `Cargo.lock`：bandmap-core 2.0.0、hp_tr_finder 0.5.0
