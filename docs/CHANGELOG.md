# Changelog

格式参照 Keep a Changelog，版本与 git tag（`v*`）一一对应；0.8.1 之前未逐条记录。

## Unreleased

### Added

- `--mode aligned`：全比对质量 metric（mm2 仓库 `gsmm2-aligned-metric` bin 的移植，mm2 三 metric bin 中最后一个迁入）。直接消费 raw mapping（不重对齐、不做 identity 过滤），多片段按 query 区间去重叠截断后合并统计，输出 45 列（segs/covlen/identity/homo ins-del/longIndel/per-base ACGT 分解）。保留 mm2 的「每条输入 read 恰好一行」语义：unmapped / `--noMar` 过滤的 read 输出全 0 占位行（identity 列为 `NaN`）；reverse hit 的 target 侧取反向互补切片、cigar 反序后按 read 方向展开（与 mm2 `TseqAndRecord::alignment_pair` 逐列一致，per-base 桶按 read 原始方向碱基计数）。golden 以 mm2 原版二进制输出为基线（fasta + BAM fixture），`--refAnchored` 对该 mode 无效。
- `src/kits/regions.rs`：gskits `ds::region`（`Regions::ovlp_ratio` / `gaps`）的本地移植；`src/kits/dna.rs` 移植 `SEQ_NT4_TABLE`。

### Changed

- `Metric::compute` 返回类型 `MetricOutcome` → `Vec<String>`，删除 `MetricOutcome`（`Rows(vec![])` 与 `Skipped` 行为等价）；`--noMar` 命中的 read 改走与 unmapped 相同的空-mappings 路径（hp 系 mode 行为不变，`aligned` 借此产出占位行）。
- `Metric` trait 新增关联常量 `IDENTITY_THRESHOLD: Option<f32>`（默认 `Some(0.75)`，即原 `main.rs` 硬编码行为；`aligned` 覆盖为 `None`），identity 过滤由 main 层按 mode 门控。

## 0.10.3 — 2026-11-22

行为锁定的整体重构：全程 golden 字节级验收（fasta 6/6 + BAM 10/10），对外 CLI、输出格式、列集合均不变。

### Added

- CI（ci.yml）新增门禁：`cargo fmt --check`、`cargo clippy -D warnings`、release 构建后运行 `scripts/compare_golden.sh` / `compare_bam_golden.sh` 字节级 golden 验收。
- 测试补强：hp-v2 / hp-tr-v2 端到端测试从纯 println 改为排序后逐行断言（基线以 Cargo.lock 锁定的 bandmap-core 捕获），`test_data` 路径改用 `CARGO_MANIFEST_DIR` 消除 cwd 依赖；新增 `get_target_substr` 序列末尾窗口回归测试。
- `--noSeco` 兼容旗标：仅被 CLI 接受并在 log 留痕，不进 `OupParams`。本管线中重比对产生的 secondary hit 在 `compute_metric` 保留过滤中本就被无条件丢弃（恒等价于 `--noSeco` 开启），旗标无额外行为，用于 mm2 版命令行兼容；顺带修正 `compute_metric` 注释中对已删除的 `discard_secondary` 过滤的指称。

### Changed

- 抽取三个 metric 复制的 ~80 行重对齐前缀为 `metrics/realign.rs::realign_hit`（mode 差异收敛为 `RealignConfig`），`HP_TR_REG` / `HP_REG` 正则表由 `MotifSpec` 统一派生。
- `TMetric` trait 收紧为 `Metric`：列名 `COLUMNS` 单点声明，`MetricContext` 注入取代 9 个 accessor，「该 read 无输出」由 `MetricOutcome::Skipped` 显式表达（取代空 `metric_str` 哨兵），metric struct 只保留聚合状态。
- `GlobalData` 的 enum 手工分派改为类型化 struct + `OnceCell` 字段；删除无消费者的 `TargetName2Seq` / `TargetRegion2Motif4Hp` 索引。
- 删除死代码与死配置：`MapParams`/`MapArgs`（CLI 面不变）、`OupParams::discard_secondary`（从未生效，secondary 在 main 层无条件丢弃）、`oup_identity_threshold` / `oup_coverage_threshold`（恒 -1.0 且从未被查询）、`pass_through_tags` 及 `ReadInfo`/`BamRecordExt` 的整套 tag 提取设备（`ReadInfo` 仅保留实际消费的 name/seq/np/rq/qual）、无引用依赖 `ndarray`、空跑测试 `test_hp_tr_metric_case`（其 query 在 bandmap 下比对不到任何 hit）。
- 魔法数常量化（线程模型 / identity 阈值 / 排序加分 / 重对齐打分），命名修正（`msa_entrence`→`msa_entrance` 等），`kits/dna.rs` 的 unsafe transmute 改安全实现。

### Fixed

- `get_target_substr` 窗口 `[len-1, len)` 顶到序列末尾时起点收缩循环越界 panic（Bug 8）。
- 默认输出文件名 `rsplit_once(".")` 在 query 路径无扩展名时 panic、目录名含「.」时切错位置，改用 `Path::with_extension`（Bug 9）。
- bam2-msa Shannon 熵静默丢弃 N 等非 -/A/C/G/T 碱基（Bug 10）。
- bam2-msa `compute_max_ins_of_each_ref_position` 的 end 边界误用 `.max()`，调用方传入的 `rend` 上限被静默放大回 record 末尾（Bug 6，当前调用点未触发）；`draw_msa` 对 pos2ins 查不到的 ref 位点跳过该柱不再 panic（Bug 6）。
- metric 输出与 msa 输出显式 flush 并检查错误（drop-flush 会吞掉磁盘满等错误）；primary target 提取不再对无 primary 的 hit 列表越界 panic（Bug 7）。

## 0.10.2 — 2026-11-22

### Added

- `--noMar` 实装：query 比对产生不止一个 chain（存在 supplementary 或 secondary mapping）时，整条 query 过滤掉、不产出输出行（复用「空 metric_str 不落盘」机制，静默跳过）。判定先于 same-target / identity 过滤（跨 target 的 alt chain 也参与判定），不做 identity 阈值过滤；语义对应迁移前 mm2 `align_worker` 的 `records.len() > 1 -> continue`。新增 `has_multiple_chains` 纯函数与 5 用例单测，split-read 端到端验证通过。

### Fixed

- `--short-aln` 实装：`post_process_param` 此前是死代码（定义于 `cli.rs` 但无任何调用点），传入后六个短比对参数（`kmer=5`/`wins=3`/`min_cnt=3`/`min_dp_max=20`/`min_chain_score=5`/`min_ksw_len=0`）全部不生效，短序列因默认 k=15/w=10 大量漏比对。修复为 `main()` 在 `Cli::parse()` 后调用 `cli.post_process_param()`；仅填充用户未显式指定的 `Option` 参数，显式值不被覆盖（[bug-report.md](../bug-report.md) Bug 2）。
- `--noMar` 传入后程序启动即 panic：`OupParams` setter 的互斥断言与 `to_oup_params` 硬编码的 `set_discard_secondary(true)` 冲突（[bug-report.md](../bug-report.md) Bug 1）。

### Changed

- 移除 `OupParams::set_discard_secondary` / `set_discard_supplementary` / `set_discard_multi_align_reads` 的互斥断言，三个 flag 相互独立；`--noMar` 与 `--noSupp` 组合合法但无实际效果（存在 supp 时 `--noMar` 已过滤整条 read）。
- README mode 表移除已删除的 `hp-tr` / `hp` 两行（v0.8.2 起 v1 metric 文件与 mode 已不存在，clap 在参数解析阶段即拒绝）；`--noMar` / `--noSupp` 描述与共用流水线清单同步实际语义。

## 0.10.1 — 2026-11-21

### Changed

- CI 迁移：GitLab CI → GitHub Actions，拆分为 `ci.yml`（build + test）与 `release.yml`（v* tag 触发 crates.io 发布与二进制 release）；README 增加 CI/crates.io/license 徽章。

## 0.10.0 — 2026-11-21

### Changed

- refactor!: rust-htslib / gskits 依赖替换为 noodles（BAM/FASTX 读取），gskits 子集就地 vendor 进 `src/kits`；移植时顺带修复 gskits `sp` tag 实读 `sd` 的笔误。
- bam2-msa 的 BGZF 解压线程数等价于原 `set_threads` 语义。

## 0.9.0 — 2026-11-20

### Changed

- feat!: 比对引擎由 mm2 0.25.4（minimap2 2.28 FFI 封装）替换为 bandmap-core 2.0.0（纯 Rust minimap2 算法复刻）。一致性验收：合成数据 6/6 逐字节一致；真实数据 54 万条 reads 行级差异率 ≤0.0025%（归因于 chimeric read 拆分链合的引擎差异），详见 [mm2_bandmap_migration.md](mm2_bandmap_migration.md)。

## 0.8.3 — 2026-11-20

- 版本维护，无功能变更。

## 0.8.2 — 2026-11-20

### Fixed

- 移除未接入 main 分发的 v1 `hp_metric.rs` / `hp_tr_metric.rs` 及对应 `Mode` 变体（无效 mode 现由 clap 在参数解析阶段拒绝）。
- 重实现 `query_seq_sender`：修复多 fastq 输入时 qname 后缀被追加到序列而非 read 名的 mm2 上游 bug（fa/bam 分支行为与 mm2 保持一致）。
- hp-tr-v3：motif 起点前一个 ref 位点发生 deletion 时 query 窗口起始错位、called 系统性少算 1 个 unit（`build_ref2qpos` 增补 `is_deleted` 语义，q_start 按 deletion 分支取值），补边界 deletion 回归测试。
- hp-v2：取窗口对 insertion 前算后不算的不对称（trailing insertion 被丢弃），逻辑抽成 `count_window` 纯函数并补 5 个回归测试。

## 0.8.1 — 2026-11-19

### Added

- feat: HpTrMetricV3（`(ref repeats, called repeats)` 重复数分布统计，mode `hp-tr-v3`）。

### Fixed

- 比对端点落在 motif 区域内部导致窗口塌缩时跳过该 hit（`shrink_to_motif_boundary` 返回 `None` 而非 panic）。
