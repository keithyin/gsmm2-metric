# gsmm2-metric 代码评审与重构方案

> 评审范围：`src/` 全部 3668 行（main.rs、cli.rs、aligner.rs、global_data.rs、aligned_pairs.rs、msa.rs、metrics/*）。
> 结论先行：**架构方向是对的（trait 泛型流水线），但三个活跃 metric 之间存在约 80 行 × 3 的复制粘贴前缀逻辑，trait 过胖、输出无 schema、死代码未清理**。重构应以「抽出共用重对齐前缀」为核心，分 5 个阶段增量进行，全程用 golden test 锁行为。

---

## 一、现状评价

### 1.1 做得好的地方

| 方面 | 证据 |
|---|---|
| 流水线抽象 | `metric_entrance::<M>` 用泛型把「读入 → 比对 → 算 metric → 落盘」三段线程模型统一，新增 mode 只需实现 `TMetric`（main.rs:76-131） |
| 并发模型清晰 | crossbeam 有界 channel + `thread::scope`，读入/计算/写出职责分离 |
| 纯函数抽取与单测 | `count_repeats`、`build_ref2qpos`、`shrink_to_motif_boundary`、`get_target_substr` 都是可独立测试的纯函数，且带回归测试（含边界案例，如 shrink 收缩成空返回 `None` 而非 panic） |
| 惰性全局数据 | `GlobalData` 用 `OnceCell` 按需初始化反转序列 / motif 索引，避免无谓的预计算 |
| 文档 | README 对 mode 差异、列语义、边界行为的记录非常完整 |

### 1.2 主要问题

#### P1 复用性：重对齐前缀逻辑三处复制（最核心）

`hp_tr_metric_v2.rs:114-188`、`hp_tr_metric_v3.rs:179-254`、`hp_metric_v2.rs:107-158` 三段 `compute_metric` 的开头几乎逐行相同：

```
取 target 序列（GlobalData match 解包） → [可选] shrink_to_motif_boundary
→ get_target_substr → 截取 read 子串 → [refAnchored] 双向反向互补
→ do_align_4_homo → None/仍 reverse 则 warn + continue → single_seq_hp_tr_finder 找 motif
```

差异只有两点：v2/v3 做 boundary shrink 而 hp-v2 不做；正则表不同（`HP_TR_REG` unit 1-4 vs `HP_REG` unit 1）。这 80 行 × 3 的复制意味着：
- 任何一处 bug 修复（如 `raw_target_end == 0` 的防御、warn 日志字段）都要改 3 处；
- v3 修了 `get_target_substr` 的越界防御，v2/hp-v2 靠的是同一个函数才间接获益；
- 新增第 4 个 mode 时大概率继续复制第 4 份。

伴生的重复：
- **Counter** 四份：`hp_tr_metric_v3::Counter`、`hp_tr_metric_v2::AlignCounter`、`hp_metric_v2::{Counter, CounterWithAlignSpan}`，外加死代码 `hp_metric.rs` 里一份；
- **正则表** 两份且语义重复：`HP_TR_REG`（v2 文件）与 `HP_REG`（hp_v2 文件），后者只是前者去掉 unit 2-4；
- **build_metric_str** 两份近乎逐行相同（v3 与 hp-v2，只差是否输出 `ref` 列）；
- **测试样板** 四份 × 100+ 行：build_aligner 参数构造 + 两条几百碱基的硬编码序列字面量在 3 个文件的 `#[cfg(test)]` 里重复出现。

#### P2 扩展性

1. **`TMetric` trait 过胖且职责混杂**（mod.rs:12-43）。一个 metric 实现要被迫写 8 个 getter/setter 样板（qname/tname/global_data/mappings/metric_str），trait 同时管「身份」「数据源」「输出」三件事。`metric_str: Option<String>` + `set_metric_str()` 默认方法 + `get_metric_str().as_ref().unwrap()` 的链条，把「计算」和「序列化」用可变状态耦死。
2. **配置传递路径不统一**。`reference_anchored` 走 `compute_metric` 的参数；`discard_supplementary`/identity 阈值在 main.rs 里硬编码过滤；重对齐参数（`a=4 b=10 q=4...`）写死在 `do_align_4_homo`。未来想按 mode 调整任何一个，都要动不同层。
3. **motif 索引两套体系并存**。HP 用预计算的 `TargetRegion2Motif4Hp`（IntervalTree），HP-TR 用 `TargetRegion2Motif4HpTr`（flatten 后的 HashMap），而 read 级别的 motif 查找又走 `single_seq_hp_tr_finder` 现算。三种获取 motif 的方式没有统一抽象。
4. **`GlobalDataKey/Value` enum 手工分派**（global_data.rs）。加一种全局数据要改 enum + `new()` + `get()` 三处 match，是典型的「用 enum 模拟字段」。类型安全的 `OnceCell<T>` 直接做 struct 字段即可。

#### P3 维护性

1. **死代码**：`hp_metric.rs`（436 行）、`hp_tr_metric.rs`（303 行）未接入 main 分发，执行即报错；且报错文案漏了 `hp-tr-v3`（main.rs:70）。死代码还拖着 `BASE_COMP`（mod.rs）等只为它存在的符号。
2. **输出无 schema**：所有 metric 输出是手工 `join("\t")` 的 String，列名在 `csv_header()` 里手写一遍、取值在 `build_metric_str()` 里再手写一遍，两处靠人肉对齐。同名列跨 mode 语义不同（`eq` 在 v2/v3/hp-v2 是三码事），只能靠 README 约束，代码层面无类型保护。文件扩展名 `.csv` 实际是 TSV 也是历史包袱。
3. **错误处理粗糙**：`panic!("")`（空消息）、`expect(&format!(...))` 反模式（main.rs:35、207）、多层 `unwrap()`、`get_target_substr` 用「返回空串」当错误信号（调用方还要再判一次空）。
4. **魔法数**：`threads - 4`、`assert!(tot_threads >= 10)`、`identity > 0.75`、`best_n = 10000`、`pri_ratio = 0.2`、`abs(cnt - expected) >= 10`（hp-v2 的调试打印阈值）散落各处无常量无注释。
5. **类型混用**：usize/i64 反复 cast（aligned_pairs、各 metric 的坐标计算），是 off-by-one 风险聚集地；`pre_ins_cnt` 是 i32 而计数是 usize。
6. **测试有效性不均**：v3 的 `test_str_called_counting` / `test_count_repeats` 有真实断言，质量好；但 v2 / hp-v2 的测试只 `println!` 不断言，且依赖 `test_data` 相对路径（`cargo test` 必须在 crate 根目录跑）。
7. **聚合语义陷阱**：metric_core 的 key 只有 motif 字符串，同 read 内同名 motif 不同位点被静默合并（README 已注明）。这是设计选择，但代码里没有任何注释/类型提醒新维护者。

---

## 二、重构方案

### 总原则

1. **行为锁定优先**：先建 golden test（固定输入 → 固定 TSV 输出字节级对比），每个阶段跑通后再进入下一阶段。输出格式的任何变化（如加位点列）必须是显式声明的独立变更。
2. **先抽共性、再收接口、最后清理**：复制粘贴是当前最大的维护成本，优先消灭。
3. **不引入新依赖**（最多加 `serde`/`thiserror`，可选）。

### Phase 0：安全网（半天）

- 在 `test_data/` 中固化 2-3 个小输入（fastq + fasta，覆盖 fwd/reverse/refAnchored 三种路径），对 `hp-tr-v2`、`hp-tr-v3`、`hp-v2` 各生成一份期望输出，写成 `tests/golden.rs`。
- 把现有 4 份测试样板抽成 `#[cfg(test)] pub mod test_support`：`fn build_test_aligner(preset) -> Vec<NoMemLeakAligner>` + fixture 序列改从 `test_data/*.fasta` 读。测试断言缺失的（v2、hp-v2）顺手补上基本断言。

### Phase 1：抽出共用「重对齐前缀」（核心，1-2 天）

新建 `src/metrics/realign.rs`，把三处复制的前缀逻辑收敛为一个纯流程函数：

```rust
/// mode 无关的重对齐前缀流水线。
pub struct RealignConfig {
    pub shrink_to_motif_boundary: bool, // v2/v3 = true, hp-v2 = false
    pub ref_anchored: bool,
    pub motif_regs: Arc<Vec<HashMap<String, Regex>>>, // 由 mode 提供自己的正则表
}

pub struct RealignedHit {
    pub read_seq: String,       // 已按 refAnchored 归一
    pub target_substr: String,  // 已做 homopolymer 收缩
    pub align: Mapping,         // do_align_4_homo 的结果，保证 Forward
    pub region2motif: Vec<((usize, usize), Arc<String>)>, // 在 target_substr 上找到的 motif
}

/// 返回 None 表示该 hit 应跳过（收缩成空 / 对齐失败 / 仍 reverse），并负责打对应 warn。
pub fn realign_hit(
    read_info: &ReadInfo,
    mapping: &Mapping,
    global_data: &GlobalData,
    cfg: &RealignConfig,
) -> Option<RealignedHit>;
```

三个 metric 的 `compute_metric` 收缩为：`for hit in mappings { let Some(h) = realign_hit(...) else { continue }; <mode 专属统计>(h) }`。

配套统一：
- **MotifSpec** 取代重复的 lazy_static：每个 mode 声明 `{ unit_len_min, unit_len_max, min_repeats }`，运行期构建正则（或用 `once_cell` 缓存 per-spec 的表）。`HP_TR_REG`/`HP_REG` 归一为同一构造函数的两个实例。
- warn 日志统一在一处打，字段（qname、qstart/qend、tstart/tend、strand）一致。

**收益**：消灭 ~160 行重复；`do_align_4_homo` 的参数成为 `RealignConfig` 一部分（仍给默认值），为将来按 mode 调参留口。

### Phase 2：收紧 `TMetric` trait（1 天）

把 trait 从「胖接口 + setter 生命周期」改为「上下文注入 + 有类型输出」：

```rust
/// 一次性携带计算所需的全部上下文，取代 8 个 setter/getter。
pub struct MetricContext<'a> {
    pub read: &'a ReadInfo,
    pub target_name: &'a str,
    pub mappings: Vec<Mapping>,        // 已经过 main.rs 的 primary/supp/identity 过滤
    pub global_data: &'a GlobalData,
    pub config: &'a PipelineConfig,    // reference_anchored 等，取代散落的参数
}

pub trait Metric: Send + Default {
    fn header() -> &'static [&'static str];
    fn compute(&mut self, ctx: &MetricContext);
    fn rows(&self) -> Vec<Row>;        // 见下
}

/// 类型化输出行，取代手工 join("\t")。
pub enum Cell { S(String), U(usize) }
pub struct Row(Vec<Cell>);
```

- main.rs 的 `dump_metric_worker` 改为：写 header（由 `Metric::header()` 生成）→ 逐行 `row.join('\t')`。header 与取值不再两处手写（列顺序由 `rows()` 的构造顺序单点决定）。
- 各 metric struct 只保留自己的聚合状态，删掉 `qname/tname/global_data/align_infos/metric_str` 五个字段和全部 unwrap。
- 阶段内保持输出字节不变（header 仍是 TSV，`.csv` 扩展名暂保留并在 header 注释/文档标注「实为 TSV」）。

### Phase 3：统一计数与输出语义（1 天，含一处显式行为增强）

- 合并 `Counter` 家族为一个 `RepeatCounter`（`counts: [HashMap<usize, usize>; 2]` + eq + align_span），v2 的 `AlignCounter` 独立保留（它统计的是算子谱，语义确实不同，但放同一个模块并命名清楚）。
- 输出列改造（**显式行为变更，单独一次 commit + README 更新**）：
  - `hp-tr-v3` / `hp-v2` 增加 `tstart`/`tend` 两列（motif 在 target 上的坐标），解决「同名 motif 不同位点被合并」的下游无法区分问题；不改列名与既有列顺序，新增列追加在末尾，下游按列名解析可兼容。
  - 或保守方案：不改格式，仅在 `RepeatCounter` 与 `Row` 上加 doc comment + debug_assert 提醒聚合语义。
- `parse_motif` 移入公共模块（`hp_tr_tools` 或 `motif.rs`），hp-v2 删掉 `motif.as_bytes()[1]` / `motif[3..]` 的位置硬编码，统一走 `parse_motif`，顺带修掉它对 `(X)N` 以外形式的脆弱性。

### Phase 4：GlobalData 与基础设施（半天）

- `GlobalData` 改为普通 struct 字段，删 enum 分派：

```rust
pub struct GlobalData {
    pub target_seqs: HashMap<Arc<String>, Arc<String>>,
    seqs_with_rev: OnceCell<HashMap<Arc<String>, [Arc<String>; 2]>>,
    region2motif_hptr: OnceCell<HashMap<Arc<String>, Arc<RefPos2Motifs>>>,
    region2motif_hp: OnceCell<HashMap<Arc<String>, Arc<HpIntervalTree>>>,
}
impl GlobalData {
    pub fn seqs_with_rev(&self) -> &...;        // get_or_init
    pub fn region2motif_hptr(&self) -> &...;
}
```

  类型即文档，加新数据源不再需要动三处 match。`TargetRegion2Motif4Hp` 与 `4HpTr` 的索引结构差异（IntervalTree vs flatten HashMap）在 Phase 1 之后可以评估统一为 IntervalTree 一种（v2/v3 目前只按端点查一次，flatten HashMap 是够用的；统一与否按性能实测决定，不强求）。
- 错误处理：`thiserror` 定义 `MetricError`，`realign_hit` 返回 `Result<Option<RealignedHit>, MetricError>` 或至少替换 `panic!("")` / 空 expect；warn 保留为「跳过该 hit」的正常路径。
- 魔法数常量化：`MIN_THREADS = 10`、`OUTPUT_THREADS = 4`、`MIN_IDENTITY = 0.75`、`BEST_N`、`PRI_RATIO`，集中到 `PipelineConfig`。

### Phase 5：清理（半天）

- 删除 `hp_metric.rs`、`hp_tr_metric.rs` 及 `BASE_COMP`、`Mode::{HpTr, Hp}`（或保留 Mode 变体但让报错文案列出全部合法 mode）。git 历史可找回，不值得留 739 行死代码。
- main.rs 的报错改为 `Mode::iter()` 动态生成合法列表。
- 修 v2 文件里遗留的 `use hp_tr_metric::HpTrMetric` 测试导入（测试实际在测旧 metric，Phase 5 后改为测 `HpTrMetricV2`）。
- （可选）`aligned_pairs` 的 `Box<dyn Iterator>` 换成 `impl Iterator` 返回已可达成（当前写法已经用了 `impl Trait` 签名，内部 Box 是实现细节，可留）；补 op=0（`Match`，eqx required）的 panic 消息说明依赖 `with_cigar()`。

### 明确不做的事

- **不动 msa.rs / bam2-msa**：口径完全不同，与 metric 流水线无共享逻辑，强行统一收益低。
- **不改三段线程模型**（读入/计算/落盘、threads-4）：结构正确，只做常量化。
- **不改默认输出文件名**（`.csv` 扩展名）：下游兼容性优先，文档说明即可。
- **不动重对齐打分参数默认值**：属实验调参，不属于重构。

---

## 三、阶段收益预估与顺序依赖

| Phase | 删减/收敛 | 依赖 |
|---|---|---|
| 0 安全网 | +测试样板去重（~300 行测试代码收敛为 1 份） | 无 |
| 1 重对齐前缀 | metric 文件合计 -160~200 行；bug 修复面从 3 处变 1 处 | Phase 0 |
| 2 trait 收紧 | 每 metric -40 行样板；unwrap 链清零 | Phase 1（compute 里直接调 realign_hit） |
| 3 计数/输出统一 | Counter -60 行；列语义类型化 | Phase 2 |
| 4 GlobalData/错误处理 | global_data.rs 约 -50 行；panic 清零 | Phase 1（realigned 已消费 GlobalData） |
| 5 清理 | -700+ 行死代码 | Phase 0（golden test 保证误删可发现） |

总计约 1000 行的净减/收敛（当前 3668 行，预计降到 ~2700），核心收益不在行数，而在于：

1. **新增一个 metric mode 的成本**从「复制 80 行前缀 + 8 个 setter + 手写两遍列名」降到「实现 `Metric::compute` + `rows` + 声明 `MotifSpec`」；
2. 共用逻辑（边界收缩、refAnchored、重对齐、motif 发现）有了唯一实现，语义差异（shrink 开关、正则表、列集合）全部变成声明式配置；
3. 输出 schema 有了类型承载，跨 mode 同名列异义的坑从 README 层下沉到代码层。

## 四、风险与回退

- **最大风险**是 Phase 1/2 改动触碰三个 metric 的统计口径。缓解：Phase 0 的 golden test 覆盖 fwd / reverse / refAnchored / supplementary 四种输入路径后再动手；每个 phase 独立 commit，可二分回退。
- **Phase 3 加位点列**是唯一的对外可见变更，需与下游消费方确认（或先走保守方案不改格式）。
- `hp-v2` 的「diff>=10 打对齐串到 log」调试行为在重构中保留（挪进 mode 专属统计函数），避免丢失排障能力。
