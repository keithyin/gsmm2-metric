# gsmm2-metric Code Review — Bug Report

Review 范围：`src/` 全部 Rust 文件（main.rs、cli.rs、aligner.rs、global_data.rs、aligned_pairs.rs、msa.rs、metrics/*），并交叉核对了依赖 `mm2 0.25.6`、`hp_tr_finder 0.4.0`、`minimap2 0.1.23+minimap2.2.28`、`rust-htslib 0.49.0` 的源码（Cargo.lock 锁定版本）。注：报告撰写时为 13 个文件；2026-11-20 移除未接线的 v1 `hp_metric.rs` / `hp_tr_metric.rs` 后为 11 个，涉及条目已就地更新。

验证情况：`cargo check --all-targets` 通过；其中 Bug #1 已用 release 二进制 + 最小输入实测复现。

---

## 高严重度

### Bug 1：`--noMar` 传入后程序启动即 panic（已实测复现）

- 位置：`src/cli.rs:283-292`（`OupArgs::to_oup_params`）
- 复现：`gsmm2-metric --mode hp-v2 -q q.fa -t ref.fasta --noMar`
- 实际输出：
  ```text
  panicked at mm2-0.25.6/src/params.rs:243:13:
  assertion failed: self.discard_secondary == false
  ```

原因：`to_oup_params` 无条件先调 `set_discard_secondary(true)`，随后调 `set_discard_multi_align_reads(self.discard_multi_mapping_reads)`。mm2 的 `OupParams` 里这三个互斥开关有断言约束（`discard_multi_align_reads == true` 时要求 `discard_secondary == false` 且 `discard_supplementary == false`），而本项目的 `set_discard_secondary(true)` 是硬编码的，因此只要用户传 `--noMar`，builder 链在 `build_aligner`（aligner.rs:28）里就直接 panic，任何 metric mode 都无法启动。

雪上加霜的是：即使不 panic，`--noMar` 在本管线中也完全没有实现。mm2 只在 `align_worker`（另一条 BAM 输出管线）里检查 `discard_multi_align_reads`，而 gsmm2-metric 走的 `align_single_query_to_targets`（mm2 lib.rs:257-287）根本不看这个字段；`OupParams::modify_aligner` 也只处理 `discard_secondary`。

修复建议：要么移除 `--noMar` 参数（README 未承诺它有独立的实现路径），要么在 `compute_metric` 里自行实现"primary 唯一性过滤"（例如：与 primary query 区间重叠 > 50% 的其他 hit 丢弃），并把 `to_oup_params` 改为不再触发互斥断言的写法。

### Bug 2：`--short-aln` 完全无效（`post_process_param` 是死代码）

- 位置：定义在 `src/cli.rs:29-53`（`Cli::post_process_param`），但全仓库没有任何调用点（grep 确认仅此一处出现）。
- 后果：`--short-aln` 声称会下调 `kmer/wins/min_cnt/min_dp_max/min_chain_score/min_ksw_len`，实际传入后什么都不做，短序列比对会因默认 k=15/w=10 等参数而大量漏比对。README 明确宣传了这个功能，属于文档与行为不符的功能性 bug。

修复建议：在 `main()` 的 `Cli::parse()` 之后调用 `let mut cli = Cli::parse(); cli.post_process_param();`。

### Bug 3：多 fastq 输入时 qname 后缀被追加到「序列」而不是「名字」（依赖 mm2 的 bug）—— ✅ 已修复

- 位置：mm2-0.25.6 `src/lib.rs:213-218`（`query_seq_sender` 的 fq 分支）：
  ```rust
  } else if filename.ends_with("fq") || filename.ends_with("fastq") {
      for mut record in fa_iter {
          if let Some(suffix) = &qname_suffix {
              record.seq.push_str(suffix);   // <- 应该是 record.name
          }
  ```
- 对比 fa 分支（lib.rs:193-195）是正确的 `record.name.push_str(suffix)`。
- 后果：当 `-q a.fq b.fq` 传入多个 fastq 时，每条 read 的序列尾部被污染成 `...<SEQ>__0`，而 qname 没有被重写——序列被破坏（尾部多出 `__0` 会参与比对与 HP 统计），且两个文件的 read 名字无法区分。
- 另外注意 CLI help（cli.rs:131-138）写的是后缀 `___0`，实现是 `__{idx}`（双下划线），文档也不一致。

**修复情况（2026-11-19）**：已确认 mm2 0.26.1 中该 bug 仍未修复，因此在 gsmm2-metric 仓库内绕开依赖修复——在 `src/main.rs` 中重新实现了 `query_seq_sender`（fa / bam 分支行为与 mm2 逐行一致），fq/fastq 分支改为 `record.name.push_str(suffix)`，`metric_entrance` 的调用点改用本地实现。mm2 升级到修复版本后可删除本地实现换回 `mm2::query_seq_sender`。`cargo check --all-targets` / `cargo build --release` 通过；剩余 warning 均为修复前已存在的。双 fastq 输入的端到端验证尚未执行（命令被取消）。

---

## 中严重度

### Bug 4：hp-tr-v3 — motif 起点前一个 ref 位点发生 deletion 时，query 窗口起始错位（called 系统性少算 1 个 unit）

- 位置：`src/metrics/hp_tr_metric_v3.rs:276-280`（q_start 公式）与 `src/metrics/hp_tr_metric_v3.rs:100-103`（`build_ref2qpos` 的 deletion 分支）。

#### 1. `build_ref2qpos` 的写入约定

数组下标 i 对应局部 ref 位点 `ts+i`。三种 aligned_pairs 事件的处理（`src/metrics/hp_tr_metric_v3.rs:86-106`）：

| 事件 | aligned_pairs 产出 | 写入 | 游标更新 |
|---|---|---|---|
| match | `(Some(q), Some(r))` | `ref2qpos[r-ts] = q` | `rcur = r+1`，`qcur = q+1` |
| insertion | `(Some(q), None)` | `ref2qpos[rcur-ts] = q`（写到下一个 ref 位点上，若后面跟 match 会被该 match 覆盖） | `qcur = q+1` |
| deletion | `(None, Some(r))` | `ref2qpos[r-ts] = qcur` | `rcur = r+1` |

关键点：**deletion 分支写入的值不是「上一个 match 的 query 位点」，而是「下一个 match 的 query 位点」**——deletion 分支不推进 qcur，它写入的 qcur 正是下一个 match 的 q。即对被删掉的 ref 位点 r，`ref2qpos[r-ts]` 与「r 之后第一个 match 的 ref 位点」最终取到同一个 query 值。

#### 2. q_start 公式的本意

```rust
let q_start = if start > ts {
    ref2qpos[start - 1 - ts] + 1
} else { ... };
```

这个 `+1` 是为了吸收 `start-1` 与 `start` 之间的 leading insertion（flanking 口径，见源码 272-275 行注释）。两种正常情形都能自洽：

- **无 indel**：`ref2qpos[start-ts] == ref2qpos[start-1-ts] + 1`，`q_start` 恰好落在 motif 首碱基对齐的 query 位点，正确；
- **`start-1` 是 match、其后有 leading insertion**：插入碱基占据 query 区间 `[q_{start-1}+1, q_{start}-1]`（`ref2qpos[start-ts]` 最终被 start 处的 match 覆盖为第一个 match 的 q），`q_start = q_{start-1}+1` 把这些插入全部包进窗口，正确。

#### 3. 失效场景：ref 位点 `start-1` 被 delete

deletion 分支执行 `ref2qpos[start-1-ts] = qcur`，而此时 `qcur` 正是下一个 match（即 ref 位点 `start`）的 query 位点 `q`。于是：

```
ref2qpos[start-1-ts] == ref2qpos[start-ts] == q   （motif 首碱基对齐的 query 位点）
q_start = ref2qpos[start-1-ts] + 1 = q + 1
```

motif 的第一个 query 碱基 `q` 被排除在窗口之外，窗口整体右移 1。

#### 4. 数值例子（`(AC)4`，deletion 在 motif 起点前 1bp）

设 ts=0（局部坐标），motif 区域 `[start=2, end=10)` 内容为 `ACACACAC`，read 在 ref 位点 1 处丢 1 个碱基（motif 区间内容不受影响，仍是完整的 `ACACACAC`）。比对示意（`-` 为 ref 侧 deletion 柱）：

```text
ref pos :  0  1  2  3  4  5  6  7  8  9 10 11
ref     :  G  G  A  C  A  C  A  C  A  C  T  T
           |  -  |  |  |  |  |  |  |  |  |  |
qry pos :  0     1  2  3  4  5  6  7  8  9 10
qry     :  G     A  C  A  C  A  C  A  C  T  T
                 ^^^^^^^^^^^^^^^ motif 区域 [2,10)
```

aligned_pairs 逐步产出与 `build_ref2qpos` 写入：

| 步骤 | aligned_pairs | 写入 | 游标变化 |
|---|---|---|---|
| 1 | `(Some(0), Some(0))` | `ref2qpos[0]=0` | `rcur=1, qcur=1` |
| 2 | `(None, Some(1))` | **`ref2qpos[1]=qcur=1`** | `rcur=2` |
| 3 | `(Some(1), Some(2))` | `ref2qpos[2]=1` | `rcur=3, qcur=2` |
| 4..9 | `(Some(2..8), Some(3..9))` | `ref2qpos[3..9]=2..8` | ... |
| 10 | `(Some(9), Some(10))` | `ref2qpos[10]=9` | ... |
| 11 | `(Some(10), Some(11))` | `ref2qpos[11]=10` | ... |

得到完整映射（注意下标 1、2 处**两个值都是 1**——deletion 分支写的 qcur 与下一个 match 写的 q 相同）：

```text
ref2qpos = [0, 1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
              ^^ ^^
              deletion 分支写入    match 分支写入
```

窗口计算（start=2, end=10, ts=0）：

```text
q_start = ref2qpos[2-1-0] + 1 = ref2qpos[1] + 1 = 2   （期望 1，多偏了 1）
q_end   = ref2qpos[10-0]   - 1 = ref2qpos[10] - 1 = 8 （正确）
```

窗口投影到 query 上的对比：

```text
qry pos :  0     1  2  3  4  5  6  7  8  9 10
qry     :  G     A  C  A  C  A  C  A  C  T  T
期望窗口:        [A  C  A  C  A  C  A  C]        q∈[1,8]  "ACACACAC"
实际窗口:           [C  A  C  A  C  A  C]        q∈[2,8]  "CACACAC"
                    ^^ 首个 A 被切掉了
```

`count_repeats(b"AC", b"CACACAC")` 贪婪逐窗走查：

| i | 窗口 | 结果 |
|---|---|---|
| 0 | `CA` | ≠AC，i→1 |
| 1 | `AC` | ✓ called=1，i→3 |
| 3 | `AC` | ✓ called=2，i→5 |
| 5 | `AC` | ✓ called=3，i→7（越界退出） |

得到 `called = 3`，且 `eq = 6 != 7` → `is_misc = true`。正确答案（窗口 `ACACACAC`）应为 `called = 4`、`pure`。

类似 `AAA-T / AAAAC` 的对照记法（上 ref 下 read，`-` 是 ref 侧 deletion 柱），本例等价于：

```text
GG-ACACACACTT
GG ACACACACTT   <- read 在 motif 起点前少 1 根 ref 柱，窗口起点随之右移 1
```

#### 5. 错位量与 deletion 长度无关；deletion 伸进 motif 内部时同样错

- deletion 长度为 L（ref 位点 `[start-L, start)` 全删）时，这些位点在 deletion 分支里全部写入同一个 `qcur = q0`，`q_start` 恒为 `q0+1`——**错位恒定 1 个碱基**，对 unit_len=2 的 AC motif 是「少 1 个完整 unit 起点」，对 unit_len=1 的 HP motif 直接 called−1。
- deletion 跨过 motif 起点（ref 位点 `start..start+k` 也被删）时同样少 1：此时 `ref2qpos[start-ts]` 也等于「下一个 match（ref 位点 `start+k`）的 query 位点」`q'`，而 `q_start = q'+1` 把 motif 第一个幸存 query 碱基也排除了。即只要 `start-1` 处是 deletion，无论删到哪，窗口起点都统一多偏 1。

#### 6. 对输出的影响

- `called` 系统性 −1（向下偏），且 tag 由 `pure` 翻成 `mixed`（窗口长度 7、eq=6，`eq != query_substr.len()` 恒成立）；
- `Counter` 的 `eq`、`align_span` 同步少算（eq 少 `unit_len`，align_span 少 1）；
- **不会 panic**：整段 motif 被删时 `q_start > q_end` 的 guard（`hp_tr_metric_v3.rs:287-289`）会跳过该区域；危害纯粹是数值错误。带边界 indel 的 read 在真实 ONT 数据里并不罕见，错误行会混入最终 called 分布。

#### 7. q_end 侧复核（确认正确，问题只在 q_start 一侧）

q_end 用 `ref2qpos[end-ts] - 1`，对 `end` 处 deletion 恰好自洽：

- **ref 位点 `end` 被 delete**：deletion 分支写 `ref2qpos[end-ts] = qcur`。若 deletion 紧跟 match（ref 位点 `end-1` 是 motif 末位），`qcur = q_{end-1}+1`，`q_end = q_{end-1}` 正好是 motif 最后一个 query 碱基；若 `end-1` 与 `end` 之间还有 trailing insertion（cigar 形如 `M I D`），insertion 分支把 qcur 一路推进到 `q_ins_last+1`，deletion 再写 `ref2qpos[end-ts] = q_ins_last+1`，`q_end = q_ins_last`——trailing insertion 和 motif 末位碱基都被包进窗口，正确；
- **ref 位点 `end` 是 match**：match 分支覆盖 `ref2qpos[end-ts] = q_{end}`，`q_end = q_{end}-1` 把 `end-1` 与 `end` 之间的 trailing insertion 包进来，正确。

即 q_end 侧在 deletion / insertion 两种顺序组合下都自洽，无需改动。

#### 8. 为什么现有单测没抓到

- `test_build_ref2qpos`（hp_tr_metric_v3.rs:376-413）：只测 `build_ref2qpos` 自身的 deletion 约定（断言 `ref2qpos == [5, 5]`），不经过窗口公式，测不出「约定 × 公式」的组合效应；
- `test_str_called_counting`（hp_tr_metric_v3.rs:488-573）：三个 case（`(AC)6` 扩张、`(AC)3` 收缩、`(AC)4` 带替换）在 motif 边界上都没有 indel，窗口构造始终走正常分支；
- `test_count_repeats`（hp_tr_metric_v3.rs:577-592）：直接测计数函数，与窗口构造完全解耦。

#### 9. 修复建议

在 `build_ref2qpos` 里额外输出一个 `is_deleted: Vec<bool>`（或返回 `(Vec<usize>, Vec<bool>)`），q_start 取：

```rust
let q_start = if start > ts {
    if is_deleted[start - 1 - ts] {
        ref2qpos[start - ts]                      // motif 首碱基（或首个幸存碱基）对齐的 query 位点
    } else {
        ref2qpos[start - 1 - ts] + 1              // 吸收 start-1 与 start 之间的插入
    }
} else { align_info.query_start as usize };
```

该写法对第 5 节两种情形（deletion 止于 `start-1`、deletion 伸进 motif 内部）都能恢复正确窗口。补单测建议：在 `test_str_called_counting` 的 cases 里加两个 case：(a) read 在 ref 位点 173（`(AC)4` 起点 174 的前一位）删 1bp，期望 `called=4`、`pure`；(b) 删 ref 位点 174-175（deletion 伸进 motif 内部），期望值按 minimap2 实际对齐结果确定，断言点都是「窗口起点不再右移 1」。

口径说明：修复后若 cigar 为 `M D I M`（deletion 在 `start-1`、随后紧跟 insertion），`ref2qpos[start-ts]` 最终被 start 处的 match 覆盖，q_start 取第一个 match 的 q，该边界插入被排除在窗口外——这是「motif 首碱基必含」优先的口径选择；与第 2 节无 deletion 时「边界插入必含」的口径不完全对称，如需完全统一需在 README 里明确定义。

**修复情况（2026-11-20）**：已按第 9 节方案修复——`build_ref2qpos` 改为返回 `(Vec<usize>, Vec<bool>)`，每个分支同步写入 `is_deleted`（deletion=true，match/insertion=false，保持「最后写入者生效」语义）；q_start 按 `is_deleted[start-1-ts]` 分支取值。新增回归测试 `test_boundary_deletion_window`（两个 case：deletion 在 motif 起点前 1bp → called=4 pure；deletion 跨入 motif 起点 → called=4 pure），并更新 `test_build_ref2qpos` 断言。注意第二个 case 的实测发现：ref[172]='A' 与 motif 首碱基 ref[174]='A' 相同，deletion 落在 [172,173] 与 [173,174] 产出同一条 read（比对歧义），比对器选择 [172,173] 后 read 的 motif 内容仍是完整 `ACACACAC`，故正确期望为 called=4 pure——修复前后结果分别为 (4,pure) 与 (3,mixed)，仍能有效区分。`metrics::hp_tr_metric_v3` 全部 5 个单测通过；`cargo check --all-targets` 通过，warning 均为修复前已存在的。

### Bug 5：hp-v2 — 窗口两端对 insertion 的处理不对称（前算后不算）—— ✅ 已修复

- 位置：`src/metrics/hp_metric_v2.rs:171-211`。（撰写时的 v1 `src/metrics/hp_metric.rs:101-132` 是同一段拷贝逻辑，该文件已于 2026-11-20 随 v1 整体移除，本条现仅涉及 hp-v2。）

机制：取窗口用 `skip_while(rpos != Some(start-1)).skip(1)` + `take_while((pre_rpos+1) < end)`：

- **起点侧**：跳到 `start-1` 那根柱之后才开始，因此位于 `start-1` 与 `start` 之间的 **leading insertion 会被计入** `cnt`/`is_misc`；
- **终点侧**：`take_while` 在处理完 `rpos == end-1` 的柱后立即停止，因此位于 `end-1` 与 `end` 之间的 **trailing insertion 被排除**。

后果：同一条 read 在 HP 区域边界上的插入，落在起点侧会抬高 `called`（甚至把本应 pure 的位点标成 mixed），落在终点侧则被完全忽略。hp-v2 与 hp-tr-v3（v3 的 flanking 公式两侧都吸收边界插入，见 README 第 73-74 行的描述）对同一事件的口径不一致；hp-v2 自身前后也不一致。

修复建议：把终点侧改成与起点侧对称——多消费 `end` 柱**之前**的 insertion 柱。例如在 `take_while` 之后手动补齐：继续迭代直到遇到 `rpos == Some(end)` 的柱（不含），或者改用 v3 的 `ref2qpos` 式窗口投影统一实现。

**修复情况（2026-11-20）**：已按最小修复方案处理——窗口计数逻辑抽成模块级纯函数 `count_window`（`src/metrics/hp_metric_v2.rs`），起点侧跳过逻辑不变（等价 `skip_while(rpos != Some(start-1)).skip(1)`），终点侧改为遇到首个 `rpos >= end` 的柱即停止，`end-1` 与 `end` 之间的 trailing insertion 柱在停止前计入 `cnt`/`is_misc`/`align_span`，与起点侧及 v3 口径对称；`eq` 仍只统计 match 柱，deletion 柱语义不变。无边界插入时与旧 `take_while` 逻辑完全一致（零回归）。新增 5 个手造 aligned_pairs 回归测试（trailing insertion 计入 / 非 motif 插入置 mixed / leading insertion 行为不变 / 无边界插入等价旧逻辑 / 起点前 deletion+match 跳过），`cargo test` 18 passed / 0 failed；warning 均为修复前已存在的。

### Bug 6：bam2-msa — 参考末端位点全部被 delete 时 `draw_msa` panic；`compute_max_ins` 的 end 边界方向错误

- 位置：`src/msa.rs:283-286`（`major_pos2major_starting_point.get(&rpos).unwrap()`）、`src/msa.rs:89-177`（`compute_max_ins_of_each_ref_position`）。

两个独立问题：

1. **潜在 panic**：`pos2ins` 只在两条路径写入：`rpos_ > start` 时写 `rpos_-1`，以及 break 路径写当前 `rpos_cursor`。若某参考末端位点 P 从未被任何 read match（所有覆盖它的 read 都以 `(None, Some(P))` deletion 结尾且 P 是各自 `reference_end` 的前一位），P 就不会进入 `pos2ins`；而 `draw_msa` 处理这些 read 时 `rpos_cursor` 会等于 P，`get(&P).unwrap()` 直接 panic。
2. **边界方向错误（当前未触发的隐患）**：msa.rs:104-106

   ```rust
   let mut end = rend.unwrap_or(record_ext.reference_end() as i64);
   end = end.max(record_ext.reference_end() as i64);
   ```

   若调用方传入 `rend` 上限，`.max()` 会把上限重新放大回 record 末尾，上界被静默忽略（`start` 一侧用 `.max()` 是对的，`end` 一侧应当用 `.min()`）。当前唯一调用点传 `(Some(0), None)`，所以没有实际触发，但函数签名允许两个 None，属于埋雷。

另附一个观察（非直接 bug）：msa.rs:120-122 的 sentinel 写的是 `aligned_pair_full.last().unwrap()[0]`，变量名叫 `last_ref_pos` 但 `[0]` 取的是 **qpos**。它之所以能工作，是因为 `qpos+1 == query_end` 恰好触发下面 `qpos >= query_end` 的 break 分支，间接把最后一个 ref 位点记进 `pos2ins`；如果末位是 deletion（`[0] == None`）则 sentinel 不会追加。这段逻辑正确性依赖巧合，建议重构成显式「读末尾必须记录最后一根 ref 柱」的逻辑并补注释。

### Bug 7：main.rs — hits 非空但不含 primary 时索引越界 panic（当前被 mm2 的隐式约定兜底）

- 位置：`src/main.rs:175-185`：

```rust
let target_name = if hits.len() > 0 {
    hits.iter().filter(|v| v.is_primary).take(1)
        .map(|v| v.target_name.clone()).collect::<Vec<_>>()[0]   // 无 primary 则越界
        .take().unwrap()
```

- 现状：`mm2::set_primary_supp_alignment_according_2_align_score` 内部有 `assert!(hits.first().is_primary)`，保证只要 hits 非空就存在 primary，所以当下不会触发。但这是对依赖内部实现的隐式依赖：一旦 mm2 升级放宽该约定（或此处过滤条件再叠加别的谓词），就是 worker 线程里的越界 panic（进而 `sender.send(...).unwrap()` 连锁 panic）。
- 修复建议：改为 `if let Some(...)` 取第一个 primary，取不到时回退到第一个 hit 的 target 或空串。

---

## 低严重度 / 健壮性问题

### Bug 8：`get_target_substr` 的 start 收缩循环有理论上的越界 panic

- 位置：`src/metrics/hp_tr_tools.rs:90-96`。当窗口恰好为 `[len-1, len)`（`start = end-1 = len-1`）时，`ori_seq_bytes[start + 1]` 访问 `seq[len]` 越界。需要重对齐窗口为 1bp 且顶到序列末尾才触发，实际概率很低，但建议把条件改成 `while start + 1 < end` 并补防御性检查。

### Bug 9：默认输出文件名对 `rsplit_once(".")` 的 unwrap 依赖

- 位置：`src/cli.rs:171`、`src/msa.rs:26`。
- 两种触发方式：query 路径不含 `.`（panic）；或**目录名**含 `.` 而 query 文件本身不含扩展名（例如 `./data.v2/query` 会在 `./data.v2` 处切开，输出写进错误位置）。建议改用 `Path::extension()`/`with_extension()`。

### Bug 10：bam2-msa 的 Shannon 熵把非 ACGT 碱基静默丢弃

- 位置：`src/msa.rs:292-297`：`base_counter.counter.entry(base).and_modify(|cnt| *cnt += 1)` 只会更新预初始化的 `-/A/C/G/T` 五个 key，read 序列中的 `N`（rust-htslib `seq().as_bytes()` 可以产生）既不报错也不计数，导致熵的分母与 fasta 内容口径不一致。建议用 `or_insert(0)` 或显式报错。

### Bug 11：hp_tr_metric(v1) 会对未使用的 strand 输出全 0 行 —— ✅ 已随 v1 移除而消除

- 原位置：`src/metrics/hp_tr_metric.rs:124-141`（已删除）。

**处置情况（2026-11-20）**：v1 `HpTrMetric` 本就未接入 main（无可达的 mode 分支），连同 `HpMetric` 一起被移除（`Mode::HpTr` / `Mode::Hp` 枚举变体同步删除），本条不再存在。v2/v3 的 counts 是 HashMap，空 map 不产生行，没有此问题。

### Bug 12：无效 mode 的报错文案过时 —— ✅ 已修复

- 原位置：`src/main.rs` 的 `_ =>` 报错分支与 `src/cli.rs` 的 `Mode` 枚举。

**修复情况（2026-11-20）**：随 v1 清理一并处理——`Mode::HpTr` / `Mode::Hp` 变体从枚举中删除，`--mode hp-tr` / `--mode hp` 现在由 clap 在参数解析阶段直接拒绝（报 unknown value 并列出合法值），不再有「接受后落入过时文案报错」的体验问题；`main()` 里不可达的 `_ =>` 分支也一并移除，`match cli.mode` 对全部 4 个变体显式穷尽匹配。

### 其他小问题（不单独列 bug）

- `src/main.rs:81` `assert!(tot_threads >= 10)`：README 有说明，但对小样本调试场景相当不友好（`--threads 4` 直接崩）；且实际线程占用是 `1 + (tot-4) + 1 = tot-2`，预留 4 个是过度保守，可考虑放宽到 `tot_threads >= 4` 并让 `metric_threads = tot_threads - 2`。
- `src/main.rs:206` `BufWriter` 在作用域结束时 drop-flush，flush 错误被静默吞掉；磁盘满时可能"成功"退出而丢尾部数据。建议显式 `writer.flush()?` 并检查。
- `src/main.rs:35`/`207`：`expect(&format!(...))` 是反模式（`expect` 接受 `&str`，用 `.unwrap_or_else(|_| panic!(...))` 或 `expect_owned` 更好），仅风格问题。

---

## 已核查、确认不是问题的点（避免重复排查）

- `aligned_pairs()`（aligned_pairs.rs）对 reverse strand 的 query 游标从 `query_end-1` 倒走：与 mm2 的约定一致（`query_start/query_end` 是原始 read 坐标，mm2 lib.rs:324 注释确认），ins/del/match 的坐标展开均正确。
- eqx 依赖：主管线 `align_single_query_to_targets` 与重比对 `do_align_4_homo` 都显式传了 `MM_F_EQX (0x4000000)`，`aligned_pairs` 里的 `panic!("eqx required")` 不会触发；`MappingExt::identity` 的分母也不含 clip（minimap2-rs 默认 `cigar_clipping=false`）。
- BAM 输入时 `ReadInfo::from_bam_record` 取的是参考方向序列，mm2 对其重新比对，后续 `read_info.seq[qs..qe]` 与 `aligned_pairs` 的 qpos 口径自洽。
- `GlobalData` 的 `OnceCell` 懒初始化在多线程 metric worker 下是安全的（不同 key 各自独立，`get_or_init` 保证单次执行）。
- hp-tr-v2 的 `pre_ins_cnt` 两侧效应（紧贴 motif 起点前的插入计入、末尾后的插入丢失）README 已作为已知行为记录，未列为 bug。
- v2/v3 输出按 motif 字符串聚合、不区分位点，README 已注明是设计如此。

## 修复优先级建议

1. Bug 1（--noMar panic）与 Bug 2（--short-aln 无效）：一行级修复，立刻生效。
2. ~~Bug 4（v3 q_start deletion 错位）~~：已于 2026-11-20 修复（见上文「修复情况」），并补充了边界 deletion 回归测试。
3. ~~Bug 5（hp-v2 前后不对称）~~：已于 2026-11-20 修复（见上文「修复情况」），终点侧已与起点侧及 v3 口径对称。
4. ~~Bug 3~~：已于 2026-11-19 在本仓库内修复（见上文「修复情况」），剩余事项：跑一遍双 fastq 端到端验证；mm2 升级到修复版本后清理本地实现。
5. ~~Bug 11~~ / ~~Bug 12~~：已于 2026-11-20 随 v1 metric 移除一并处理（见上文）。
6. 其余按防御式编程顺手加固。
