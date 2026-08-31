# gsmm2-metric

对 ONT 测序数据（fastq / fasta / bam）按 **read 粒度**统计 homopolymer(HP) 与 tandem repeat(TR) 位点上的比对特征与重复数分布。

不同的 `--mode` 决定「拿到重对齐结果之后统计什么」，比对、motif 发现、重对齐这几步是**共用**的。本文档重点说明各 mode 的统计逻辑差异。

## 构建与运行

```bash
make build            # cargo build --release -> target/release/gsmm2-metric

gsmm2-metric --mode hp-tr-v3 \
    -q query.fastq -t target.fasta \
    --threads 12 --preset map-ont \
    -m 2 -M 5 -o 2,24 -e 1,0 \
    --out result.HpTrV3-fact.csv
```

- `--threads` 至少 10：其中 `threads - 4` 个线程算 metric，其余用于读入、比对、落盘。
- 默认输出文件名 `{query 去掉扩展名}.{Mode:?}-fact.csv`，同时写一个同名 `.log`。**文件扩展名是 csv，实际分隔符是制表符**（历史原因），下游按 TSV 解析。
- `--short-aln`：短 query / 短 target 场景下调低 kmer 与 min_cnt 等阈值（见 [src/cli.rs](src/cli.rs)）。
- `--noSupp` / `--noMar`：丢弃 supplementary / 多比对 read。
- `--refAnchored`：reverse 比对先把 read 子串和 target 子串都反向互补再重对齐，使统计统一锚定到 reference 方向。**该开关同时影响 hp-tr-v2 / hp-tr-v3 / hp-v2 三个 mode 的取窗口逻辑**。

## Mode 一览

| mode | metric 实现 | 统计对象 | 输出列 |
|---|---|---|---|
| `hp-tr` | [hp_tr_metric.rs](src/metrics/hp_tr_metric.rs) | 未接入 `main` 分发，会直接报错退出 | — |
| `hp-tr-v2` | [hp_tr_metric_v2.rs](src/metrics/hp_tr_metric_v2.rs) | **比对算子谱**：motif 区域内的 eq/diff/ins/del | `qname rname motif eq diff ins del` |
| `hp-tr-v3` | [hp_tr_metric_v3.rs](src/metrics/hp_tr_metric_v3.rs) | **重复数分布**：ref repeats vs called repeats | `qname rname motif tag ref called num eq alignspan` |
| `hp` | [hp_metric.rs](src/metrics/hp_metric.rs) | 未接入 `main` 分发，会直接报错退出 | — |
| `hp-v2` | [hp_metric_v2.rs](src/metrics/hp_metric_v2.rs) | **HP 重复数分布**（仅 unit 长度 1） | `qname rname motif tag called num eq alignspan` |
| `bam2-msa` | [msa.rs](src/msa.rs) | 不走 `TMetric`，由 bam 生成 MSA fasta + 逐列 Shannon 熵 | `.fasta` / `.shannon.txt` |

## 共用流水线

所有 metric mode 都走 [main.rs](src/main.rs) 的 `metric_entrance::<M>`：

1. read 比对到 target（minimap2，`--preset` / `-m -M -o -e` 控制打分）。
2. 保留 primary（以及未开启 `--noSupp` 时的 supplementary），只保留与 primary **同一个 target** 的 hit，并且 `identity > 0.75`（[main.rs:189-192](src/main.rs#L189-L192)）。
3. 每个 hit 依次执行下面的「mode 无关前缀」：
   1. **motif 边界收缩**：如果比对端点落在某个 HP/TR 区域内部，就把端点收缩到该区域的边界，避免 flank 干扰重对齐。查的是 `GlobalDataKey::TargetRegion2Motif4HpTr`（按 reference 位点建索引，见 [global_data.rs:100-111](src/global_data.rs#L100-L111)）。仅 `hp-tr-v2` / `hp-tr-v3` 做这一步，`hp-v2` 不做。
   2. `get_target_substr()`：取 target 子串；如果切片边界落在同种碱基内部，再往内收缩 1~n 个碱基，保证窗口不切断 homopolymer（[hp_tr_tools.rs:16-49](src/metrics/hp_tr_tools.rs#L16-L49)）。
   3. `--refAnchored` 且 hit 是 reverse 时，read 子串和 target 子串一起反向互补。
   4. `do_align_4_homo()`：以 target 子串为 reference 重新做一次 Smith-Waterman 式精比对（固定参数 `a=4 b=10 q=4 q2=24 e=2 e2=1 k=11 w=1`，只取 primary）。重对齐失败或仍然是 reverse 的结果直接跳过并打 warn。
   5. 在 **target 子串**（即 reference 侧）上跑 `single_seq_hp_tr_finder()` 找 motif。正则表 `HP_TR_REG` = unit 长度 1~4、最少 3 次重复（[hp_tr_metric_v2.rs:37-46](src/metrics/hp_tr_metric_v2.rs#L37-L46)），`hp-v2` 只用 unit 长度 1。`hp_tr_finder` 内部会做 `dedup_overlaps`：一个区域只保留最长的那条注释（`(AT)23` 会吃掉同区域的 `(TA)22`、`(TATA)11`）。

**分歧从第 6 步开始**：如何把「motif 区间」映射到 query，以及统计什么。

## hp-tr-v2：比对算子谱（per-base）

[hp_tr_metric_v2.rs:180-213](src/metrics/hp_tr_metric_v2.rs#L180-L213)

1. 对重对齐结果调 `aligned_pairs()`，得到逐碱基的 `(qpos, rpos, AlignOp)` 序列。
2. 把 `region2motif.flatten()` 展开成「ref 位点 → 覆盖它的 motif 列表」，对每个有 `rpos` 的比对柱查表。
3. 把该柱的 `AlignOp`（`Equal`/`Diff`/`Ins`/`Del`）累加到**每一个**覆盖该位点的 motif 上；过滤条件 `motif.0.0 >= target_start && motif.0.1 <= target_end`，即 motif 必须完整落在比对范围内。
4. 额外维护 `pre_ins_cnt`：连续插入串（`rpos == None`）的长度，在遇到下一个有 `rpos` 的柱时，整体作为 `Ins` 记给覆盖该柱的 motif，随后清零。

结果：每个 motif 一行，四个计数。**它不解析 motif 字符串**，因此不知道 reference 上有几个 repeat，也无法直接回答「这条 read 测到几重复」。

> `pre_ins_cnt` 的两端效应：紧贴 motif 起点之前的插入会被算进该 motif；而 motif 末尾之后、下一个 motif 柱之前的插入串，如果后面没有再出现被覆盖的柱，就不会被计入。

## hp-tr-v3：重复数分布（ref vs called）

[hp_tr_metric_v3.rs:80-108](src/metrics/hp_tr_metric_v3.rs#L80-L108)、[hp_tr_metric_v3.rs:247-296](src/metrics/hp_tr_metric_v3.rs#L247-L296)

1. `build_ref2qpos()`：遍历 `aligned_pairs()`，构造 `ref 局部坐标 -> query 坐标` 的数组（长度固定为 `target_end - target_start`）。indel 规则：
   - **deletion**（ref 有 query 无）→ 指向下一个未消费的 query 位点；
   - **insertion**（query 有 ref 无）→ 回落到当前 ref 游标 `rcur`，若后续有 match 会被覆盖；
   - query 游标从 `align_info.query_start()` 起算，所以**前导 deletion 映射到第一个 match 的 query 位点，而不是 0**（这条不变量有单测 [test_build_ref2qpos](src/metrics/hp_tr_metric_v3.rs#L370-L408)）。
2. 逐条 motif 区间 `[start, end)`：跳过不满足 `start >= ts && end <= te` 的。
3. 用**两侧 flanking 公式**把 ref 区间投影成 query 窗口：
   `q_start = ref2qpos[start-1] + 1`、`q_end = ref2qpos[end] - 1`（motif 顶到比对首/末端时回落到 `query_start` / `query_end-1`）。这样 motif 边界两侧的插入/删除会被**吸收进窗口**，而不是被漏掉；`q_start > q_end`（该区域在 query 侧整个被删掉）则跳过。
4. `parse_motif()` 把 `(AC)4` 拆成 unit=`AC`、`ref_repeats=4`。
5. `count_repeats()` 在 query 窗口里**贪婪、非重叠**地数完整 unit：从左往右滑，命中一个 unit 就跳 `unit_len`，否则滑一格。产出
   - `called` = 完整匹配 unit 数；
   - `eq` = `called × unit_len`（被完整 unit 覆盖的碱基数）；
   - `is_misc` = `eq != 窗口长度`，即存在前导 / 尾随 / 夹缝碱基 → 标 `mixed`，否则 `pure`。
6. 聚合进 `Counter`：`counts[pure|mixed][called] += 1`，同时累加 `eq`、`align_span`（= query 窗口长度）。

结果：每个 motif 按 `(tag, called)` 展开成多行，附带 `ref`（reference 重复数）。这是**基因型/等位基因长度**视角的指标。

## hp-v2：HP 专用重复数分布

[hp_metric_v2.rs:157-238](src/metrics/hp_metric_v2.rs#L157-L238)

输出格式与 `hp-tr-v3` 几乎一致（少一个 `ref` 列，因为 HP motif 名本身就带 reference 长度，且正则表只含 unit 长度 1）。差别在窗口取法和计数方式：

| | hp-tr-v3 | hp-v2 |
|---|---|---|
| motif 边界收缩（第 3.1 步） | 做 | **不做**，直接用原始比对范围 |
| 取 query 窗口 | `ref2qpos` 数组 + flanking 公式 | 直接在 `aligned_pairs()` 上用 `take_while`/`skip_while` 截取两端之间的比对柱 |
| 计数方式 | 非重叠贪婪匹配整个 unit | 逐 query 碱基与 HP 的 `target_base` 比较，相等则 `cnt += 1` |
| `eq` 含义 | `called × unit_len` | 有 `rpos` 的比对柱中碱基相等的数量（前导/尾随插入不计入 eq） |
| `is_misc` | `eq != 窗口长度` | 窗口内出现过任何不等于 `target_base` 的碱基 |
| 异常输出 | — | `abs(cnt - expected_cnt) >= 10` 时，把对齐后的 target/query 两条串打到 `.log`，便于肉眼查错 |
| motif 字符串 | `parse_motif()` 通用解析 | 位置硬编码：`motif.as_bytes()[1]` 取碱基、`motif[3..]` 取重复数，只对 `(X)N` 形式成立 |

因为 unit 长度为 1，两者对纯 HP 的 `called` 结果应当一致；差异集中在 `eq` 的定义和窗口边界处理上。

## 列名语义对照（同名不同义）

| 列 | hp-tr-v2 | hp-tr-v3 | hp-v2 |
|---|---|---|---|
| `eq` | 比对柱中 `Equal` 的碱基数 | `called × unit_len` 累加（被完整 unit 覆盖的碱基） | 比对柱中碱基与 HP base 相等的数量 |
| `diff` / `ins` / `del` | 比对柱中 `Diff` / `Ins` / `Del` 的碱基数（`ins` 含 `pre_ins_cnt`） | 无此列 | 无此列 |
| `alignspan` | 无此列 | query 窗口长度累加 | 截取的比对柱数量 |
| `tag` | 无 | `pure` / `mixed` | `pure` / `mixed` |
| `ref` | 无 | reference 上的重复数 | 无（在 motif 名里） |
| `called` | 无 | 该 read 实测到的完整 unit 数 | 该 read 实测到的 HP 碱基数 |
| `num` | 无 | **该 read 内**同名 motif 出现该 `(tag, called)` 组合的次数 | 同左 |

## 实测示例

同一位点 `(AC)4`（reference `CP009361.1:174-182 = ACACACAC`），构造三种 read：与 reference 一致 / 扩成 `(AC)6` / 缩成 `(AC)3`。

`hp-tr-v2`（比对算子谱）：

```text
qname	rname	motif	eq	diff	ins	del
read_ref	CP009361.1	(AC)4	8	0	0	0
read_ac6	CP009361.1	(AC)4	8	0	4	0
read_ac3	CP009361.1	(AC)4	6	0	0	2
```

`hp-tr-v3`（重复数分布）：

```text
qname	rname	motif	tag	ref	called	num	eq	alignspan
read_ref	CP009361.1	(AC)4	pure	4	4	1	8	8
read_ac6	CP009361.1	(AC)4	pure	4	6	1	12	12
read_ac3	CP009361.1	(AC)4	pure	4	3	1	6	6
```

读法：v2 说「这里多了 4 个碱基的插入 / 少了 2 个碱基的删除」，v3 说「reference 是 4 重复，read 测到 6 / 3 重复」。二者互补：v2 适合评估比对质量与错误构成，v3 适合看 STR 的等位基因长度分布。

## 其它注意点

- **聚合 key 只有 motif 字符串**，输出里没有位点坐标列。同一 read 内不同位置的**同名** motif 会被合并到一行累加（`hp-tr-v2` 的 `eq` 等直接相加；`hp-tr-v3` / `hp-v2` 落进同一张直方图）。例如上面 `read_ref` 的 `(T)3` 在 v2 里 `eq=30`，其实是该 read 上 10 个 `(T)3` 位点各贡献 3 个碱基的总和——同一位点在 v3 输出里 `num=10`，正好可以对照验证。
- **每个 read 一个 metric 实例**，行内 `num` 只统计这一条 read 内的出现次数；跨 read 的重复数分布需要在下游按 `(rname, motif, tag, called)` 再做一次汇总。
- `hp-tr` / `hp` 两个 mode 的实现（按 strand 分 `[fwd, rev]` 两组计数、不做重对齐）仍保留在仓库里，但 [main.rs](src/main.rs) 的 `match` 没有分发它们，执行会报 `not a valid mode`（该报错文案也漏写了 `hp-tr-v3`）。
- `bam2-msa` 的统计口径完全不同：先用 `compute_max_ins_of_each_ref_position()` 求出每个 ref 位点上**所有 read 的最大插入长度**，据此把 ref 坐标展开成一条带空隙的「major 坐标轨」，再把每条 read 摆进 MSA；同时逐列输出 A/C/G/T/`-` 的 Shannon 熵。
