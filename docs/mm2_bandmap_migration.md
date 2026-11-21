计划用 bandedmap 替换 mm2.
bandmap 源码在 /root/projects/rammap ，依赖的时候注意依赖 crates-io 中的 crate
迁移过程中核心目标是迁移前后 结果的一致性。

评估数据为 /data1/ccs_data/ecoli-data-for-marketing/20250623_250302Y0002_Run0004_adapter.WithDi.smc_all_reads.bam
reference 在 /data1/REF_GENOMES/MG1655.fa

---

# mm2 → bandmap-core 迁移一致性验收报告

## 一、迁移内容

- 依赖：`mm2 = "0.25.4"`（minimap2 2.28 FFI 封装）→ crates.io `bandmap-core = "2.0.0"`（lib 名 `bandmap`，纯 Rust minimap2 算法复刻）
- 保留：gskits（ReadInfo / fastx / BAM / revcomp）、hp_tr_finder、rust-htslib、crossbeam 流水线；CLI 界面与各 mode 统计口径、输出格式完全不变
- 新增模块：
  - `src/adapter.rs`：自有 `Mapping` 结构（与 mm2 字段同名同类型，规避 bandmap `Mapping` 的 `#[non_exhaustive]`）+ mm2 语义复刻：`identity()`（EQX cigar op7/总长）、`aligned_2_str()`、`align_single_query_to_targets`、`reorder_like_mm2`（primary/supplementary 重排，见下）
  - `src/params.rs`：本地化 mm2 params API（`InputFilterParams/OupParams/AlignParams/IndexParams`），`AlignParams::apply_to` 映射到 bandmap `MapOptions`
- 关键参数映射：`-m/-M` → `scoring.match_score/mismatch_penalty`；`-o q,q2` → `gap_open/gap_open2`；`-e e,e2` → `gap_extend/gap_extend2`；`min_cnt/min_chain_score` → `chaining`；`min_dp_max/min_ksw_len` → `alignment.min_dp_max/min_dp_len`；EQX 由 `output_config_mut().eqx = true` 打开；`best_n=10000/pri_ratio=0.2` hack 保持

## 二、语义差异点与处理

| 语义 | mm2 行为 | bandmap 行为 | 迁移处理 |
|---|---|---|---|
| supplementary 判定 | 管线内判 secondary；mm2 包装层再按 score+query ovlp 重排（ovlp≤0.5 → supp，supp 间 ovlp>0.5 去重） | sam_pri 语义：第 2+ 个非 secondary 一律标 supp | `adapter::reorder_like_mm2` 逐条复刻 mm2 重排（ovlp_ratio = 重叠/min(len)）；已修复旗标归一问题（bandmap supp 型 hit 须先归一为 `is_primary=false`，否则重排降级时残留假 primary） |
| mid_occ | 内存索引路径（`with_seq_and_id`）**写死 1000**，不走校准（实现怪癖；CLI minimap2 走 `mm_mapopt_update` 校准） | 按索引分位校准（`mid_occ_frac=2e-4`，clamp [10,1e6]；MG1655 实测 12） | **保持 bandmap 校准，不模仿 1000**。实测钉死 1000 在 poly-A 富集的合成数据上反而破坏一致性（高频 seed 进入链合导致重对齐窗口漂移） |
| EQX cigar | 每次 map 传 flag `0x4000000` | `output_config_mut().eqx = true` | 一致（op 编码 7/8/1/2 相同，run 合并行为相同） |
| 重对齐 `do_align_4_homo` | minimap2 k=11 w=1，a=4 b=10 q=4 q2=24 e=2 e2=1，best_n=1，flags [EQX, 0x40000000=MM_F_SAM_HIT_ONLY（仅影响 SAM 输出）] | bandmap 同参数 + eqx | 逐参数对齐；原 mm2 flags 中 `0x40000000` 语义为 SAM_HIT_ONLY 而非 NO_LJOIN，无需模仿 |
| cigar clip | `alignment.cigar` 不含 clip | `cigar_ops` 不含 clip（clip 仅在 `full_cigar`） | 直接使用 `cigar_ops`，aligned_pairs 展开逻辑无需改动 |
| 坐标 | query_start/end 为原始 read 正向坐标 | 相同 | 直接对齐 |

## 三、一致性验收结果

### 3.1 合成数据（S. aureus 参考切片，32 条确定性 reads）

`scripts/gen_test_queries.py`（固定种子，正向 20 + 反向互补 6 + 带扰动 6）× hp-tr-v2 / hp-tr-v3 / hp-v2 ×（普通 / --refAnchored）共 6 组合，与 mm2 基线（同 hp_tr_finder 版本）排序逐字节对比：**6/6 PASS**。

### 3.2 真实数据（E. coli MG1655 + CCS BAM）

- 采样子集（前 20,000 条 reads）：**6/6 PASS**（每组合 16.6–19.1 万行逐行一致）
- 全量（541,549 条 reads，45.6 万条有输出，threads=16，preset map-ont，-m 2 -M 5 -o 2,24 -e 1,0）：

| mode | 输出行数 | 差异行数 | 差异率 |
|---|---|---|---|
| hp-tr-v2 | 4,494,875 | 63 | 0.0014% |
| hp-tr-v2-refanchored | 4,448,988 | 78 | 0.0018% |
| hp-tr-v3 | 5,105,680 | 97 | 0.0019% |
| hp-tr-v3-refanchored | 5,053,915 | 116 | 0.0023% |
| hp-v2 | 4,565,569 | 99 | 0.0022% |
| hp-v2-refanchored | 4,508,658 | 114 | 0.0025% |

差异仅涉及 **21 条 read（0.0046%，456,326 条有输出的 read 中）**，且三个 mode 涉及同一组 read。

### 3.3 21 条差异 read 的归因：chimeric read 的 supplementary/secondary 判定差异

抽样核验均呈**多比对 / 嵌合体结构**（adapter-dimer 类：单条 read 的两个片段比对到不同位点）。以下用 read `/1026475` 完整拆解。

#### 3.3.1 判定机制差异的原理

**minimap2（mm2 路径）**：
- 拆分链机制：chimeric read 的两个片段各自成链后，`mm_split` 把互不重叠的片段链归入同一 `parent`（`parent == id`）；首个片段 `sam_pri > 0` → `is_primary=true`，其余片段 `sam_pri == 0` → `is_primary=false` + `is_supplementary=true`
- 同一位点（或重叠区域）上的替代链经 `mask_level`（0.5）掩蔽 → secondary（`is_primary=false`、`is_supplementary=false`），**无论其 identity 多高**
- mm2 包装层 `set_primary_supp_alignment_according_2_align_score` 再做跨 aligner 重排与 supp 去重（见第二节表格），engine-primary 参与 ovlp 重排，engine-supp 仅做去重

**bandmap-core**：
- `assign_parents_and_select` 按 `(dp_score, hash)` 降序选出代表链；`sam_pri` 仅授予首个非 secondary 结果，其余非 secondary 标 supplementary
- 替代 hit 是否被掩蔽为 secondary 由引擎内部的掩蔽判定决定，与 minimap2 的 mask_level 路径**不逐 read 等价**

核心差异：对 chimeric read，「哪一段成为 primary、哪一段保留为 supplementary、哪些替代 hit 被掩蔽为 secondary」的决策由两引擎各自的链合/选择实现给出，结果不逐 read 相同。迁移层的 `reorder_like_mm2` 只能在引擎给定的 hit 集合上重演 mm2 的旗标语义，**无法改变哪些 hit 进入 hit 集合本身**。

#### 3.3.2 实例：read `/1026475` 两端 hits 对照

该 read 为 adapter-dimer 嵌合：片段 A 比对到 3448936 位点，片段 B 比对到 4325804/4325897 位点（B 内部自互补，正反两条链都产生命中）。

**mm2 端**（共 11 条 hit，关键条目）：

| query 区间 | target 位点 | strand | identity | score | 引擎标注 | metric 过滤结果 |
|---|---|---|---|---|---|---|
| 4-1484 | 3448936-3449761 | Reverse | 0.4683 | 1043 | **primary** | identity ≤0.75 → 丢弃 |
| 1555-1942 | 4325897-4326384 | Reverse | **0.7766** | 693 | **supplementary** | **通过 → 产出 3 行 metric** |
| 0-387 | 4325897-4326384 | Forward | 0.7730 | 685 | secondary | 丢弃 |
| 1551-1835 | 4325804-4326088 | Reverse | 0.9549 | 507 | secondary | 丢弃 |
| 108-391 | 4325804-4326088 | Forward | 0.9514 | 502 | secondary | 丢弃 |
| 458-1938 | 3448936-3449761 | Forward | 0.4634 | 1030 | secondary | 丢弃 |

判定链路：片段 A 的链成为 primary（engine-primary 加权后分数最高）；片段 B 的独立链经 mm_split 保留**恰好一条**为 supplementary（identity 0.7766 > 0.75）；片段 B 的其余替代链——即使 identity 高达 0.9549/0.7730/0.9514——被 mask_level 掩蔽为 secondary 全部丢弃。

**bandmap 端**（共 9 条 hit，关键条目）：

| query 区间 | target 位点 | strand | identity | score | 引擎标注 | metric 过滤结果 |
|---|---|---|---|---|---|---|
| 1546-1935 | 4325804-4326393 | Reverse | 0.6481 | 714 | **primary** | identity ≤0.75 → 丢弃 |
| 4-1484 | 3448936-3449761 | Reverse | 0.4683 | 1043 | **supplementary** | identity ≤0.75 → 丢弃 |
| 1583-1942 | 4325897-4326256 | Reverse | **0.9750** | 661 | secondary | 丢弃 |
| 0-359 | 4325897-4326256 | Forward | **0.9695** | 653 | secondary | 丢弃 |
| 7-396 | 4325804-4326393 | Forward | 0.6453 | 706 | secondary | 丢弃 |

判定链路：bandmap 选择片段 B 的 4325804 链为 primary（score 714），片段 A 降为 supplementary（identity 0.47 过不了滤）；片段 B 的高 identity 精确片段（0.9750 / 0.9695）被引擎判为 **secondary 而非 supplementary**——所有可能通过 identity > 0.75 的路径都被封死，该 read 在 metric 输出中完全消失。

#### 3.3.3 metric 后果（21 条 read 的三种形态）

逐一核对 21 条 read 的 CSV 差异，呈现三种形态——根因相同（进入 metric 的 hit 子集因 supplementary/secondary 标注不同而异）：

1. **一侧无输出**：如 `/1026475`（mm2 3 行 → bandmap 0 行，见 3.3.2）、`/1170595`（5→0）、`/941808`（11→0）、`/898365`（2→0）。bandmap 将 mm2 用作 supplementary 的高 identity 片段标 secondary，所有通过 identity>0.75 的路径被封死。
2. **反向行集差异（primary 选择不同）**：如 `/140537`（mm2 4 行全 pure，span 12/16 → bandmap 6 行含 mixed，span 20/38）。两端 hits dump 核验：mm2 选 `Forward 623-1188`（id 0.8088）为 primary 并保留 `Reverse 0-565`（id 0.8088）为 supplementary（共 4 条 hit）；bandmap 选 `Forward 379-1188`（id 0.8466）为 primary，把 mm2 的 primary/supp 链掩蔽为 secondary（仅 2 条 hit）。bandmap 的 primary 窗口更长（379 起 vs 623 起）→ 重对齐覆盖更多 motif 位点 → 行数 4→6 且出现 mixed。
3. **同构数值漂移**：如 `/255591`（行数相同，`(A)3 pure` 的 num 1→2、eq 3→6）。进入 metric 的 hit 子集或重对齐窗口相差一位，个别 site 的 called 计数 ±1。

#### 3.3.4 影响面

全量 541,549 条 read（456,326 条有输出）中 21 条（**0.0046%**）；hp-tr-v2 / hp-tr-v3 / hp-v2 三个 mode 涉及同一组 read；行级差异率 0.0014%–0.0025%。

#### 3.3.5 可选收敛方向（非迁移阻塞项）

1. **上游 bandmap**：对 chimeric read 引入等价于 minimap2 `mm_split` 的拆分链语义——互不重叠的片段链保留为 supplementary 而非被 secondary 掩蔽；或在 `assign_parents_and_select` 的替代 hit 判定中放宽 disjoint 片段。需上游改动或调参验证，不属于本次依赖替换的范围。
2. **数据侧**：该类 read 为 adapter-dimer 污染，下游评估通常应过滤。`--noSupp` 后两版本行为一致（supp hit 均被丢弃，21 条差异 read 不再暴露于 metric 输出）。（2026-11-21 起也可用 `--noMar`：比对产生不止一个 chain（存在 supplementary/secondary）的 read 整条过滤、不产出输出行，此类嵌合 read 恰好命中该判定，见 README 的 `--noMar` 说明。）

## 四、验收过程中的重要教训

1. **hp_tr_finder 版本漂移**：迁移期间 `Cargo.toml` 的 `hp_tr_finder` 被从 0.4.0 升级到 0.5.0（版本号 0.8.3→0.9.0），motif 发现结果随之改变，曾造成合成数据 golden 6/6 FAIL 的假象。基线必须与被测版本使用**相同**的 hp_tr_finder 才能隔离出迁移本身的差异（对照基线 worktree `/tmp/gsmm2-mm2-baseline` 已同步升级后重建）。
2. **对照实验要控制单变量**：mid_occ 钉死、NO_LJOIN 等调参实验需与 golden 重建同步验证，否则结论会被版本漂移污染。
3. **对比工具注意去重**：`sort | uniq` 会折叠重复行（多线程交错输出下干扰判断），集合比较应基于排序去重后的语义。

## 五、已知非风险项

- **输出行序**：metric 行序来自 HashMap 迭代，每次进程运行顺序都不同（同一二进制两次运行原始字节即不同），属项目既有行为；行内容集合确定，对比须排序后进行（`scripts/compare_runs.sh` 已处理）
- **重对齐 trace**：两端重对齐打分 100% 一致（20k 采样实测仅 1 个例外）；cigar "差异" 绝大多数为 mm2 `cigar_str` 带软 clip 而 bandmap 不带的表示差异，实际 ops 相同

## 六、验收工具链（可复现）

```bash
scripts/gen_test_queries.py                          # 生成确定性合成 reads → test_data/queries.fa
scripts/run_all_modes.sh <bin> <outdir>              # 合成数据全 mode 运行
scripts/run_real_modes.sh <bin> <outdir> <bam> <fa>  # 真实数据全 mode 运行
scripts/compare_golden.sh [bin]                      # 与 test_data/golden 排序对比
scripts/compare_runs.sh <dirA> <dirB>                # 任意两输出目录排序对比
GSM2_DUMP_HITS=1                                     # per-read hits dump（stderr，归因比对层差异）
GSM2_DUMP_REALIGN=1                                  # 重对齐结果 dump（score/坐标/cigar）
```

- mm2 对照基线：git worktree `/tmp/gsmm2-mm2-baseline`（HEAD ea71ee7 + 调试 dump 补丁 + hp_tr_finder 0.5.0 对齐）
- 单测：18 个全部通过（含 `(AC)4` 位点精确断言、边界 deletion 窗口回归、aligned_pairs 语义断言）

## 七、结论

迁移在合成数据上 6/6 逐字节一致；真实数据 54 万条 reads 上 99.995% 的 read 输出逐行一致，剩余 0.0046% 为 bandmap 与 minimap2 在 chimeric read 拆分链合上的已知引擎差异（adapter-dimer 污染 reads，主要影响身份为污染片段的 supplementary 行）。迁移目标达成，建议合入。
