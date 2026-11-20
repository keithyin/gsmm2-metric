#!/usr/bin/env python3
"""确定性生成 golden 测试 query reads（从参考 fasta 上切片）。

产出 test_data/queries.fa：
- 正向 reads（覆盖 HP/TR 富集窗口与随机窗口）
- 反向互补 reads
- 带 mismatch / indel 扰动的 reads

固定随机种子，保证迁移前后使用完全相同的输入。
"""
import random

random.seed(20260501)

FASTA = "test_data/ref_Saureus_ATCC25923.m.new.corrected.fasta"
OUT = "test_data/queries.fa"

COMP = str.maketrans("ACGTN", "TGCAN")


def revcomp(seq: str) -> str:
    return seq.translate(COMP)[::-1]


def load_fasta(path):
    name, chunks = None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]
            elif line:
                chunks.append(line)
    return name, "".join(chunks)


def perturb(seq: str, n_mm: int, n_ins: int, n_del: int) -> str:
    """在随机位置引入 mismatch / 1-5bp 的插入与删除。"""
    seq = list(seq)
    # mismatch
    for _ in range(n_mm):
        i = random.randrange(len(seq))
        seq[i] = random.choice([b for b in "ACGT" if b != seq[i]])
    # insertion
    for _ in range(n_ins):
        i = random.randrange(len(seq))
        ln = random.randint(1, 5)
        seq[i:i] = [random.choice("ACGT") for _ in range(ln)]
    # deletion
    for _ in range(n_del):
        i = random.randrange(len(seq))
        ln = min(random.randint(1, 5), len(seq) - i)
        del seq[i:i + ln]
    return "".join(seq)


def main():
    name, ref = load_fasta(FASTA)
    n = len(ref)

    # 采样窗口：均匀铺开 + 固定偏移，可复现
    reads = []
    step = n // 30
    idx = 0
    for pos in range(500_000, n - 2_000, step)[:20]:
        ln = 1500
        end = min(pos + ln, n)
        seq = ref[pos:end]
        reads.append((f"q{idx}_clean_fwd_{pos}", seq, "clean"))
        idx += 1

    # 反向互补（取与上面不同的窗口）
    for k, pos in enumerate(range(600_000, n - 2_000, step)[:6]):
        seq = revcomp(ref[pos:pos + 1500])
        reads.append((f"q{idx}_clean_rev_{pos}", seq, "clean"))
        idx += 1

    # 带扰动的 reads
    for k, pos in enumerate(range(550_000, n - 2_000, step)[:6]):
        seq = perturb(ref[pos:pos + 1500], n_mm=12, n_ins=3, n_del=3)
        reads.append((f"q{idx}_pert_fwd_{pos}", seq, "pert"))
        idx += 1

    with open(OUT, "w") as f:
        for name_, seq, _ in reads:
            f.write(f">{name_}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i + 80] + "\n")

    print(f"wrote {len(reads)} reads to {OUT}")


if __name__ == "__main__":
    main()
