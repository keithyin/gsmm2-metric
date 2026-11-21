#!/usr/bin/env python3
"""从 test_data/queries.fa 取前 N 条 reads 生成 fastq fixture(qual 全 'I' = phred40)。

用法: python3 scripts/gen_fq_fixture.py [n_reads] [out.fq]
"""
import sys

FA = "test_data/queries.fa"
N_DEFAULT = 5
PHRED40 = "I"


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else N_DEFAULT
    out = sys.argv[2] if len(sys.argv) > 2 else "test_data/queries_fq/fixture.fq"

    records = []
    name = None
    with open(FA) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name is not None and len(records) >= n:
                    break
                name = line[1:]
                records.append([name, ""])
            elif name is not None:
                records[-1][1] += line
    records = records[:n]

    with open(out, "w") as f:
        for name, seq in records:
            f.write(f"@{name}\n{seq}\n+\n{PHRED40 * len(seq)}\n")

    print(f"wrote {out}: {len(records)} reads")


if __name__ == "__main__":
    main()
