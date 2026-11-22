#!/usr/bin/env bash
# 用指定版本的 gsmm2-metric 二进制生成 BAM / FQ 路径的输出（golden 基线 / 迁移后对比）。
# 用法: scripts/run_bam_modes.sh <binary> <outdir>
set -euo pipefail

BIN=${1:?usage: run_bam_modes.sh <binary> <outdir>}
OUTDIR=${2:?usage: run_bam_modes.sh <binary> <outdir>}

THREADS=12
ARGS=(--threads $THREADS --preset map-ont -m 2 -M 5 -o 2,24 -e 1,0)
BAM=test_data/queries_bam/fixture.bam

mkdir -p "$OUTDIR"
cd "$(dirname "$0")/.."   # project root

# metric 模式 ×6
for mode in hp-tr-v2 hp-tr-v3 hp-v2; do
    for ra in "" "--refAnchored"; do
        if [ -z "$ra" ]; then
            name="$mode"
        else
            name="${mode}-refanchored"
        fi
        "$BIN" --mode "$mode" "${ARGS[@]}" $ra \
            -q "$BAM" \
            -t test_data/ref_Saureus_ATCC25923.m.new.corrected.fasta \
            --out "$OUTDIR/${name}.csv"
        echo "done: $name"
    done
done

# aligned（全比对质量 metric；BAM 输入路径；不做 --refAnchored 交叉，理由同 run_all_modes.sh）
"$BIN" --mode aligned "${ARGS[@]}" \
    -q "$BAM" \
    -t test_data/ref_Saureus_ATCC25923.m.new.corrected.fasta \
    --out "$OUTDIR/aligned.csv"
echo "done: aligned"

# np/rq 过滤
"$BIN" --mode hp-v2 "${ARGS[@]}" \
    --np-range 1:2 --rq-range 0.9:1.0 \
    -q "$BAM" \
    -t test_data/ref_Saureus_ATCC25923.m.new.corrected.fasta \
    --out "$OUTDIR/hp-v2-filtered.csv"
echo "done: hp-v2-filtered"

# bam2-msa（生成 <out-base>.fasta 与 .fasta.shannon.txt）
"$BIN" --mode bam2-msa --threads 12 \
    -q "$BAM" \
    -t test_data/ref_Saureus_ATCC25923.m.new.corrected.fasta \
    --out "$OUTDIR/bam2-msa.csv"
echo "done: bam2-msa"

# fastq 输入（覆盖 noodles-fastq 重写路径）
"$BIN" --mode hp-v2 "${ARGS[@]}" \
    -q test_data/queries_fq/fixture.fq \
    -t test_data/ref_Saureus_ATCC25923.m.new.corrected.fasta \
    --out "$OUTDIR/fq-hp-v2.csv"
echo "done: fq-hp-v2"
