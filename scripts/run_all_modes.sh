#!/usr/bin/env bash
# 用指定版本的 gsmm2-metric 二进制生成各 mode 的输出（golden 基线 / 迁移后对比）。
# 用法: scripts/run_all_modes.sh <binary> <outdir>
set -euo pipefail

BIN=${1:?usage: run_all_modes.sh <binary> <outdir>}
OUTDIR=${2:?usage: run_all_modes.sh <binary> <outdir>}

THREADS=12
ARGS=(--threads $THREADS --preset map-ont -m 2 -M 5 -o 2,24 -e 1,0)

mkdir -p "$OUTDIR"
cd "$(dirname "$0")/.."   # project root

for mode in hp-tr-v2 hp-tr-v3 hp-v2; do
    for ra in "" "--refAnchored"; do
        if [ -z "$ra" ]; then
            name="$mode"
        else
            name="${mode}-refanchored"
        fi
        "$BIN" --mode "$mode" "${ARGS[@]}" $ra \
            -q test_data/queries.fa \
            -t test_data/ref_Saureus_ATCC25923.m.new.corrected.fasta \
            --out "$OUTDIR/${name}.csv"
        echo "done: $name"
    done
done

# aligned（全比对质量 metric）：不做 --refAnchored 交叉——该 flag 对其无效
# （raw mapping 不重对齐），且 mm2 原版无此概念，golden 与 mm2
# gsmm2-aligned-metric 输出对齐。
"$BIN" --mode aligned "${ARGS[@]}" \
    -q test_data/queries.fa \
    -t test_data/ref_Saureus_ATCC25923.m.new.corrected.fasta \
    --out "$OUTDIR/aligned.csv"
echo "done: aligned"
