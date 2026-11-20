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
