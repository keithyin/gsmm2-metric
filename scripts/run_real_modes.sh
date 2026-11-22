#!/usr/bin/env bash
# 用指定版本的 gsmm2-metric 二进制在真实数据上生成各 mode 输出。
# 用法: scripts/run_real_modes.sh <binary> <outdir> <query.bam> <target.fa> [threads]
set -euo pipefail

BIN=${1:?usage: run_real_modes.sh <binary> <outdir> <query.bam> <target.fa> [threads]}
OUTDIR=${2:?}
QUERY=${3:?}
TARGET=${4:?}
THREADS=${5:-16}

ARGS=(--threads "$THREADS" --preset map-ont -m 2 -M 5 -o 2,24 -e 1,0)

mkdir -p "$OUTDIR"

for mode in hp-tr-v2 hp-tr-v3 hp-v2; do
    for ra in "" "--refAnchored"; do
        if [ -z "$ra" ]; then
            name="$mode"
        else
            name="${mode}-refanchored"
        fi
        "$BIN" --mode "$mode" "${ARGS[@]}" $ra \
            -q "$QUERY" \
            -t "$TARGET" \
            --out "$OUTDIR/${name}.csv"
        echo "done: $name"
    done
done

# aligned（全比对质量 metric）：不做 --refAnchored 交叉（对该 mode 无效）
"$BIN" --mode aligned "${ARGS[@]}" \
    -q "$QUERY" \
    -t "$TARGET" \
    --out "$OUTDIR/aligned.csv"
echo "done: aligned"
