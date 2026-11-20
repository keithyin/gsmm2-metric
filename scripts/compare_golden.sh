#!/usr/bin/env bash
# 端到端一致性验收：用当前构建的二进制重跑全部 mode，与 test_data/golden 逐字节对比。
# 用法: scripts/compare_golden.sh [binary]
set -uo pipefail

BIN=${1:-./target/release/gsmm2-metric}
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUTDIR=$(mktemp -d)
trap 'rm -rf "$OUTDIR"' EXIT

bash "$ROOT/scripts/run_all_modes.sh" "$BIN" "$OUTDIR" > /dev/null

fail=0
for f in "$ROOT"/test_data/golden/*.csv; do
    name=$(basename "$f")
    # metric 输出行序来自 HashMap 迭代（每次进程运行顺序都不同，属既有行为），
    # 因此 header 行原样比较，其余行排序后逐字节比较内容集合。
    if diff <(head -1 "$f") <(head -1 "$OUTDIR/$name") > /dev/null && \
       diff <(tail -n +2 "$f" | sort) <(tail -n +2 "$OUTDIR/$name" | sort) > /dev/null; then
        echo "PASS  $name"
    else
        n_diff=$(diff <(tail -n +2 "$f" | sort) <(tail -n +2 "$OUTDIR/$name" | sort) | grep -c '^[<>]')
        echo "FAIL  $name  ($n_diff data lines differ)"
        diff <(tail -n +2 "$f" | sort) <(tail -n +2 "$OUTDIR/$name" | sort) | head -20
        fail=1
    fi
done

exit $fail
