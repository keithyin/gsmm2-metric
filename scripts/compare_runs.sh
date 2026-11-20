#!/usr/bin/env bash
# 对比两个输出目录中同名 csv 的内容一致性（排序后逐字节，容忍 HashMap 行序随机性）。
# 用法: scripts/compare_runs.sh <dirA> <dirB>
set -uo pipefail

DIRA=${1:?usage: compare_runs.sh <dirA> <dirB>}
DIRB=${2:?usage: compare_runs.sh <dirA> <dirB>}

fail=0
for f in "$DIRA"/*.csv; do
    name=$(basename "$f")
    other="$DIRB/$name"
    if [ ! -f "$other" ]; then
        echo "FAIL  $name  (missing in $DIRB)"
        fail=1
        continue
    fi
    if diff <(tail -n +2 "$f" | sort) <(tail -n +2 "$other" | sort) > /dev/null && \
       diff <(head -1 "$f") <(head -1 "$other") > /dev/null; then
        n_lines=$(tail -n +2 "$f" | wc -l)
        echo "PASS  $name  ($n_lines data lines)"
    else
        n_diff=$(diff <(tail -n +2 "$f" | sort) <(tail -n +2 "$other" | sort) | grep -c '^[<>]')
        echo "FAIL  $name  ($n_diff data lines differ)"
        diff <(tail -n +2 "$f" | sort) <(tail -n +2 "$other" | sort) | head -30
        fail=1
    fi
done

exit $fail
