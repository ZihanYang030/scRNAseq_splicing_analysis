#!/usr/bin/env bash
set -euo pipefail
: "${PROJECT_DIRECTORY:?Set PROJECT_DIRECTORY}"
OUT=$PROJECT_DIRECTORY/SRRs_L23ITCTX_Oligo_all.txt
awk 'NF && !seen[$1]++ {print $1}' \
  "$PROJECT_DIRECTORY/L23ITCTX_Acc_List.txt" \
  "$PROJECT_DIRECTORY/Oligo_Acc_List.txt" > "$OUT"
wc -l "$OUT"
