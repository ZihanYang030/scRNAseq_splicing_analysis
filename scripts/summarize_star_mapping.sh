#!/usr/bin/env bash
set -euo pipefail
: "${PROJECT_DIRECTORY:?Set PROJECT_DIRECTORY}"
RUN_NAME=${1:?Usage: summarize_star_mapping.sh RUN_NAME}
STAR_OUT_DIR=$PROJECT_DIRECTORY/star_out_by_gsm/$RUN_NAME
SUMMARY=$PROJECT_DIRECTORY/results/${RUN_NAME}_STAR_mapping_summary.tsv
mkdir -p "$PROJECT_DIRECTORY/results"
printf 'SRR\tinput_reads\tunique_mapped_reads\tunique_mapped_percent\tSJ_count\n' > "$SUMMARY"
shopt -s nullglob
for LOG in "$STAR_OUT_DIR"/*.Log.final.out; do
    SRR=$(basename "$LOG" .Log.final.out)
    SJ="$STAR_OUT_DIR/${SRR}.SJ.out.tab"
    input_reads=$(awk -F'|' '/Number of input reads/ {gsub(/^[ \t]+|[ \t]+$/, "", $2); print $2}' "$LOG")
    unique_reads=$(awk -F'|' '/Uniquely mapped reads number/ {gsub(/^[ \t]+|[ \t]+$/, "", $2); print $2}' "$LOG")
    unique_pct=$(awk -F'|' '/Uniquely mapped reads %/ {gsub(/^[ \t]+|[ \t]+$/, "", $2); print $2}' "$LOG")
    sj_count=$(wc -l < "$SJ")
    printf '%s\t%s\t%s\t%s\t%s\n' "$SRR" "$input_reads" "$unique_reads" "$unique_pct" "$sj_count" >> "$SUMMARY"
done
head "$SUMMARY"
