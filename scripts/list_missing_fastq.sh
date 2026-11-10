#!/usr/bin/env bash
set -euo pipefail
: "${PROJECT_DIRECTORY:?Set PROJECT_DIRECTORY}"
LIST=${1:?Usage: list_missing_fastq.sh ACCESSION_LIST OUTPUT_LIST}
OUT=${2:?Usage: list_missing_fastq.sh ACCESSION_LIST OUTPUT_LIST}
FASTQ_SHARED=$PROJECT_DIRECTORY/fastq
: > "$OUT"
while read -r ACC; do
    [[ -z "$ACC" ]] && continue
    if [[ -f "$FASTQ_SHARED/${ACC}_1.fastq.gz" && -f "$FASTQ_SHARED/${ACC}_2.fastq.gz" ]]; then
        continue
    elif [[ -f "$FASTQ_SHARED/${ACC}.fastq.gz" ]]; then
        continue
    fi
    echo "$ACC" >> "$OUT"
done < "$LIST"
wc -l "$OUT"
