#!/usr/bin/env bash
set -euo pipefail
: "${PROJECT_DIRECTORY:?Set PROJECT_DIRECTORY}"
RUN_NAME=${1:?Usage: link_fastq_by_run.sh RUN_NAME ACCESSION_LIST}
LIST=${2:?Usage: link_fastq_by_run.sh RUN_NAME ACCESSION_LIST}
FASTQ_SHARED=$PROJECT_DIRECTORY/fastq
LINK_DIR=$PROJECT_DIRECTORY/fastq_by_gsm/$RUN_NAME
mkdir -p "$LINK_DIR"
while read -r ACC; do
    [[ -z "$ACC" ]] && continue
    if [[ -f "$FASTQ_SHARED/${ACC}_1.fastq.gz" && -f "$FASTQ_SHARED/${ACC}_2.fastq.gz" ]]; then
        ln -sfn "$FASTQ_SHARED/${ACC}_1.fastq.gz" "$LINK_DIR/${ACC}_1.fastq.gz"
        ln -sfn "$FASTQ_SHARED/${ACC}_2.fastq.gz" "$LINK_DIR/${ACC}_2.fastq.gz"
    elif [[ -f "$FASTQ_SHARED/${ACC}.fastq.gz" ]]; then
        ln -sfn "$FASTQ_SHARED/${ACC}.fastq.gz" "$LINK_DIR/${ACC}.fastq.gz"
    else
        echo "Missing FASTQ: $ACC" >&2
    fi
done < "$LIST"
find "$LINK_DIR" -maxdepth 1 -type l | wc -l
