#!/bin/bash

set -euo pipefail


: "${PROJECT_DIRECTORY:?Set PROJECT_DIRECTORY to your data directory}"
SRA_DIR=$PROJECT_DIRECTORY/sra
LOG_DIR=$PROJECT_DIRECTORY/logs

mkdir -p "$SRA_DIR" "$LOG_DIR"

LIST=${1:-$PROJECT_DIRECTORY/logs/missing_sra_from_full_list.txt}
if [[ ! -s "$LIST" ]]; then
    echo "ERROR: Missing or empty accession list: $LIST" >&2
    exit 1
fi

prefetch --option-file "$LIST" \
  -O "$SRA_DIR" \
  --max-size u \
  2>&1 | tee "$LOG_DIR/prefetch_missing_sra.log"

