#!/bin/bash

set -euo pipefail


: "${PROJECT_DIRECTORY:?Set PROJECT_DIRECTORY to your data directory}"
LIST=${1:-$PROJECT_DIRECTORY/completed_or_missing_runs.txt}
THREADS=${THREADS:-4}

SRA_DIR=$PROJECT_DIRECTORY/sra
FASTQ_DIR=$PROJECT_DIRECTORY/fastq
TMP_DIR=$PROJECT_DIRECTORY/tmp/fasterq_generic
LOG_DIR=$PROJECT_DIRECTORY/logs/fasterq_generic

mkdir -p "$SRA_DIR" "$FASTQ_DIR" "$TMP_DIR" "$LOG_DIR"

while read -r ACC; do
    [[ -z "$ACC" ]] && continue

    echo "========== Processing $ACC =========="

    paired_ok=false
    single_ok=false

    if [[ -f "$FASTQ_DIR/${ACC}_1.fastq.gz" && -f "$FASTQ_DIR/${ACC}_2.fastq.gz" ]]; then
        if gzip -t "$FASTQ_DIR/${ACC}_1.fastq.gz" && gzip -t "$FASTQ_DIR/${ACC}_2.fastq.gz"; then
            paired_ok=true
        else
            echo "Corrupted paired FASTQ detected; removing."
            rm -f "$FASTQ_DIR/${ACC}_1.fastq.gz" "$FASTQ_DIR/${ACC}_2.fastq.gz"
        fi
    fi

    if [[ -f "$FASTQ_DIR/${ACC}.fastq.gz" ]]; then
        if gzip -t "$FASTQ_DIR/${ACC}.fastq.gz"; then
            single_ok=true
        else
            echo "Corrupted single-end FASTQ detected; removing."
            rm -f "$FASTQ_DIR/${ACC}.fastq.gz"
        fi
    fi

    if [[ "$paired_ok" == true || "$single_ok" == true ]]; then
        echo "$ACC already complete and valid. Skipping."
        continue
    fi

    echo "Running prefetch if SRA is absent."
    if [[ ! -s "$SRA_DIR/$ACC/$ACC.sra" && ! -s "$SRA_DIR/$ACC.sra" ]]; then
        prefetch "$ACC" -O "$SRA_DIR" --max-size u \
            2>&1 | tee "$LOG_DIR/${ACC}_prefetch.log"
    fi

    if [[ -s "$SRA_DIR/$ACC/$ACC.sra" ]]; then
        SRA_INPUT="$SRA_DIR/$ACC/$ACC.sra"
    elif [[ -s "$SRA_DIR/$ACC.sra" ]]; then
        SRA_INPUT="$SRA_DIR/$ACC.sra"
    else
        echo "ERROR: Missing SRA file for $ACC after prefetch."
        continue
    fi

    mkdir -p "$TMP_DIR/$ACC"
    echo "Running fasterq-dump."
    fasterq-dump "$SRA_INPUT" \
        -O "$FASTQ_DIR" \
        -t "$TMP_DIR/$ACC" \
        -e "$THREADS" \
        -p \
        2>&1 | tee "$LOG_DIR/${ACC}_fasterq.log"

    if ls "$FASTQ_DIR/${ACC}"*.fastq >/dev/null 2>&1; then
        echo "Compressing FASTQ output."
        if command -v pigz >/dev/null 2>&1; then
            pigz -p "$THREADS" "$FASTQ_DIR/${ACC}"*.fastq
        else
            gzip "$FASTQ_DIR/${ACC}"*.fastq
        fi
    else
        echo "WARNING: No FASTQ produced for $ACC."
    fi

    rm -rf "$TMP_DIR/$ACC"
    echo "Finished $ACC"
done < "$LIST"

echo "All jobs completed."
