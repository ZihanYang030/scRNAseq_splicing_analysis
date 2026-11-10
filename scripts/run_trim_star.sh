#!/bin/bash

set -euo pipefail


TRIMMOMATIC_BIN=${TRIMMOMATIC_BIN:-trimmomatic}
THREADS=${THREADS:-4}

: "${PROJECT_DIRECTORY:?Set PROJECT_DIRECTORY to your data directory}"

: "${RUN_NAME:?ERROR: RUN_NAME is not set}"
: "${LIST:?ERROR: LIST is not set}"

FASTQ_LINK_DIR="$PROJECT_DIRECTORY/fastq_by_gsm/$RUN_NAME"
TRIMMED_DIR="$PROJECT_DIRECTORY/trimmed_by_gsm/$RUN_NAME"
STAR_INDEX_DIR="$PROJECT_DIRECTORY/star_index"
STAR_OUT_DIR="$PROJECT_DIRECTORY/star_out_by_gsm/$RUN_NAME"

mkdir -p "$TRIMMED_DIR" "$STAR_OUT_DIR" "$PROJECT_DIRECTORY/logs"

if [[ ! -f "$LIST" ]]; then
    echo "ERROR: LIST not found: $LIST"
    exit 1
fi

if [[ ! -d "$FASTQ_LINK_DIR" ]]; then
    echo "ERROR: FASTQ link directory not found: $FASTQ_LINK_DIR"
    exit 1
fi

if [[ ! -f "$STAR_INDEX_DIR/Genome" || ! -f "$STAR_INDEX_DIR/SA" || ! -f "$STAR_INDEX_DIR/SAindex" ]]; then
    echo "ERROR: STAR index incomplete: $STAR_INDEX_DIR"
    exit 1
fi

process_one_srr() {
    local SRR="$1"
    echo "Processing ${SRR}"

    if [[ -s "$STAR_OUT_DIR/${SRR}.SJ.out.tab" && -s "$STAR_OUT_DIR/${SRR}.Log.final.out" ]]; then
        echo "Skipping ${SRR}: STAR output already exists."
        return 0
    fi

    RAW_R1="$FASTQ_LINK_DIR/${SRR}_1.fastq.gz"
    RAW_R2="$FASTQ_LINK_DIR/${SRR}_2.fastq.gz"
    RAW_SINGLE="$FASTQ_LINK_DIR/${SRR}.fastq.gz"

    if [[ -f "$RAW_R1" && -f "$RAW_R2" ]]; then
        LAYOUT="paired"
    elif [[ -f "$RAW_SINGLE" ]]; then
        LAYOUT="single"
    else
        echo "ERROR: FASTQ not found for ${SRR}"
        return 1
    fi

    rm -f "$STAR_OUT_DIR/${SRR}."*

    if [[ "$LAYOUT" == "paired" ]]; then
        TRIM_R1="$TRIMMED_DIR/${SRR}_1_paired.fastq.gz"
        TRIM_R2="$TRIMMED_DIR/${SRR}_2_paired.fastq.gz"
        UNP_R1="$TRIMMED_DIR/${SRR}_1_unpaired.fastq.gz"
        UNP_R2="$TRIMMED_DIR/${SRR}_2_unpaired.fastq.gz"

        if [[ ! -s "$TRIM_R1" || ! -s "$TRIM_R2" ]]; then
            "$TRIMMOMATIC_BIN" PE -threads "$THREADS" \
                "$RAW_R1" "$RAW_R2" \
                "$TRIM_R1" "$UNP_R1" \
                "$TRIM_R2" "$UNP_R2" \
                HEADCROP:10 TRAILING:15 AVGQUAL:20 MINLEN:30
        fi

        STAR \
            --runThreadN "$THREADS" \
            --genomeDir "$STAR_INDEX_DIR" \
            --readFilesIn "$TRIM_R1" "$TRIM_R2" \
            --readFilesCommand zcat \
            --outFileNamePrefix "$STAR_OUT_DIR/${SRR}." \
            --outSAMtype None
    else
        TRIM_SINGLE="$TRIMMED_DIR/${SRR}_trimmed.fastq.gz"

        if [[ ! -s "$TRIM_SINGLE" ]]; then
            "$TRIMMOMATIC_BIN" SE -threads "$THREADS" \
                "$RAW_SINGLE" \
                "$TRIM_SINGLE" \
                HEADCROP:10 TRAILING:15 AVGQUAL:20 MINLEN:30
        fi

        STAR \
            --runThreadN "$THREADS" \
            --genomeDir "$STAR_INDEX_DIR" \
            --readFilesIn "$TRIM_SINGLE" \
            --readFilesCommand zcat \
            --outFileNamePrefix "$STAR_OUT_DIR/${SRR}." \
            --outSAMtype None
    fi

    grep -E "Number of input reads|Uniquely mapped reads number|Uniquely mapped reads %" \
        "$STAR_OUT_DIR/${SRR}.Log.final.out" || true
    wc -l "$STAR_OUT_DIR/${SRR}.SJ.out.tab" || true
}

FAILED=0
while IFS= read -r SRR || [[ -n "$SRR" ]]; do
    SRR=${SRR%$'\r'}
    [[ -z "$SRR" ]] && continue
    if ! process_one_srr "$SRR"; then
        FAILED=$((FAILED + 1))
    fi
done < "$LIST"

if [[ "$FAILED" -gt 0 ]]; then
    echo "ERROR: ${FAILED} SRR(s) failed."
    exit 1
fi

echo "Done."
date
