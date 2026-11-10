#!/bin/bash

set -euo pipefail


: "${PROJECT_DIRECTORY:?Set PROJECT_DIRECTORY to your data directory}"
LIST=${1:-$PROJECT_DIRECTORY/SRRs_L23ITCTX_Oligo_all.txt}
THREADS=${THREADS:-4}
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGION=$(python3 "$SCRIPT_DIR/regions.py" samtools bam_window)
STAR_INDEX_DIR=$PROJECT_DIRECTORY/star_index
REGION_BAM_DIR=$PROJECT_DIRECTORY/BAM/Pakap_region_L23ITCTX_Oligo_all
TMP_FULL_BAM_DIR=$PROJECT_DIRECTORY/BAM/tmp_full_bam_Pakap_region_all

mkdir -p "$REGION_BAM_DIR" "$TMP_FULL_BAM_DIR" "$PROJECT_DIRECTORY/logs"

if [[ ! -s "$LIST" ]]; then
    echo "ERROR: Missing or empty accession list: $LIST" >&2
    exit 1
fi

find_trimmed_dir() {
    local SRR="$1"
    if [[ -s "$PROJECT_DIRECTORY/trimmed_by_gsm/L23ITCTX_matched/${SRR}_1_paired.fastq.gz" ]]; then
        echo "$PROJECT_DIRECTORY/trimmed_by_gsm/L23ITCTX_matched"
        return 0
    fi
    if [[ -s "$PROJECT_DIRECTORY/trimmed_by_gsm/Oligo_matched/${SRR}_1_paired.fastq.gz" ]]; then
        echo "$PROJECT_DIRECTORY/trimmed_by_gsm/Oligo_matched"
        return 0
    fi
    return 1
}

while IFS= read -r SRR || [[ -n "$SRR" ]]; do
    SRR=${SRR%$'\r'}
    [[ -z "$SRR" ]] && continue
    REGION_BAM="$REGION_BAM_DIR/${SRR}.Pakap_region.bam"
    REGION_BAI="$REGION_BAM_DIR/${SRR}.Pakap_region.bam.bai"

    if [[ -s "$REGION_BAM" && -s "$REGION_BAI" ]]; then
        echo "Skipping $SRR: regional BAM already exists."
        continue
    fi

    TRIMMED_DIR=$(find_trimmed_dir "$SRR") || {
        echo "ERROR: Cannot find trimmed FASTQ for $SRR"
        continue
    }

    R1="$TRIMMED_DIR/${SRR}_1_paired.fastq.gz"
    R2="$TRIMMED_DIR/${SRR}_2_paired.fastq.gz"

    SRR_TMP="$TMP_FULL_BAM_DIR/$SRR"
    rm -rf "$SRR_TMP"
    mkdir -p "$SRR_TMP"
    PREFIX="$SRR_TMP/${SRR}."

    STAR \
        --runThreadN "$THREADS" \
        --genomeDir "$STAR_INDEX_DIR" \
        --readFilesIn "$R1" "$R2" \
        --readFilesCommand zcat \
        --outFileNamePrefix "$PREFIX" \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI AS nM

    FULL_BAM="${PREFIX}Aligned.sortedByCoord.out.bam"
    if [[ ! -s "$FULL_BAM" ]]; then
        echo "ERROR: STAR did not produce BAM for $SRR"
        rm -rf "$SRR_TMP"
        continue
    fi

    samtools index "$FULL_BAM"
    samtools view -b "$FULL_BAM" "$REGION" > "$REGION_BAM"
    samtools index "$REGION_BAM"
    rm -rf "$SRR_TMP"

    echo "Saved regional BAM: $REGION_BAM"
done < "$LIST"

echo "Done."
date
