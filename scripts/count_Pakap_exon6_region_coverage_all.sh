#!/bin/bash
set -euo pipefail


: "${PROJECT_DIRECTORY:?Set PROJECT_DIRECTORY to your data directory}"
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
BAM_DIR=$PROJECT_DIRECTORY/BAM/Pakap_region_L23ITCTX_Oligo_all
E6A_REGION=$(python3 "$SCRIPT_DIR/regions.py" samtools exon_6a)
E6B_REGION=$(python3 "$SCRIPT_DIR/regions.py" samtools exon_6b)
UTR3_REGION=$(python3 "$SCRIPT_DIR/regions.py" samtools utr_3)
E7_REGION=$(python3 "$SCRIPT_DIR/regions.py" samtools exon_7)
OUT=$PROJECT_DIRECTORY/results/Pakap_exon6a_6b_3UTR_coverage_by_cell_all.tsv

mkdir -p "$PROJECT_DIRECTORY/results"
echo -e "SRR\treads_E6a\treads_E6b\treads_3UTR\treads_E7" > "$OUT"

shopt -s nullglob
bams=("$BAM_DIR"/*.Pakap_region.bam)
if [[ ${#bams[@]} -eq 0 ]]; then
    echo "ERROR: No regional BAM files found in $BAM_DIR" >&2
    exit 1
fi
for BAM in "${bams[@]}"; do
    SRR=$(basename "$BAM" .Pakap_region.bam)

    E6A=$(samtools view -c "$BAM" "$E6A_REGION")
    E6B=$(samtools view -c "$BAM" "$E6B_REGION")
    UTR3=$(samtools view -c "$BAM" "$UTR3_REGION")
    E7=$(samtools view -c "$BAM" "$E7_REGION")

    echo -e "${SRR}\t${E6A}\t${E6B}\t${UTR3}\t${E7}" >> "$OUT"
done

echo "Saved: $OUT"
