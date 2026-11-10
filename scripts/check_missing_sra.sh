#!/bin/bash
set -euo pipefail

: "${PROJECT_DIRECTORY:?Set PROJECT_DIRECTORY to your data directory}"
ACC_LIST=${1:-$PROJECT_DIRECTORY/SRR_Acc_List.txt}
SRA_DIR=$PROJECT_DIRECTORY/sra
OUT_DIR=$PROJECT_DIRECTORY/logs

MISSING_LIST=$OUT_DIR/missing_sra_from_full_list.txt
PRESENT_LIST=$OUT_DIR/present_sra_from_full_list.txt

mkdir -p "$OUT_DIR"
: > "$MISSING_LIST"
: > "$PRESENT_LIST"

echo "Checking SRA files..."
echo "Accession list: $ACC_LIST"
echo "SRA directory:   $SRA_DIR"

while read -r SRR; do
    [[ -z "$SRR" ]] && continue

    SRA_FILE_1="$SRA_DIR/$SRR/$SRR.sra"
    SRA_FILE_2="$SRA_DIR/$SRR.sra"

    if [[ -s "$SRA_FILE_1" || -s "$SRA_FILE_2" ]]; then
        echo "$SRR" >> "$PRESENT_LIST"
    else
        echo "$SRR" >> "$MISSING_LIST"
    fi
done < "$ACC_LIST"

echo "Done."
echo "Present SRA count:"
wc -l "$PRESENT_LIST"
echo "Missing SRA count:"
wc -l "$MISSING_LIST"
echo "First 10 missing SRRs:"
head "$MISSING_LIST"
