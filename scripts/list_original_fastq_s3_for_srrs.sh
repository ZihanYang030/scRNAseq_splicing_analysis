#!/bin/bash

set -euo pipefail

: "${PROJECT_DIRECTORY:?Set PROJECT_DIRECTORY to your data directory}"
cd "$PROJECT_DIRECTORY"
mkdir -p "$PROJECT_DIRECTORY/metadata"

SRR_LIST=${1:-$PROJECT_DIRECTORY/SRR_Acc_List.txt}
OUT_TSV="$PROJECT_DIRECTORY/metadata/PRJNA772116_original_fastq_s3_names_from_aws.tsv"
FAILED="$PROJECT_DIRECTORY/metadata/PRJNA772116_original_fastq_s3_failed_srr.txt"

: > "$OUT_TSV"
: > "$FAILED"
echo -e "SRR\ts3_path\tfile_name\tread_pair\tbucket" >> "$OUT_TSV"

if ! command -v aws >/dev/null 2>&1; then
    echo "ERROR: aws command not found. Install AWS CLI and add it to PATH."
    exit 1
fi

infer_read_pair() {
    local fname="$1"
    if [[ "$fname" =~ (^|[_\.-])R1([_\.-]|$) ]]; then
        echo "R1"
    elif [[ "$fname" =~ (^|[_\.-])R2([_\.-]|$) ]]; then
        echo "R2"
    else
        echo ""
    fi
}

BUCKETS=$(seq 1 30)
i=0
total=$(grep -c '^SRR' "$SRR_LIST" || true)

while read -r SRR; do
    [[ -z "$SRR" ]] && continue
    [[ ! "$SRR" =~ ^SRR ]] && continue

    i=$((i + 1))
    echo "[$i/$total] Checking $SRR"
    found=0

    for b in $BUCKETS; do
        bucket="sra-pub-src-$b"
        prefix="s3://${bucket}/${SRR}/"
        listing=$(aws s3 ls "$prefix" --no-sign-request 2>/dev/null || true)

        if [[ -n "$listing" ]]; then
            found=1
            echo "$listing" | while read -r date time size fname; do
                [[ -z "${fname:-}" ]] && continue
                if [[ "$fname" == *.fastq* || "$fname" == *.fq* ]]; then
                    s3_path="${prefix}${fname}"
                    read_pair=$(infer_read_pair "$fname")
                    echo -e "${SRR}\t${s3_path}\t${fname}\t${read_pair}\t${bucket}" >> "$OUT_TSV"
                fi
            done
            break
        fi
    done

    if [[ "$found" -eq 0 ]]; then
        echo "$SRR" >> "$FAILED"
    fi
done < "$SRR_LIST"

ls -lh "$OUT_TSV"
wc -l "$FAILED"
