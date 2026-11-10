#!/usr/bin/env python3

import csv
import os
from pathlib import Path
from regions import load_regions

base_dir = Path(os.environ["PROJECT_DIRECTORY"])
regions = load_regions()

runs = [
    {"run_name": "L23ITCTX_matched", "srr_list": base_dir / "L23ITCTX_Acc_List.txt", "class_label": "Glutamatergic", "subclass_label": "L2/3 IT CTX"},
    {"run_name": "Oligo_matched", "srr_list": base_dir / "Oligo_Acc_List.txt", "class_label": "Non-neuronal", "subclass_label": "Oligo"},
]

chrom = regions["gene_span"]["chrom"]
#gene_start = regions["start"]
gene_start = regions["gene_span"]["start"]
gene_end = regions["gene_span"]["end"]
buffer = 1000

fieldnames = [
    "SRR", "sample_name", "class_label", "subclass_label", "chrom",
    "intron_start", "intron_end", "strand_code", "intron_motif",
    "annotated", "unique_reads", "multi_reads", "max_overhang",
]

for cfg in runs:
    run_name = cfg["run_name"]
    list_path = cfg["srr_list"]
    #star_out_dir = / base_dir / "star_out_by_gsm" / run_name
    star_out_dir = base_dir / "star_out_by_gsm" / run_name
    results_dir = base_dir / "results" / run_name
    results_dir.mkdir(parents=True, exist_ok=True)

    out_path = results_dir / "all_Pakap_manual_region_junctions.csv"
    missing_path = results_dir / "missing_SJ_out_tab_current.txt"

    with open(list_path, "r") as f:
        srr_list = [line.strip().split()[0] for line in f if line.strip()]

    all_rows = []
    missing = []

    for srr in srr_list:
        sj_path = star_out_dir / f"{srr}.SJ.out.tab"
        if not sj_path.exists():
            missing.append(srr)
            continue

        with open(sj_path, "r") as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) != 9:
                    continue

                sj_chrom = parts[0]
                intron_start = int(parts[1])
                intron_end = int(parts[2])

                if sj_chrom != chrom:
                    continue
                if intron_start < gene_start - buffer:
                    continue
                if intron_end > gene_end + buffer:
                    continue

                all_rows.append({
                    "SRR": srr,
                    "sample_name": srr,
                    "class_label": cfg["class_label"],
                    "subclass_label": cfg["subclass_label"],
                    "chrom": sj_chrom,
                    "intron_start": intron_start,
                    "intron_end": intron_end,
                    "strand_code": parts[3],
                    "intron_motif": parts[4],
                    "annotated": parts[5],
                    #"annotated": parts[6],
                    "unique_reads": parts[6],
                    "multi_reads": parts[7],
                    "max_overhang": parts[8],
                })

    with open(out_path, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_rows)

    if missing:
        with open(missing_path, "w") as out:
            for srr in missing:
                out.write(srr + "\n")

    print(f"Run: {run_name}")
    print(f"saved junctions: {out_path}")
    print(f"Total junction rows: {len(all_rows)}")
    print(f"Missing SJ.out.tab files: {len(missing)}")
