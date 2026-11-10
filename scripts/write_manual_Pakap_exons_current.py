#!/usr/bin/env python3

import csv
import os
from pathlib import Path
from regions import load_regions

base_dir = Path(os.environ["PROJECT_DIRECTORY"])
regions = load_regions()
run_names = ["L23ITCTX_matched", "Oligo_matched"] # earlier defined in calc_Pakap_manual_EIP_current, matched with desired cell tyeps May 2, 2026

gene_id = "Pakap"
gene_name = "Pakap"
chrom = regions["exon_1"]["chrom"]
strand = regions["exon_1"]["strand"]


def exon(number, region_name, exon_id=None):
    row = regions[region_name]
    return {"exon_number": number, "exon_id": exon_id or f"Pakap_manual_E{number}",
            "start": row["start"], "end": row["end"]}


E1 = exon("1", "exon_1")
E2 = exon("2", "exon_2")
E3 = exon("3", "exon_3")
E4 = exon("4", "exon_4")
E5 = exon("5", "exon_5")
E6A_INTERNAL = exon("6a", "exon_6a", "Pakap_manual_E6a_internal")
E6_TERMINAL_FULL = exon("6_terminal", "exon_6_terminal", "Pakap_manual_E6_terminal")
E7 = exon("7", "exon_7")
E8 = exon("8", "exon_8")
E9 = exon("9", "exon_9")
E10 = exon("10", "exon_10")
E11 = exon("11", "exon_11")

structures = {
    "skip_exon6": [E1, E2, E3, E4, E5, E7, E8, E9, E10, E11],
    "include_exon6a_internal_to_E7": [E1, E2, E3, E4, E5, E6A_INTERNAL, E7, E8, E9, E10, E11],
    "terminal_exon6a_6b_3UTR": [E1, E2, E3, E4, E5, E6_TERMINAL_FULL],
}

fieldnames = [
    "chrom", "source", "feature", "start", "end", "strand",
    "gene_id", "gene_name", "transcript_id", "exon_number",
    "exon_id", "manual_structure", "attributes"
]


def make_rows(structure_name, exon_defs):
    rows = []
    for exon in exon_defs:
        rows.append({
            "chrom": chrom, "source": "manual", "feature": "exon",
            "start": exon["start"], "end": exon["end"],
            "strand": strand, "gene_id": gene_id, "gene_name": gene_name,
            "transcript_id": f"manual_{structure_name}",
            "exon_number": exon["exon_number"], "exon_id": exon["exon_id"],
            "manual_structure": structure_name,
            "attributes": f'transcript_id "manual_{structure_name}"; gene_name "{gene_name}"; exon_number "{exon["exon_number"]}"; exon_id "{exon["exon_id"]}";',
        })
    return rows


for run_name in run_names:
    #manual_dir = /base_dir / "results" / run_name / "manual_Pakap_current_exons"
    manual_dir = base_dir / "results" / run_name / "manual_Pakap_current_exons"
    manual_dir.mkdir(parents=True, exist_ok=True)
    all_rows = []

    for structure_name, exon_defs in structures.items():
        rows = make_rows(structure_name, exon_defs)
        all_rows.extend(rows)
        out_file = manual_dir / f"Pakap_manual_exons_{structure_name}.csv"
        with open(out_file, "w", newline="") as out:
            writer = csv.DictWriter(out, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)
        print(f"Saved: {out_file}")

    combined_file = manual_dir / "Pakap_manual_exons_current_combined.csv"
    with open(combined_file, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_rows)
    print(f"Saved combined manual exon file: {combined_file}")
