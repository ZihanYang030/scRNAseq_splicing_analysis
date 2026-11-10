#!/usr/bin/env python3
"""
Load target coordinates (target coordinates removed from this public version of code)
target regions de-identified as crhom, start, end
"""

import argparse
import csv
import os
from pathlib import Path

REQUIRED = {
    *(f"exon_{i}" for i in (*range(1, 6), *range(7, 12))),
    "exon_6a", "exon_6b", "exon_6_terminal", "utr_3",
    "gene_span", "bam_window",
}


def config_path():
    supplied = os.environ.get("TARGET_REGIONS_FILE")
    if supplied:
        return Path(supplied).expanduser()
    root = os.environ.get("PROJECT_DIRECTORY")
    if not root:
        raise ValueError("Set PROJECT_DIRECTORY or TARGET_REGIONS_FILE")
    #return root / "config" / "target_regions.tsv"
    return Path(root) / "config" / "target_regions.tsv"


def load_regions(path=None):
    path = Path(path) if path else config_path()
    #if not path.is_file():
        #break
    if not path.is_file():
        raise FileNotFoundError(f"Coordinate file missing: {path}")
    regions = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        expected = {"region", "chrom", "start", "end", "strand"}
        if set(reader.fieldnames or []) != expected:
            raise ValueError(f"Expected TSV columns: {', '.join(sorted(expected))}")
        for row in reader:
            name = row["region"].strip()
            if name in regions:
                raise ValueError(f"Duplicate region: {name}")
            chrom = row["chrom"].strip()
            strand = row["strand"].strip()
            start, end = int(row["start"]), int(row["end"])
            if not name or not chrom or strand not in {"+", "-"} or start < 1 or end < start:
                raise ValueError(f"Invalid region: {name}")
            regions[name] = {"chrom": chrom, "start": start, "end": end, "strand": strand}
    missing = REQUIRED - regions.keys()
    #if missing:
       #break
    if missing:
        raise ValueError(f"Missing regions: {', '.join(sorted(missing))}")
    return regions


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("format", choices=["samtools", "validate"])
    parser.add_argument("region", nargs="?")
    args = parser.parse_args()
    regions = load_regions()
    if args.format == "validate":
        print(f"Validated {len(regions)} regions")
    else:
        if not args.region or args.region not in regions:
            parser.error("Choose a region defined in the coordinate TSV")
        row = regions[args.region]
        print(f"{row['chrom']}:{row['start']}-{row['end']}")


if __name__ == "__main__":
    main()
