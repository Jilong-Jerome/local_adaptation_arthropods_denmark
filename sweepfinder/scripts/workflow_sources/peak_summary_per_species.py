#!/usr/bin/env python3
"""
merge_selection_peaks.py  –  v2  (individual-peak rows)

Workflow
========
1. Merge nearby peaks *within* each (species, population, chromosome) if gap ≤ --distance.
2. Merge *overlapping* peaks *across* populations on the same chromosome to form
   broader selection regions.
3. Emit one row **per original peak** that falls inside a merged region.

Input
-----
TSV with at least these columns:
    species, population, chromosome, start, end, highest_LR, lowest_alpha

Output
------
TSV, one line per (merged_region × original_peak), columns:
    region_id          unique ID for merged region
    species, chromosome, start, end   merged-region boundaries
    num_populations    how many populations have peaks in this region
    population_id
    peak_start, peak_end              coordinates of the contributing peak
    overlap_fraction                  (peak_len / region_len, inclusive)
    highest_LR, lowest_alpha          values from that peak

Usage
-----
python merge_selection_peaks.py peaks.tsv merged_peaks.tsv \
        --distance 5000
"""

import argparse
from collections import defaultdict
from typing import List, Tuple

import pandas as pd


# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("infile",  help="Input TSV with peak data")
    p.add_argument("outfile", help="Output TSV after merging")
    p.add_argument(
        "-d", "--distance", type=int, default=5000,
        help="Max gap (bp) to merge nearby peaks within a population [default: 5000]"
    )
    return p.parse_args()


# ---------------------------------------------------------------------
# Step 1 – merge within population
# ---------------------------------------------------------------------
def merge_within_population(df: pd.DataFrame, distance: int) -> pd.DataFrame:
    """Collapse nearby peaks per species/population/chromosome."""
    merged_records = []
    group_cols = ["species", "population", "chromosome"]

    for (sp, pop, chrom), grp in df.groupby(group_cols):
        grp = grp.sort_values("start")
        cur_start, cur_end = None, None
        cur_lrs, cur_alphas = [], []

        for row in grp.itertuples(index=False):
            if cur_start is None:
                cur_start, cur_end = row.start, row.end
                cur_lrs     = [row.highest_LR]
                cur_alphas  = [row.lowest_alpha]
                continue

            # close enough → merge
            if row.start <= cur_end + distance:
                cur_end = max(cur_end, row.end)
                cur_lrs.append(row.highest_LR)
                cur_alphas.append(row.lowest_alpha)
            else:
                merged_records.append({
                    "species": sp, "population": pop, "chromosome": chrom,
                    "start": cur_start, "end": cur_end,
                    "highest_LR": max(cur_lrs),
                    "lowest_alpha": min(cur_alphas)
                })
                # start new interval
                cur_start, cur_end = row.start, row.end
                cur_lrs     = [row.highest_LR]
                cur_alphas  = [row.lowest_alpha]

        # last interval
        merged_records.append({
            "species": sp, "population": pop, "chromosome": chrom,
            "start": cur_start, "end": cur_end,
            "highest_LR": max(cur_lrs),
            "lowest_alpha": min(cur_alphas)
        })

    return pd.DataFrame(merged_records)


# ---------------------------------------------------------------------
# Step 2 – merge across populations
# ---------------------------------------------------------------------
def merge_across_populations(df: pd.DataFrame) -> pd.DataFrame:
    """
    Build merged regions across populations and return one row
    **per original contributing peak**.
    """
    rows = []
    region_counter = defaultdict(int)

    for (sp, chrom), grp in df.groupby(["species", "chromosome"]):
        # sort all within-population intervals
        intervals = grp.sort_values("start").to_dict("records")
        clusters = []   # list of dicts representing merged regions

        # greedy clustering of overlapping intervals
        for iv in intervals:
            added = None
            for cl in clusters:
                if iv["start"] <= cl["end"] and iv["end"] >= cl["start"]:
                    added = cl
                    break
            if added:
                added["start"] = min(added["start"], iv["start"])
                added["end"]   = max(added["end"],   iv["end"])
                added["records"].append(iv)
                added["pops"].add(iv["population"])
            else:
                clusters.append({
                    "start": iv["start"],
                    "end":   iv["end"],
                    "records": [iv],               # original peaks
                    "pops": {iv["population"]}    # set of populations
                })

        # emit one output row per original peak
        for cl in clusters:
            region_counter[(sp, chrom)] += 1
            region_id = f"{sp}_{chrom}_R{region_counter[(sp, chrom)]}"
            region_len = cl["end"] - cl["start"] + 1  # inclusive

            for rec in cl["records"]:
                peak_len = rec["end"] - rec["start"] + 1
                rows.append({
                    "region_id": region_id,
                    "species": sp,
                    "chromosome": chrom,
                    "start": cl["start"],
                    "end": cl["end"],
                    "num_populations": len(cl["pops"]),
                    "population_id": rec["population"],
                    "peak_start": rec["start"],
                    "peak_end":   rec["end"],
                    "overlap_fraction": round(peak_len / region_len, 5),
                    "highest_LR": round(rec["highest_LR"], 5),
                    "lowest_alpha": f"{rec['lowest_alpha']:.8g}"
                })

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------
# main
# ---------------------------------------------------------------------
def main() -> None:
    args = parse_args()
    df = pd.read_csv(args.infile, sep="\t")

    required = {"species", "population", "chromosome",
                "start", "end", "highest_LR", "lowest_alpha"}
    if not required.issubset(df.columns):
        raise ValueError(f"Input is missing required columns: {required}")

    stage1 = merge_within_population(df, args.distance)
    stage2 = merge_across_populations(stage1)

    stage2.to_csv(args.outfile, sep="\t", index=False)
    print(f"✓ Wrote {len(stage2)} rows to {args.outfile}")


if __name__ == "__main__":
    main()

