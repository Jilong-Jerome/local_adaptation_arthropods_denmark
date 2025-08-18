#!/usr/bin/env python3
"""
SweepFinder peak detection using site-based rolling window,
constrained by a maximum genomic span. Outputs include number
of sites per peak.
"""

import argparse
import pandas as pd
import numpy as np

# ------------------- User defaults ---------------------
HIGH_THR   = 50       # LR threshold to trigger a peak
LOW_THR    = 5        # min rolling LR mean to extend a peak
ROLL_N     = 5        # number of sites in rolling window
STEP_N     = 1        # step size in number of sites
MAX_SPAN   = 10000    # max span in bp allowed within window
# -------------------------------------------------------

def extract_constrained_peaks(block, high_thr, low_thr, roll_n, step_n, max_span):
    pos   = block["location"].to_numpy()
    lr    = block["LR"].to_numpy()
    alpha = block["alpha"].to_numpy()

    assigned = np.full(len(block), False)
    seeds = np.argsort(-lr)  # descending LR

    for idx in seeds:
        if assigned[idx] or lr[idx] < high_thr:
            continue

        peak_sites = set([idx])

        # ---- Extend to the right ----
        i = idx + step_n
        while i + roll_n <= len(block):
            window = np.arange(i, i + roll_n)
            if assigned[window].any():
                break
            span = pos[window[-1]] - pos[window[0]]
            if span > max_span:
                break
            if lr[window].mean() < low_thr:
                break
            peak_sites.update(window)
            i += step_n

        # ---- Extend to the left ----
        i = idx - step_n
        while i - roll_n + 1 >= 0:
            window = np.arange(i - roll_n + 1, i + 1)
            if assigned[window].any():
                break
            span = pos[window[-1]] - pos[window[0]]
            if span > max_span:
                break
            if lr[window].mean() < low_thr:
                break
            peak_sites.update(window)
            i -= step_n

        # assign all peak sites
        for i in peak_sites:
            assigned[i] = True

        region = block.iloc[list(peak_sites)]
        yield {
            "species":       region["sp"].iloc[0],
            "population":    region["pop"].iloc[0],
            "chromosome":    region["chrom"].iloc[0],
            "start":         int(region["location"].min()),
            "end":           int(region["location"].max()),
            "highest_LR":    region["LR"].max(),
            "lowest_alpha":  region["alpha"].min(),
            "region_size":   int(region["location"].max() - region["location"].min()),
            "num_sites":     len(region)
        }

def build_summary(df, high_thr, low_thr, roll_n, step_n, max_span):
    summary = []
    for _, block in df.groupby(["sp", "pop", "chrom"]):
        block = block.sort_values("location").reset_index(drop=True)
        summary.extend(extract_constrained_peaks(block, high_thr, low_thr, roll_n, step_n, max_span))
    return pd.DataFrame(summary, columns=[
        "species", "population", "chromosome",
        "start", "end", "highest_LR", "lowest_alpha",
        "region_size", "num_sites"
    ])

def main():
    parser = argparse.ArgumentParser(description="SweepFinder peak summary with site- and span-based rolling window")
    parser.add_argument("input", help="SweepFinder result file (TSV or space-delimited)")
    parser.add_argument("-o", "--out", default="selection_peaks.tsv", help="Output TSV file")
    parser.add_argument("--high", type=float, default=HIGH_THR, help=f"Seed LR threshold (default: {HIGH_THR})")
    parser.add_argument("--low", type=float, default=LOW_THR, help=f"Rolling mean LR threshold (default: {LOW_THR})")
    parser.add_argument("--roll", type=int, default=ROLL_N, help=f"Rolling window size in sites (default: {ROLL_N})")
    parser.add_argument("--step", type=int, default=STEP_N, help=f"Step size in sites (default: {STEP_N})")
    parser.add_argument("--span", type=int, default=MAX_SPAN, help=f"Max genomic span per window (bp) (default: {MAX_SPAN})")
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep=r"\s+")
    summary = build_summary(df, args.high, args.low, args.roll, args.step, args.span)
    summary.to_csv(args.out, sep="\t", index=False)
    print(f"{len(summary)} peak regions written to {args.out}")

if __name__ == "__main__":
    main()

