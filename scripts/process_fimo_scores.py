#!/usr/bin/env python

import pandas as pd
import numpy as np
import sys

def process_fimo(fimo_path, output_path):
    # Load FIMO results
    df = pd.read_csv(fimo_path, sep="\t", comment="#")

    # Sanity check
    required_cols = {"sequence_name", "start", "stop", "p-value"}
    if not required_cols.issubset(df.columns):
        raise ValueError(f"FIMO file missing required columns: {required_cols - set(df.columns)}")

    # Construct region string: chr:start-end
    df["region"] = df["sequence_name"].astype(str) + ":" + df["start"].astype(str) + "-" + df["stop"].astype(str)

    # Keep only best match per region (lowest p-value)
    best_hits = df.sort_values("p-value").drop_duplicates("region")

    # Compute -log10(p-value)
    best_hits["motif_score"] = -np.log10(best_hits["p-value"])

    # Save
    best_hits = best_hits[["region", "motif_score"]]
    best_hits.to_csv(output_path, sep="\t", index=False)
    print(f"✅ Saved {len(best_hits)} motif scores to {output_path}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python process_fimo_scores.py <fimo.tsv> <output.tsv>")
        sys.exit(1)

    fimo_input = sys.argv[1]
    output_file = sys.argv[2]
    process_fimo(fimo_input, output_file)