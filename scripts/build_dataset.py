#!/usr/bin/env python

import pandas as pd
import pyranges as pr
import numpy as np
from sklearn.model_selection import train_test_split

# === Load FIMO output from log score files ===
def load_fimo_scores_from_log(path):
    df = pd.read_csv(path, sep="\t", header=0)  # region, motif_score

    # Split region string: chr:start-end
    split = df["region"].str.extract(r"^(?P<Chromosome>[^:]+):(?P<Start>\d+)-(?P<End>\d+)$")
    df = pd.concat([split, df["motif_score"]], axis=1)

    # Cast to correct types
    df["Start"] = df["Start"].astype(int)
    df["End"] = df["End"].astype(int)

    return df[["Chromosome", "Start", "End", "motif_score"]]

# === Load 500bp regions ===
def load_chip_regions(bed_path, signal_tab):
    bed = pd.read_csv(bed_path, sep="\t", header=None,
                      names=["Chromosome", "Start", "End", "region"])
    signal = pd.read_csv(signal_tab, sep="\t", header=None,
                         names=["region", "size", "covered", "sum", "mean", "mean0"])
    return bed, signal

# === Build feature table ===
def build_feature_table(fimo_log_path, bed_path, signal_tab):
    fimo_df = load_fimo_scores_from_log(fimo_log_path)
    bed_df, signal_df = load_chip_regions(bed_path, signal_tab)

    # Overlap with PyRanges
    fimo_pr = pr.PyRanges(fimo_df)
    bed_pr = pr.PyRanges(bed_df)
    joined = bed_pr.join(fimo_pr)

    # Keep best match per region
    hits = joined.df.sort_values("motif_score", ascending=False).drop_duplicates("region")
    hits = hits[["region", "motif_score"]]

    # Merge with ChIP signal
    merged = pd.merge(hits, signal_df[["region", "mean0"]], on="region")
    return merged

# === Build dataset ===
pos = build_feature_table("fimo_outputs/features/fimo_pos_log_scores.tsv",
                          "data/ctcf_positives_named.bed",
                          "features/pos_chip_signal.tab")
neg = build_feature_table("fimo_outputs/features/fimo_neg_log_scores.tsv",
                          "data/ctcf_negatives_named.bed",
                          "features/neg_chip_signal.tab")

# Add labels
pos["label"] = 1
neg["label"] = 0

# Combine and save
df = pd.concat([pos, neg], ignore_index=True)
df = df.rename(columns={"mean0": "chip_signal"})
df.to_csv("features/final_dataset.csv", index=False)
print(f"Saved features/final_dataset.csv with {len(df)} examples.")

# === Train/test split ===
train_df, test_df = train_test_split(
    df,
    test_size=0.2,
    stratify=df["label"],
    random_state=42
)

# Save train/test
train_df.to_csv("features/train_dataset.csv", index=False)
test_df.to_csv("features/test_dataset.csv", index=False)

print(f"🧪 Saved {len(train_df)} training and {len(test_df)} testing examples.")