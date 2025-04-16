# Script: process_fimo_scores.py
# Purpose: Extract -log10(p-value) from FIMO tsv and return best match per region

import pandas as pd
import numpy as np
import argparse

def process_fimo(fimo_path, output_path):
    df = pd.read_csv(fimo_path, sep='\t', comment='#')
    df['region'] = df['sequence_name']
    df['motif_score'] = -np.log10(df['p-value'])
    df = df.groupby('region').agg({'motif_score': 'max'}).reset_index()
    df.to_csv(output_path, sep='\t', index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--fimo", required=True, help="Path to FIMO output fimo.tsv")
    parser.add_argument("--out", required=True, help="Output .tsv with max motif scores")
    args = parser.parse_args()
    process_fimo(args.fimo, args.out)