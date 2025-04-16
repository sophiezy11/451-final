#!/bin/bash

# Activate environment if needed
# conda activate tfbs

# Run FIMO for positive sequences
fimo --thresh 1e-4 --oc fimo_outputs/fimo_pos_out data/CTCF.meme sequences/positives.fa

# Run FIMO for negative sequences
fimo --thresh 1e-4 --oc fimo_outputs/fimo_neg_out data/CTCF.meme sequences/negatives.fa
