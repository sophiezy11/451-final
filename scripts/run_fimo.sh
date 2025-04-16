# Script: run_fimo.sh
# Purpose: Run FIMO on a given FASTA file using a MEME motif file
# Usage: bash scripts/run_fimo.sh CTCF.meme input.fa output_dir

#!/bin/bash

set -euo pipefail

if [ "$#" -ne 3 ]; then
  echo "Usage: $0 <motif_file.meme> <input.fa> <output_dir>"
  exit 1
fi

MOTIF_FILE=$1
INPUT_FASTA=$2
OUTPUT_DIR=$3

fimo --thresh 1e-4 --verbosity 1 --oc "$OUTPUT_DIR" "$MOTIF_FILE" "$INPUT_FASTA"
echo "FIMO scan complete. Results saved to $OUTPUT_DIR"
