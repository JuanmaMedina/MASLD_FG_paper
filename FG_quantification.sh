#!/bin/bash

# Gene abundance estimation pipeline using DIAMOND
# Usage:
#   ./FG_quantification.sh sample.fna.gz FG_family.dmnd threads

# Input variables
SAMPLE_FASTA=$1   # Input reads in FASTA format (can be gzipped)
FGC_DB=$2         # DIAMOND database for the gene family
THREADS=$3        # Number of threads

# Output directory (edit as needed)
OUTPUT_DIR="./diamond_out"
mkdir -p "$OUTPUT_DIR"

# Derive sample and gene names from input
SAMPLE=$(basename "$SAMPLE_FASTA" _cat.fna.gz)
GENE=$(basename "$FGC_DB" _FINAL_SEQS.dmnd)

# Output files
OUT_TOTAL="$OUTPUT_DIR/${GENE}_k1_total.tsv"

echo "Running DIAMOND for sample: $SAMPLE against gene family: $GENE"

# Run DIAMOND translated search
DIAMOND_RESULTS=$(diamond blastx --db "$FGC_DB" --query "$SAMPLE_FASTA" --threads "$THREADS" --max-hsps 1 --evalue 1e-10 -k 1 --outfmt 6 sseqid --quiet)

# Count total matches
TOTAL=$(echo "$DIAMOND_RESULTS" | wc -l)
echo -e "${SAMPLE}\t${TOTAL}" >> "$OUT_TOTAL"

echo "Finished. Results written to $OUT_TOTAL"
