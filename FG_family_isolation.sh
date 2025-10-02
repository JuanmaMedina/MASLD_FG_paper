#!/bin/bash
# ============================================================
# FG_family_isolation.sh
#
# Pipeline to isolate target gene families from UHGG v2.0.2
# Steps:
#   1. Retrieve candidate sequences from UHGG pangenomes by gene name
#   2. Optional exclusion of spurious patterns
#   3. Verify gene annotation against representative proteomes
#   4. Length filtering (min threshold + N50+20 max)
#   5. Concatenate with MetaCyc control sequence
#   6. Multiple sequence alignment (MSA) with MAFFT
#   7. Remove sequences with ≥ GAP_TH gaps in MSA
#   8. Extract clean sequences from original seed
#   9. Build HMM profile with HMMER
#
# Requirements:
#   - seqkit >= 2.0
#   - mafft >= 7.271
#   - HMMER >= 3.4
#   - Python2 with faSomeRecords.py
#   - Python3 with MSA_drop_bygap.py
#
# Usage:
#   bash FG_family_isolation.sh GENE_NAME CONTROL.faa [NEGATIVE_PATTERN]
# Example:
#   bash FG_family_isolation.sh tdcd propk_ecoli_control.faa "bcd_1|bcd_2"
#
# ============================================================

set -euo pipefail

### INPUTS ###
GENE=$1                  # target gene name (e.g. tdcd)
CONTROL=$2               # control enzyme fasta from MetaCyc
NEG_PATTERN=${3:-""}     # optional negative match pattern

### DIRECTORIES ###
PANGENOMES=~/clean_pangenomes_UHGG/genomic_MGYG-HGUT-0*.faa
PROTEOMES=~/UHGG_proteomes/original/MGYG-HGUT-0*.faa

### OUTPUT FILENAMES ###
HEADERS="${GENE}_headers.txt"
PROTEINS="${GENE}_protein_headers.txt"
SEED="${GENE}_seed.faa"
SEED_SHORT="${GENE}_seed_short.faa"
SEED_CONTROL="${GENE}_seed_control.faa"
MSA="${GENE}_MSA.fa"
MSA_CLEAN="${GENE}_MSA_clean.fa"
CLEAN_SEQS="${GENE}_clean.faa"
HMM="${GENE}.hmm"

### PARAMETERS ###
MIN_LEN=100              # minimum protein length
GAP_TH=0.3               # max allowed fraction of gaps per sequence

echo "==== STEP 1. Retrieve candidate headers from UHGG pangenomes ===="
rm -f $HEADERS $SEED $PROTEINS
for P in $PANGENOMES; do
    grep "^>" "$P" | grep -i " $GENE" >> $HEADERS
done

if [[ -n "$NEG_PATTERN" ]]; then
    echo "Removing unwanted matches ($NEG_PATTERN)"
    grep -Ev "$NEG_PATTERN" $HEADERS > tmp && mv tmp $HEADERS
fi

echo "Headers retrieved: $(wc -l < $HEADERS)"

echo "==== STEP 2. Verify annotation against representative proteomes ===="
for PRO in $PROTEOMES; do
    grep -f <(awk '{print $2}' $HEADERS) $PRO >> $PROTEINS || true
done
echo "Proteome cross-check done."

echo "==== STEP 3. Extract seed sequences ===="
for P in $PANGENOMES; do
    python2 ~/bin/faSomeRecords.py --fasta "$P" --list $HEADERS --stdout \
        | grep -v "No sequences found" >> $SEED || true
done
echo "Seed sequences: $(grep -c "^>" $SEED)"

echo "==== STEP 4. Filter sequences by length ===="
LENGTH_MAX=$(seqkit stats -a -b $SEED | tail -n1 | awk '{print $13}' | sed 's/,//' | awk '{print $1+20}')
seqkit seq -g -m $MIN_LEN -M $LENGTH_MAX $SEED > $SEED_SHORT
echo "After length filtering: $(grep -c "^>" $SEED_SHORT)"

echo "==== STEP 5. Concatenate with MetaCyc control ===="
cat $SEED_SHORT $CONTROL > $SEED_CONTROL

echo "==== STEP 6. Multiple Sequence Alignment with MAFFT ===="
mafft --localpair --maxiterate 1000 --thread -1 --quiet $SEED_CONTROL > $MSA
rm $SEED_CONTROL
echo "MSA complete."

echo "==== STEP 7. Remove sequences with ≥${GAP_TH} gaps ===="
python3 ~/bin/MSA_drop_bygap.py $MSA $MSA_CLEAN $GAP_TH
echo "Cleaned MSA: $(grep -c '^>' $MSA_CLEAN)"

echo "==== STEP 8. Extract clean sequences from original seed ===="
grep "^>" $MSA_CLEAN | sed 's/>//' > clean_ids.txt
python2 ~/bin/faSomeRecords.py --fasta $SEED --list clean_ids.txt --stdout > $CLEAN_SEQS
rm clean_ids.txt
echo "Final clean sequences: $(grep -c '^>' $CLEAN_SEQS)"

echo "==== STEP 9. Build HMM profile with HMMER ===="
hmmbuild $HMM $MSA_CLEAN
echo "HMM built: $HMM"

echo "Pipeline finished successfully!"
