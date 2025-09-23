#!/bin/bash
# Usage: bash trim.sh <R1.fastq.gz> <output_dir> <ncores>

R1="$1"
OUTDIR="$2"
ncores="$3"

mkdir -p "$OUTDIR"

# Derive R2 filename automatically
R2="${R1/%_R1_001.fastq.gz/_R2_001.fastq.gz}"

# Extract basename without path
f=$(basename "$R1")

# Remove extensions
base="${f%.fastq.gz}"

# Shorten name by removing fields 4 and 6 (split by _)
IFS=_ read -r -a p <<< "$base"
newbase="${p[0]}_${p[1]}_${p[2]}"

# Define output files with the new short base
R1_PAIRED="${OUTDIR}/${newbase}_R1_paired.fastq.gz"
R1_UNPAIRED="${OUTDIR}/${newbase}_R1_unpaired.fastq.gz"
R2_PAIRED="${OUTDIR}/${newbase}_R2_paired.fastq.gz"
R2_UNPAIRED="${OUTDIR}/${newbase}_R2_unpaired.fastq.gz"

# Run Trimmomatic
trimmomatic PE -threads "$ncores" \
    "$R1" "$R2" \
    "$R1_PAIRED" "$R1_UNPAIRED" \
    "$R2_PAIRED" "$R2_UNPAIRED" \
    ILLUMINACLIP:"$ADAPTERS":2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
