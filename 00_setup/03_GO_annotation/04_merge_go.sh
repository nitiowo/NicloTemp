#!/bin/bash

# Merge NCBI GO annotations with eggNOG-mapper GO annotations

set -e

REPO="$(cd "$(dirname "$0")/../.." && pwd)"
STEP01_DIR="${REPO}/00_setup/03_GO_annotation/output/step01_ncbi"
STEP03_DIR="${REPO}/00_setup/03_GO_annotation/output/step03_parsed"
MERGED_DIR="${REPO}/00_setup/03_GO_annotation/output/step04_merged_go"
mkdir -p "$MERGED_DIR"

# ---- Merge Btru GO ----

BTRU_NCBI="${STEP01_DIR}/btru_ncbi_go.tsv"
BTRU_EMAP="${STEP03_DIR}/btru_emapper_go.tsv"

{
    cut -f1,2 "$BTRU_NCBI"
    cat "$BTRU_EMAP"
} | sort -u > "${MERGED_DIR}/btru_gene2go.tsv"

# Coverage stats
BTRU_TOTAL=26279
BTRU_NCBI_GENES=$(cut -f1 "$BTRU_NCBI" | sort -u | wc -l)
BTRU_EMAP_GENES=$(cut -f1 "$BTRU_EMAP" | sort -u | wc -l)
BTRU_MERGED_GENES=$(cut -f1 "${MERGED_DIR}/btru_gene2go.tsv" | sort -u | wc -l)
BTRU_MERGED_PAIRS=$(wc -l < "${MERGED_DIR}/btru_gene2go.tsv")

# Genes only found via eggNOG
BTRU_NEW=$(comm -23 \
    <(cut -f1 "$BTRU_EMAP" | sort -u) \
    <(cut -f1 "$BTRU_NCBI" | sort -u) | wc -l)

# Summary output (save to file?)
echo "  NCBI:   ${BTRU_NCBI_GENES} genes"
echo "  eggNOG: ${BTRU_EMAP_GENES} genes"
echo "  New genes from eggNOG: ${BTRU_NEW}"
echo "  Merged: ${BTRU_MERGED_GENES} / ${BTRU_TOTAL} genes ($(echo "scale=1; ${BTRU_MERGED_GENES}*100/${BTRU_TOTAL}" | bc)%)"
echo "  Total gene-GO pairs: ${BTRU_MERGED_PAIRS}"

# ---- Merge Shae GO ----

SHAE_NCBI="${STEP01_DIR}/shae_ncbi_go.tsv"
SHAE_EMAP="${STEP03_DIR}/shae_emapper_go.tsv"

{
    cut -f1,2 "$SHAE_NCBI"
    cat "$SHAE_EMAP"
} | sort -u > "${MERGED_DIR}/shae_gene2go.tsv"

SHAE_TOTAL=9428
SHAE_NCBI_GENES=$(cut -f1 "$SHAE_NCBI" | sort -u | wc -l)
SHAE_EMAP_GENES=$(cut -f1 "$SHAE_EMAP" | sort -u | wc -l)
SHAE_MERGED_GENES=$(cut -f1 "${MERGED_DIR}/shae_gene2go.tsv" | sort -u | wc -l)
SHAE_MERGED_PAIRS=$(wc -l < "${MERGED_DIR}/shae_gene2go.tsv")

SHAE_NEW=$(comm -23 \
    <(cut -f1 "$SHAE_EMAP" | sort -u) \
    <(cut -f1 "$SHAE_NCBI" | sort -u) | wc -l)

echo "  NCBI:   ${SHAE_NCBI_GENES} genes"
echo "  eggNOG: ${SHAE_EMAP_GENES} genes"
echo "  New genes from eggNOG: ${SHAE_NEW}"
echo "  Merged: ${SHAE_MERGED_GENES} / ${SHAE_TOTAL} genes ($(echo "scale=1; ${SHAE_MERGED_GENES}*100/${SHAE_TOTAL}" | bc)%)"
echo "  Total gene-GO pairs: ${SHAE_MERGED_PAIRS}"

#Keep the detailed NCBI versions
cp "$BTRU_NCBI" "${MERGED_DIR}/btru_ncbi_go_detail.tsv"
cp "$SHAE_NCBI" "${MERGED_DIR}/shae_ncbi_go_detail.tsv"
