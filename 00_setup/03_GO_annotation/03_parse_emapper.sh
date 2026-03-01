#!/bin/bash

# Parse eggNOG-mapper .annotations files into gene_id > GO mappings

set -e

REPO="$(cd "$(dirname "$0")/../.." && pwd)"
STEP01_DIR="${REPO}/00_setup/03_GO_annotation/output/step01_ncbi"
EMAP_DIR="${REPO}/00_setup/03_GO_annotation/output/step02_emapper"
PARSED_DIR="${REPO}/00_setup/03_GO_annotation/output/step03_parsed"
mkdir -p "$PARSED_DIR"

BTRU_ANNOT="${EMAP_DIR}/btru_emapper.emapper.annotations"
SHAE_ANNOT="${EMAP_DIR}/shae_emapper.emapper.annotations"

# ---- Parse Btru eggNOG ----
# eggNOG format col 1: query, col 10: GO terms

awk -F'\t' 'BEGIN{OFS="\t"}
    FNR == NR { map[$1] = $2; next }
    !/^#/ && $10 != "-" && $10 != "" {
        prot = $1
        if (prot in map) {
            gene = map[prot]
            n = split($10, gos, ",")
            for (i = 1; i <= n; i++) {
                gsub(/^ +| +$/, "", gos[i])
                if (gos[i] ~ /^GO:/) print gene "\t" gos[i]
            }
        }
    }
' "${STEP01_DIR}/btru_protid_to_geneid.tsv" "$BTRU_ANNOT" \
    | sort -u > "${PARSED_DIR}/btru_emapper_go.tsv"

BTRU_GENES=$(cut -f1 "${PARSED_DIR}/btru_emapper_go.tsv" | sort -u | wc -l)
BTRU_PAIRS=$(wc -l < "${PARSED_DIR}/btru_emapper_go.tsv")

# ---- Parse Shae eggNOG ----

awk -F'\t' 'BEGIN{OFS="\t"}
    FNR == NR { map[$1] = $2; next }
    !/^#/ && $10 != "-" && $10 != "" {
        prot = $1
        if (prot in map) {
            gene = map[prot]
            n = split($10, gos, ",")
            for (i = 1; i <= n; i++) {
                gsub(/^ +| +$/, "", gos[i])
                if (gos[i] ~ /^GO:/) print gene "\t" gos[i]
            }
        }
    }
' "${STEP01_DIR}/shae_protid_to_geneid.tsv" "$SHAE_ANNOT" \
    | sort -u > "${PARSED_DIR}/shae_emapper_go.tsv"

SHAE_GENES=$(cut -f1 "${PARSED_DIR}/shae_emapper_go.tsv" | sort -u | wc -l)
SHAE_PAIRS=$(wc -l < "${PARSED_DIR}/shae_emapper_go.tsv")

# ---- Parse Btru eggNOG KEGG KO Assignments ----
# eggNOG format col 12: KEGG KOs

awk -F'\t' 'BEGIN{OFS="\t"}
    FNR == NR { map[$1] = $2; next }
    !/^#/ && $12 != "-" && $12 != "" {
        prot = $1
        if (prot in map) {
            gene = map[prot]
            n = split($12, kos, ",")
            for (i = 1; i <= n; i++) {
                gsub(/^ +| +$/, "", kos[i])
                if (kos[i] ~ /^K[0-9]/) print gene "\t" kos[i]
            }
        }
    }
' "${STEP01_DIR}/btru_protid_to_geneid.tsv" "$BTRU_ANNOT" \
    | sort -u > "${PARSED_DIR}/btru_emapper_kegg.tsv"

BTRU_KEGG_GENES=$(cut -f1 "${PARSED_DIR}/btru_emapper_kegg.tsv" | sort -u | wc -l)
BTRU_KEGG_PAIRS=$(wc -l < "${PARSED_DIR}/btru_emapper_kegg.tsv")

# ---- Parse Shae eggNOG KEGG KO Assignments ----

awk -F'\t' 'BEGIN{OFS="\t"}
    FNR == NR { map[$1] = $2; next }
    !/^#/ && $12 != "-" && $12 != "" {
        prot = $1
        if (prot in map) {
            gene = map[prot]
            n = split($12, kos, ",")
            for (i = 1; i <= n; i++) {
                gsub(/^ +| +$/, "", kos[i])
                if (kos[i] ~ /^K[0-9]/) print gene "\t" kos[i]
            }
        }
    }
' "${STEP01_DIR}/shae_protid_to_geneid.tsv" "$SHAE_ANNOT" \
    | sort -u > "${PARSED_DIR}/shae_emapper_kegg.tsv"

SHAE_KEGG_GENES=$(cut -f1 "${PARSED_DIR}/shae_emapper_kegg.tsv" | sort -u | wc -l)
SHAE_KEGG_PAIRS=$(wc -l < "${PARSED_DIR}/shae_emapper_kegg.tsv")
