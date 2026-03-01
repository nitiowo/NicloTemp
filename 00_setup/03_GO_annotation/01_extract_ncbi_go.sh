#!/bin/bash

# Extract existing GO annotations from NCBI
# Shae needs eggnog


set -e

REPO="$(cd "$(dirname "$0")/../.." && pwd)"
OUTDIR="${REPO}/00_setup/03_GO_annotation/output/step01_ncbi"
mkdir -p "$OUTDIR"

# ---- Reference Paths ----
BTRU_GTF="${REPO}/Data/reference_data/Btru/GCA_021962125.1/Btru_genomic.gtf"
SHAE_GTF="${REPO}/Data/reference_data/Shae/GCF_000699445.3/Shae_genomic.gtf"
SHAE_GAF="${REPO}/Data/reference_data/Shae/GCF_000699445.3/GCF_000699445.3_UoM_Shae.V3_gene_ontology.gaf"

# ---- Extract btru GO from GTF ----

awk -F'\t' '$3 == "CDS" && /go_function|go_process|go_component/ {
    match($9, /gene_id "([^"]+)"/, g)
    id = g[1]
    n = split($9, arr, ";")
    for (i = 1; i <= n; i++) {
        if (arr[i] ~ /go_function/) {
            match(arr[i], /go_function "([^|]+)\|([0-9]+)/, m)
            if (m[2]) print id "\tGO:" m[2] "\tF\t" m[1]
        } else if (arr[i] ~ /go_process/) {
            match(arr[i], /go_process "([^|]+)\|([0-9]+)/, m)
            if (m[2]) print id "\tGO:" m[2] "\tP\t" m[1]
        } else if (arr[i] ~ /go_component/) {
            match(arr[i], /go_component "([^|]+)\|([0-9]+)/, m)
            if (m[2]) print id "\tGO:" m[2] "\tC\t" m[1]
        }
    }
}' "$BTRU_GTF" | sort -u > "${OUTDIR}/btru_ncbi_go.tsv"

BTRU_GENES=$(cut -f1 "${OUTDIR}/btru_ncbi_go.tsv" | sort -u | wc -l)
BTRU_PAIRS=$(wc -l < "${OUTDIR}/btru_ncbi_go.tsv")

# ---- Shae GeneID-to-gene_id Mapping ----

awk -F'\t' '$3 == "gene" {
    match($9, /gene_id "([^"]+)"/, g)
    match($9, /GeneID:([0-9]+)/, d)
    if (g[1] && d[1]) print d[1] "\t" g[1]
}' "$SHAE_GTF" > "${OUTDIR}/shae_ncbiid_to_geneid.tsv"

# ---- Extract Shae GO from GAF ----

awk -F'\t' 'BEGIN{OFS="\t"}
    FNR == NR { map[$1] = $2; next }
    !/^!/ {
        geneid = $2
        go = $5
        aspect = $9
        if (geneid in map) print map[geneid], go, aspect, ""
    }
' "${OUTDIR}/shae_ncbiid_to_geneid.tsv" "$SHAE_GAF" | sort -u > "${OUTDIR}/shae_ncbi_go.tsv"

SHAE_GENES=$(cut -f1 "${OUTDIR}/shae_ncbi_go.tsv" | sort -u | wc -l)
SHAE_PAIRS=$(wc -l < "${OUTDIR}/shae_ncbi_go.tsv")

# ---- Build protein_id to gene_id Mappings (for eggNog) ----

awk -F'\t' '$3 == "CDS" && /protein_id/ {
    match($9, /gene_id "([^"]+)"/, g)
    match($9, /protein_id "([^"]+)"/, p)
    if (g[1] && p[1]) print p[1] "\t" g[1]
}' "$BTRU_GTF" | sort -u > "${OUTDIR}/btru_protid_to_geneid.tsv"

awk -F'\t' '$3 == "CDS" && /protein_id/ {
    match($9, /gene_id "([^"]+)"/, g)
    match($9, /protein_id "([^"]+)"/, p)
    if (g[1] && p[1]) print p[1] "\t" g[1]
}' "$SHAE_GTF" | sort -u > "${OUTDIR}/shae_protid_to_geneid.tsv"
