#!/bin/bash --login

#Usage: bash trim.sh <sample_dir> <output_dir> <#_of_cores>

SAMPDIR=$1
OUTDIR=$2
ncores=$3

mkdir -p ${OUTDIR}

for f in "${SAMPDIR}"/*_R1_001.fastq.gz # for each sample

do
    case "$f" in
        *Undetermined*) continue ;;  # skip
    esac
    echo "Processing $f"
    n=${f%%_R1.fastq.gz} # strip part of file name
    trimmomatic PE -threads ${ncores} ${n}_R1.fastq.gz  ${n}_R2.fastq.gz \
    ${OUTDIR}/${n}_R1_trimmed.fastq.gz ${OUTDIR}/${n}_R1_unpaired.fastq.gz \
    ${OUTDIR}/${n}_R2_trimmed.fastq.gz ${OUTDIR}/${n}_R2_unpaired.fastq.gz \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36

done
