#!/bin/bash

# === QC for trimmed FASTQ files ===
# Klebsiella_AqSCr_RNAseq project
# Author: Marina Mendoza

IN="data/trimmed"
OUT="data/qc_trimmed"

mkdir -p "$OUT"

echo "Running FastQC on trimmed reads..."
fastqc -t 8 -o "$OUT" "$IN"/*.trimmed.fastq

echo "Summarizing with MultiQC..."
conda activate multiqc
multiqc "$OUT" -o "$OUT"

echo "Done! Reports in: $OUT"
