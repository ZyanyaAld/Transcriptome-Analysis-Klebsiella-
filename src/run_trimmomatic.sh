#!/bin/bash

# === Simple Trimmomatic batch script ===
# Klebsiella_AqSCr_RNAseq project
# Author: Marina Mendoza

IN="data/fastq"
OUT="data/trimmed"
ADAPT="data/trimmomatic/TruSeq3-PE.fa"

mkdir -p "$OUT"

for r1 in "$IN"/*_1.fastq; do
    base=$(basename "${r1%_1.fastq}")
    r2="$IN/${base}_2.fastq"

    echo "Trimming $base..."

    trimmomatic PE -threads 8 \
        "$r1" "$r2" \
        "$OUT/${base}_1.trimmed.fastq" "$OUT/${base}_1.unpaired.fastq" \
        "$OUT/${base}_2.trimmed.fastq" "$OUT/${base}_2.unpaired.fastq" \
        ILLUMINACLIP:"$ADAPT":2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50
done

echo "Done! Trimmed files saved in $OUT"
