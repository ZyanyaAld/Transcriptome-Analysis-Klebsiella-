#!/bin/bash
# Quality control for RAW (untrimmed) RNA-seq reads
# Tools: FastQC + MultiQC
# Input : data/fastq/*.fastq
# Output: data/qc/ (individual FastQC reports + summary MultiQC report)

IN="data/fastq"
OUT="data/qc"
THREADS=8

echo "Running FastQC on raw reads..."
mkdir -p "$OUT"

fastqc -t "$THREADS" -o "$OUT" "$IN"/*.fastq

echo "Aggregating results with MultiQC..."
conda init
conda activate multiqc
multiqc "$OUT" -o "$OUT"

echo "QC complete. Reports saved in: $OUT"
