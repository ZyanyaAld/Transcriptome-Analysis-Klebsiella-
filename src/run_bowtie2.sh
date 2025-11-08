#!/bin/bash
# Minimal Bowtie2 alignment script
# Inputs : data/trimmed/*_1.trimmed.fastq + *_2.trimmed.fastq
# Index  : data/genome/bowtie2_idx/kleb
# Output : data/alignments/<sample>.sorted.bam (+ .bai)

THREADS=8
IDX="data/genome/bowtie2_idx/kleb"
IN="data/trimmed"
OUT="data/alignments"

mkdir -p "$OUT"

for R1 in "$IN"/*_1.trimmed.fastq; do
  [ -e "$R1" ] || { echo "No *_1.trimmed.fastq in $IN"; exit 1; }
  base=$(basename "${R1%_1.trimmed.fastq}")
  R2="$IN/${base}_2.trimmed.fastq"

  echo "Aligning $base ..."
  bowtie2 -x "$IDX" -1 "$R1" -2 "$R2" -p "$THREADS" \
  | samtools sort -@ "$THREADS" -o "$OUT/${base}.sorted.bam"

  samtools index "$OUT/${base}.sorted.bam"
done

echo "Done. BAMs in $OUT"
