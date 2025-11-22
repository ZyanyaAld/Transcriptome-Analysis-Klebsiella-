#!/bin/bash
IN="data/alignments"
OUT="data/alignments/stats"
mkdir -p "$OUT"
for bam in "$IN"/*.sorted.bam; do
  base=$(basename "${bam%.sorted.bam}")
  samtools flagstat "$bam" > "$OUT/${base}.flagstat.txt"
done
echo "Done. Stats in $OUT"
