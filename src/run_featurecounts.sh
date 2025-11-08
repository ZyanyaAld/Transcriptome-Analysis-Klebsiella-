#!/bin/bash
# Minimal featureCounts for paired-end bacterial RNA-seq
# Input : data/alignments/*.sorted.bam
# GFF   : data/genome/*.gff  (auto-picks the first)
# Output: data/counts/counts.txt (+ summary)

GFF=$(ls data/genome/*.gff | head -n1)
OUTDIR="data/counts"
mkdir -p "$OUTDIR"

# Count at gene level (common in prokaryotes: -t gene -g ID)
# If your GFF lacks 'gene' features, switch to -t CDS -g Parent
conda activate subread
featureCounts -T 8 -p -B -C \
  -t gene -g ID \
  -a "$GFF" \
  -o "$OUTDIR/counts.txt" data/alignments/*.sorted.bam

echo "Done. Matrix: $OUTDIR/counts.txt ; Summary: $OUTDIR/counts.txt.summary"
