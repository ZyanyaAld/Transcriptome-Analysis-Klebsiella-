#!/usr/bin/env bash

# Cuantificacion con featureCounts a nivel de CDS por locus_tag.
# Usa los BAM ya alineados y ordenados:
#   data/alignments/*.sorted.bam

GFF="data/genome/GCF_008452795.1_ASM845279v1_genomic.gff"
OUTDIR="data/counts"
OUT="${OUTDIR}/counts_final.txt"

mkdir -p "${OUTDIR}"

echo "Ejecutando featureCounts..."
featureCounts -T 8 -p -B -C \
  -t CDS -g locus_tag \
  -a "${GFF}" \
  -o "${OUT}" \
  data/alignments/*.sorted.bam

echo "[OK] Conteos escritos en ${OUT}"
head -n 10 "${OUT}"
