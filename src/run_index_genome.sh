#!/usr/bin/env bash

# Indexa el GENOMA (no el CDS) para Bowtie2.
# Estructura esperada:
#   data/genome/GCF_008452795.1_ASM845279v1_genomic.fna

GENOME_FA="data/genome/GCF_008452795.1_ASM845279v1_genomic.fna"
IDX_DIR="data/genome/bowtie2_idx"
IDX_PREFIX="${IDX_DIR}/kleb"

mkdir -p "${IDX_DIR}"

echo "Indexando genoma con bowtie2-build..."
bowtie2-build "${GENOME_FA}" "${IDX_PREFIX}"

echo "Checando encabezados del índice:"
bowtie2-inspect -n "${IDX_PREFIX}" | head
echo "Índice listo en ${IDX_DIR}"
