#!/bin/bash
# Build Bowtie2 index from genome in data/genome/
# Output: data/genome/bowtie2_idx/kleb.*.bt2

mkdir -p data/genome/bowtie2_idx

# 1) Search for uncompressed FASTA
FASTA=""
for f in data/genome/*_genomic.fna data/genome/*.fna data/genome/*.fa; do
  if [ -f "$f" ]; then FASTA="$f"; break; fi
done

# 2) If there is no uncompressed file, then use a .gz (makes an uncompressed copy)
if [ -z "$FASTA" ]; then
  for g in data/genome/*_genomic.fna.gz data/genome/*.fna.gz data/genome/*.fa.gz; do
    if [ -f "$g" ]; then
      gunzip -c "$g" > data/genome/genome_temp.fna
      FASTA="data/genome/genome_temp.fna"
      break
    fi
  done
fi

# 3) Exits if nothing is found
if [ -z "$FASTA" ]; then
  echo "No FASTA found in data/genome/"
  exit 1
fi

echo "Indexing genome: $FASTA"
bowtie2-build "$FASTA" data/genome/bowtie2_idx/kleb
echo "Done."
