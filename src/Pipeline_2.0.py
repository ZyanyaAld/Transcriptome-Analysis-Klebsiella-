import os
from Bio import Entrez
import subprocess
import urllib.request
import gzip
import shutil

# ===============================
# CONFIGURATION
# ===============================
Entrez.email = "zyanyava@lcg.unam.mx"  # Replace with your email

# SRA project ID for GSE160968
srp_id = "SRP291413"

output_dir = "Klebsiella_AqSCr_RNAseq"

# Genome URLs (strain AqSCr / SAMN05730075)
fasta_url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/008/452/795/GCA_008452795.1_MJDM00000000/GCA_008452795.1_MJDM00000000_genomic.fna.gz"
gff_url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/008/452/795/GCA_008452795.1_MJDM00000000/GCA_008452795.1_MJDM00000000_genomic.gff.gz"

# Full paths to SRA Toolkit executables
prefetch_path = r"C:\Program Files\sratoolkit.3.2.1-win64\bin\prefetch.exe"        # <-- UPDATE THIS
fasterq_dump_path = r"C:\Program Files\sratoolkit.3.2.1-win64\bin\fasterq-dump.exe" # <-- UPDATE THIS

# ===============================
# CREATE DIRECTORY STRUCTURE
# ===============================
dirs = {
    "raw_sra": os.path.join(output_dir, "raw_sra"),
    "fastq": os.path.join(output_dir, "fastq"),
    "genome": os.path.join(output_dir, "genome")
}

for d in dirs.values():
    os.makedirs(d, exist_ok=True)

print(f"Created project directories in {output_dir}")

# ===============================
# 1. GET SRR IDs FROM SRA PROJECT
# ===============================
print(f"Fetching SRR IDs from SRA project {srp_id}...")

handle = Entrez.esearch(db="sra", term=srp_id)
record = Entrez.read(handle)
handle.close()

sra_ids = record['IdList']

srr_ids = []
for sra_id in sra_ids:
    summary_handle = Entrez.esummary(db="sra", id=sra_id)
    summary = Entrez.read(summary_handle)
    summary_handle.close()
    runs = summary[0]['Runs']  # XML-like string
    # Extract SRR accession numbers
    for run in runs.split(','):
        start = run.find('acc="') + len('acc="')
        end = run.find('"', start)
        srr_ids.append(run[start:end])

if not srr_ids:
    print("No SRR IDs found. Check the SRP ID.")
else:
    print(f"Clean SRR IDs: {srr_ids}")

# ===============================
# 2. DOWNLOAD SRR FILES AND CONVERT TO FASTQ
# ===============================
for srr in srr_ids:
    print(f"\nDownloading SRR {srr} to raw_sra...")
    subprocess.run([prefetch_path, srr, "-O", dirs["raw_sra"]], check=True)

    print(f"Converting {srr} to FASTQ...")
    subprocess.run([fasterq_dump_path, srr, "--split-files", "-O", dirs["fastq"]], check=True)

# ===============================
# 3. DOWNLOAD GENOME ASSEMBLY
# ===============================
print("\nDownloading genome assembly and annotation...")
fasta_out = os.path.join(dirs["genome"], "genome.fna.gz")
gff_out = os.path.join(dirs["genome"], "genome.gff.gz")

# Download FASTA
if not os.path.exists(fasta_out):
    print("Downloading genome FASTA...")
    urllib.request.urlretrieve(fasta_url, fasta_out)
else:
    print("Genome FASTA already exists, skipping.")

# Download GFF
if not os.path.exists(gff_out):
    print("Downloading genome GFF...")
    urllib.request.urlretrieve(gff_url, gff_out)
else:
    print("Genome GFF already exists, skipping.")

# ===============================
# 4. UNZIP GENOME FILES
# ===============================
def unzip_gz(file_path):
    out_file = file_path.replace(".gz","")
    if not os.path.exists(out_file):
        print(f"Unzipping {file_path}...")
        with gzip.open(file_path, 'rb') as f_in:
            with open(out_file, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
    else:
        print(f"{out_file} already exists, skipping unzip.")

unzip_gz(fasta_out)
unzip_gz(gff_out)

print("\nAll downloads and conversions completed successfully!")