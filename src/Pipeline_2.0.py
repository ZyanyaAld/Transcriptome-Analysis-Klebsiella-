"""
Pipeline 2.0 Descarga automatizada de datos RNA-seq y genoma de referencia
para *Klebsiella* sp. AQSCr durante la exposición a Cr(VI).

Este script combina la descarga de datos de secuenciación desde SRA
con la obtención estandarizada del genoma de referencia mediante el
módulo `descarga_genoma_referencia.py`.

Características:
    • Descarga automática de los archivos SRR desde un proyecto SRA.
    • Conversión a formato FASTQ mediante SRA Toolkit.
    • Búsqueda y descarga del mejor genoma de referencia (RefSeq preferido).
    • Descompresión automática de archivos FASTA y GFF.
    • Compatible con Linux, macOS y Windows.
"""

import os
import gzip
import shutil
import subprocess
from Bio import Entrez
from descarga_genoma_referencia import fetch_reference_package

# ============================================================================
# CONFIGURACIÓN GENERAL
# ============================================================================
Entrez.email = "marinams@lcg.unam.mx"  # NCBI lo requiere

# Identificador del proyecto SRA (contiene los experimentos RNA-seq)
srp_id = "SRP291413"

# Carpeta principal del proyecto
output_dir = "Klebsiella_AqSCr_RNAseq"

# Parámetros para descubrir el genoma de referencia
organism_query = "Klebsiella sp. AQSCr"   # Nombre del organismo
strain_hint = "AQSCr"                      # Nombre de cepa (opcional pero útil)
bioproject_hint = "PRJNA341863"                     # Ejemplo: "PRJNA341863" (opcional) o NONE

# Rutas del SRA Toolkit (editar según tu sistema)
prefetch_path = "prefetch"        # Usamos un ambiente conda /home/ismadlsh/.conda/envs/bio_informatics
fasterq_dump_path = "fasterq-dump"

# ============================================================================
# CREAR ESTRUCTURA DE DIRECTORIOS
# ============================================================================
dirs = {
    "raw_sra": os.path.join(output_dir, "raw_sra"),
    "fastq": os.path.join(output_dir, "fastq"),
    "genome": os.path.join(output_dir, "genome"),
}

for d in dirs.values():
    os.makedirs(d, exist_ok=True)

print(f" Directorios creados en: {output_dir}")

# ============================================================================
# 1. OBTENER LOS SRR IDs DEL PROYECTO SRA
# ============================================================================
print(f"\n Obteniendo IDs SRR del proyecto SRA: {srp_id}")

handle = Entrez.esearch(db="sra", term=srp_id)
record = Entrez.read(handle)
handle.close()
sra_ids = record["IdList"]

srr_ids = []
for sra_id in sra_ids:
    summary_handle = Entrez.esummary(db="sra", id=sra_id)
    summary = Entrez.read(summary_handle)
    summary_handle.close()
    runs = summary[0]["Runs"]
    for run in runs.split(","):
        start = run.find('acc="') + len('acc="')
        end = run.find('"', start)
        srr_ids.append(run[start:end])

if not srr_ids:
    print("  No se encontraron SRR IDs. Verifica el ID del proyecto SRA.")
else:
    print(f"IDs SRR encontrados: {srr_ids}")

# ============================================================================
# 2. DESCARGA Y CONVERSIÓN DE LOS ARCHIVOS SRR A FASTQ
# ============================================================================
for srr in srr_ids:
    print(f"\n  Descargando SRR {srr} a {dirs['raw_sra']} ...")

    # Saltar descarga si el archivo ya existe
    local_sra = os.path.join(dirs["raw_sra"], f"{srr}.sra")
    if not os.path.exists(local_sra):
        subprocess.run([prefetch_path, srr, "-O", dirs["raw_sra"]], check=True)
    else:
        print(f"  El archivo {local_sra} ya existe — se omite la descarga.")

    # Convertir a FASTQ (solo si no se ha hecho antes)
    fq1 = os.path.join(dirs["fastq"], f"{srr}_1.fastq")
    fq2 = os.path.join(dirs["fastq"], f"{srr}_2.fastq")
    if not (os.path.exists(fq1) or os.path.exists(fq2)):
        print(f"  Convirtiendo {srr} a formato FASTQ ...")
        subprocess.run(
            [fasterq_dump_path, srr, "--split-files", "-O", dirs["fastq"]],
            check=True,
        )
    else:
        print(f"  Los archivos FASTQ ya existen — se omite la conversión.")

# ============================================================================
# 3. DESCUBRIR Y DESCARGAR EL GENOMA DE REFERENCIA
# ============================================================================
print("\n Buscando el mejor ensamblaje disponible (RefSeq preferido)...")

meta, fasta_out, gff_out = fetch_reference_package(
    organism_query=organism_query,
    strain=strain_hint,
    bioproject=bioproject_hint,
    out_dir=dirs["genome"],
    prefer_refseq=True,
    also_download=(),    # Ejemplo: ("protein","cds") si se requiere
    check_md5=False,     # Cambiar a True si deseas verificar MD5
)

print("\n Ensamblaje seleccionado:")
print(f"  Accession:       {meta['accession']}")
print(f"  Organismo:       {meta['organism']}")
print(f"  Nombre:          {meta['assembly_name']}")
print(f"  Estado:          {meta['assembly_status']}")
print(f"  Categoría RefSeq:{meta['refseq_category']}")
print(f"  FTP:             {meta['ftp']}")

# ============================================================================
# 4. DESCOMPRIMIR LOS ARCHIVOS DEL GENOMA
# ============================================================================
def descomprimir_gz(ruta_archivo):
    """Descomprime archivos .gz (FASTA o GFF)."""
    salida = ruta_archivo.replace(".gz", "")
    if not os.path.exists(salida):
        print(f"Descomprimiendo {ruta_archivo} ...")
        with gzip.open(ruta_archivo, "rb") as f_in, open(salida, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    else:
        print(f"{salida} ya existe — se omite la descompresión.")
    return salida


fasta_unzip = descomprimir_gz(fasta_out)
gff_unzip = descomprimir_gz(gff_out)

print(f"\n Genoma y anotaciones listos:\n  FASTA: {fasta_unzip}\n  GFF:   {gff_unzip}")

# ============================================================================
# 5. FINALIZACIÓN DEL PIPELINE
# ============================================================================
print("\n Descargas y conversiones completadas exitosamente.")
print(f"Archivos almacenados en: {os.path.abspath(output_dir)}")