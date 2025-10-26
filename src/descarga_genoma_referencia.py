# -*- coding: utf-8 -*-
"""
Descarga de genoma de referencia desde NCBI Assembly (estándar reproducible)

Uso típico desde otro script:

    from descarga_genoma_referencia import fetch_reference_package

    meta, fasta, gff = fetch_reference_package(
        organism_query="Klebsiella sp. AQSCr",
        strain="AQSCr",
        bioproject=None,
        out_dir="results/genome",
        prefer_refseq=True,
        also_download=("protein", "cds"),  # opcional
        check_md5=True                     # opcional
    )

Requisitos:
    - biopython

Notas:
    - Este módulo consulta la base de datos Assembly (no nuccore) y prioriza RefSeq (GCF_) y
      ensamblajes con estatus "Complete Genome"/"Chromosome".
    - Descarga *_genomic.fna.gz y *_genomic.gff.gz (mínimo reproducible), además de protein/cds si se pide.
    - Escribe un archivo ASSEMBLY_PROVENANCE.txt con metadatos para trazabilidad.
"""

from __future__ import annotations
from typing import Dict, List, Optional, Tuple

import datetime
import os
import re
import time
import hashlib
import urllib.request
from urllib.error import HTTPError, URLError

from Bio import Entrez

# ──────────────────────────────────────────────────────────────────────────────
# CONFIGURACIÓN GLOBAL
# ──────────────────────────────────────────────────────────────────────────────
# Sustituye por un correo real para cumplir políticas de NCBI (o exporta NCBI_EMAIL)
Entrez.email = os.environ.get("NCBI_EMAIL", "marinams@lcg.unam.mx")
# Soporte opcional para mayor cuota de peticiones (exporta NCBI_API_KEY si la tienes)
_api_key = os.environ.get("NCBI_API_KEY")
if _api_key:
    Entrez.api_key = _api_key

# ──────────────────────────────────────────────────────────────────────────────
# UTILIDADES: HTTP con reintentos, MD5, parseo de md5checksums.txt
# ──────────────────────────────────────────────────────────────────────────────

def _http_get(url: str, timeout: int = 60, tries: int = 5, backoff: float = 0.7) -> bytes:
    """GET con reintentos exponenciales (429/5xx/errores de red)."""
    last_exc = None
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as r:
                return r.read()
        except (HTTPError, URLError, TimeoutError) as e:
            last_exc = e
            time.sleep(backoff * (2 ** i))
    raise last_exc

def _download_file(url: str, out_path: str, overwrite: bool = False) -> None:
    """Descarga con reintentos; idempotente por defecto (no sobrescribe si existe con tamaño>0)."""
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    if (not overwrite) and os.path.exists(out_path) and os.path.getsize(out_path) > 0:
        return
    data = _http_get(url)
    with open(out_path, "wb") as f:
        f.write(data)

def _parse_md5_file_text(md5_text: str) -> Dict[str, str]:
    """Parses md5checksums.txt (formato: <md5> <espacio> <ruta/archivo>)."""
    md5map: Dict[str, str] = {}
    for ln in md5_text.splitlines():
        ln = ln.strip()
        if not ln or ln.startswith("#"):
            continue
        parts = ln.split()
        if len(parts) >= 2:
            md5 = parts[0].lower()
            fname = parts[-1].split("/")[-1]  # tomar nombre, ignorar path
            if len(md5) == 32 and fname:
                md5map[fname] = md5
    return md5map

def _md5sum(path: str) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

# ──────────────────────────────────────────────────────────────────────────────
# SELECCIÓN DE ENSAMBLAJE EN NCBI ASSEMBLY
# ──────────────────────────────────────────────────────────────────────────────

def _score_assembly(doc) -> Tuple[int, int, int, datetime.datetime]:
    """Asigna una puntuación objetiva a un ensamblaje.

    Criterios (en orden):
      1) RefSeq_category: reference > representative > none
      2) AssemblyStatus: Complete Genome > Chromosome > Scaffold > Contig
      3) Namespace: GCF_ (RefSeq) > GCA_ (GenBank)
      4) Fecha de liberación (más nuevo mejor)
    """
    refcat = (doc.get("RefSeq_category") or "").lower()
    ref_score = {"reference genome": 3, "representative genome": 2}.get(refcat, 0)

    status = (doc.get("AssemblyStatus") or "").lower()
    lvl_score = {"complete genome": 4, "chromosome": 3, "scaffold": 2, "contig": 1}.get(status, 0)

    acc = doc.get("AssemblyAccession", "")
    is_refseq = 1 if acc.startswith("GCF_") else 0

    # fecha tipo "2020/09/28"
    try:
        dt = datetime.datetime.strptime(doc.get("AsmReleaseDate_GenBank", "1900/01/01"), "%Y/%m/%d")
    except Exception:
        dt = datetime.datetime(1900, 1, 1)
    return (ref_score, lvl_score, is_refseq, dt)

def _entrez_read_with_retry(fn, max_tries: int = 5, base_sleep: float = 0.5):
    """Wrapper de Entrez.read con reintentos/backoff."""
    last_exc = None
    for i in range(max_tries):
        try:
            h = fn()
            rec = Entrez.read(h)
            h.close()
            return rec
        except Exception as e:
            last_exc = e
            time.sleep(base_sleep * (2 ** i))
    raise last_exc

def find_best_assembly(
    organism_query: str,
    strain: Optional[str] = None,
    bioproject: Optional[str] = None,
    max_hits: int = 100,
    prefer_refseq: bool = True,
) -> Dict[str, str]:
    """Busca y selecciona el mejor ensamblaje en NCBI Assembly.

    Retorna un diccionario con:
      - accession, ftp, organism, assembly_name, assembly_status,
        refseq_category, taxid, bioproject
    """
    # Construir query
    terms = [organism_query, "latest[filter]"]
    if bioproject:
        terms.append(f"{bioproject}[BioProject]")
    q = " AND ".join(terms)

    # Buscar IDs (con reintentos)
    es = _entrez_read_with_retry(lambda: Entrez.esearch(db="assembly", term=q, retmax=max_hits))
    ids = es.get("IdList", [])

    if not ids:
        # reintentar más laxo
        es = _entrez_read_with_retry(lambda: Entrez.esearch(db="assembly", term=organism_query, retmax=max_hits))
        ids = es.get("IdList", [])
        if not ids:
            raise RuntimeError(f"No se encontraron ensamblajes para: {organism_query}")

    # Resumen completo
    recs = _entrez_read_with_retry(lambda: Entrez.esummary(db="assembly", id=",".join(ids), report="full"))
    docs = recs["DocumentSummarySet"]["DocumentSummary"]

    # Filtro opcional por strain
    if strain:
        s = strain.lower()
        def has_strain(d) -> bool:
            fields = [
                d.get("SpeciesName", ""),
                d.get("AssemblyName", ""),
                d.get("Organism", ""),
                d.get("Biosource", ""),
                d.get("SubmitterOrganization", ""),
                d.get("WGSProject", ""),
            ]
            return any(s in (x or "").lower() for x in fields)
        cand = [d for d in docs if has_strain(d)]
        if cand:
            docs = cand

    # Preferir RefSeq si está disponible
    if prefer_refseq:
        refseq_docs = [d for d in docs if (d.get("AssemblyAccession", "").startswith("GCF_"))]
        if refseq_docs:
            docs = refseq_docs

    # Elegir el mejor según puntuación
    best = max(docs, key=_score_assembly)

    acc = best.get("AssemblyAccession")
    ftp = best.get("FtpPath_RefSeq") or best.get("FtpPath_GenBank")
    if not ftp:
        raise RuntimeError(f"El ensamblaje {acc} no tiene ruta FTP en el resumen.")

    return {
        "accession": acc,
        "ftp": ftp.replace("ftp://", "https://"),  # descarga por HTTPS
        "organism": best.get("Organism"),
        "assembly_name": best.get("AssemblyName"),
        "assembly_status": best.get("AssemblyStatus"),
        "refseq_category": best.get("RefSeq_category") or "None",
        "taxid": best.get("Taxid"),
        "bioproject": best.get("BioprojectAccn"),
    }

# ──────────────────────────────────────────────────────────────────────────────
# DESCARGA DE PAQUETE DE REFERENCIA (HTTPS + md5checksums + reintentos)
# ──────────────────────────────────────────────────────────────────────────────

def download_reference_genome(
    ftp_url: str,
    output_dir: str,
    also_download: Tuple[str, ...] = (),  # e.g., ("protein", "cds", "gbff")
    check_md5: bool = True,
) -> Dict[str, str]:
    """
    Descarga *_genomic.fna.gz y *_genomic.gff.gz (y opcionales) usando md5checksums.txt por HTTPS.
    Verifica MD5 si está disponible. Idempotente por defecto.
    Retorna dict con rutas locales: {"fasta": ..., "gff": ..., "protein": ..., "cds": ..., "gbff": ...}
    """
    os.makedirs(output_dir, exist_ok=True)

    # Asegurar HTTPS y preparar índice md5
    base = ftp_url.replace("ftp://", "https://").rstrip("/")
    md5_url = f"{base}/md5checksums.txt"

    md5_text = _http_get(md5_url).decode("utf-8", "ignore")
    md5map = _parse_md5_file_text(md5_text)

    # Derivar nombres esperados a partir del prefijo del directorio
    prefix = os.path.basename(base)  # p.ej., GCF_008452795.1_ASM845279v1
    required_fna = f"{prefix}_genomic.fna.gz"
    required_gff = f"{prefix}_genomic.gff.gz"

    if required_fna not in md5map:
        raise RuntimeError(f"No se encontró {required_fna} en md5checksums.txt")
    if required_gff not in md5map:
        raise RuntimeError(f"No se encontró {required_gff} en md5checksums.txt")

    to_get: Dict[str, str] = {
        "fasta": required_fna,
        "gff":   required_gff,
    }
    optional_candidates = {
        "protein": f"{prefix}_protein.faa.gz",
        "cds":     f"{prefix}_cds_from_genomic.fna.gz",
        "gbff":    f"{prefix}_genomic.gbff.gz",
    }
    for label, fname in optional_candidates.items():
        if label in also_download and fname in md5map:
            to_get[label] = fname  # solo si existe en este ensamblaje

    downloaded: Dict[str, str] = {}
    for label, fname in to_get.items():
        url = f"{base}/{fname}"
        out_path = os.path.join(output_dir, fname)
        _download_file(url, out_path, overwrite=False)
        if check_md5:
            calc = _md5sum(out_path)
            exp = md5map.get(fname, "").lower()
            if exp and calc.lower() != exp:
                raise RuntimeError(f"MD5 no coincide para {fname}: {calc} != {exp}")
        downloaded[label] = out_path

    return downloaded

def _write_provenance(meta: Dict[str, str], out_dir: str) -> None:
    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, "ASSEMBLY_PROVENANCE.txt"), "w", encoding="utf-8") as f:
        for k, v in meta.items():
            f.write(f"{k}: {v}\n")

def fetch_reference_package(
    organism_query: str,
    strain: Optional[str] = None,
    bioproject: Optional[str] = None,
    out_dir: str = "results/genome",
    prefer_refseq: bool = True,
    also_download: Tuple[str, ...] = (),
    check_md5: bool = True,
) -> Tuple[Dict[str, str], str, str]:
    """Descubre el mejor ensamblaje y descarga FASTA+GFF (y opcionales).

    Retorna (meta, fasta_path, gff_path)
    """
    meta = find_best_assembly(
        organism_query=organism_query,
        strain=strain,
        bioproject=bioproject,
        prefer_refseq=prefer_refseq,
    )

    downloaded = download_reference_genome(
        ftp_url=meta["ftp"],
        output_dir=out_dir,
        also_download=also_download,
        check_md5=check_md5,
    )

    _write_provenance(meta, out_dir)
    return meta, downloaded["fasta"], downloaded["gff"]

# ──────────────────────────────────────────────────────────────────────────────
# MODO STANDALONE (CLI)
# ──────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    import argparse
    import textwrap

    parser = argparse.ArgumentParser(
        description="Descubre y descarga el genoma de referencia (FASTA+GFF) desde NCBI Assembly",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent(
            """
            Ejemplos:
              python descarga_genoma_referencia.py -q "Klebsiella sp. AQSCr" -s AQSCr -o results/genome
              python descarga_genoma_referencia.py -q "Klebsiella pneumoniae" --bioproject PRJNA341863 \
                     --also protein cds --check-md5
            """
        ),
    )
    parser.add_argument("-q", "--organism-query", required=True, help="Consulta del organismo (p.ej., 'Klebsiella sp. AQSCr')")
    parser.add_argument("-s", "--strain", default=None, help="Cadena/cepa para filtrar (opcional)")
    parser.add_argument("--bioproject", default=None, help="PRJNA... para restringir (opcional)")
    parser.add_argument("-o", "--out-dir", default="results/genome", help="Directorio de salida")
    parser.add_argument("--no-refseq", action="store_true", help="No preferir RefSeq (GCF_)")
    parser.add_argument("--also", nargs="*", default=(), choices=["protein", "cds", "gbff"], help="Descargas adicionales")
    parser.add_argument("--check-md5", action="store_true", help="Verificar checksums MD5 si está disponible")

    args = parser.parse_args()

    meta, fasta, gff = fetch_reference_package(
        organism_query=args.organism_query,
        strain=args.strain,
        bioproject=args.bioproject,
        out_dir=args.out_dir,
        prefer_refseq=not args.no_refseq,
        also_download=tuple(args.also),
        check_md5=args.check_md5,
    )

    print("\n Ensamblaje seleccionado:")
    print(f"  Accession:       {meta['accession']}")
    print(f"  Organism:        {meta['organism']}")
    print(f"  Assembly name:   {meta['assembly_name']}")
    print(f"  Status:          {meta['assembly_status']}")
    print(f"  RefSeq category: {meta['refseq_category']}")
    print(f"  FTP:             {meta['ftp']}")
    print("\nArchivos descargados:")
    print(f"  FASTA: {fasta}")
    print(f"  GFF:   {gff}")
    if args.also:
        print("  Extra: " + ", ".join(args.also))
    print("\nProvenance: ASSEMBLY_PROVENANCE.txt en el directorio de salida")
