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
        check_md5=True                       # opcional
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
import urllib.request
from ftplib import FTP
from urllib.parse import urlparse

from Bio import Entrez

# ──────────────────────────────────────────────────────────────────────────────
# CONFIGURACIÓN GLOBAL
# ──────────────────────────────────────────────────────────────────────────────
# Sustituye por un correo real para cumplir políticas de NCBI
Entrez.email = os.environ.get("NCBI_EMAIL", "you@example.com")

# ──────────────────────────────────────────────────────────────────────────────
# UTILIDADES DE DESCARGA
# ──────────────────────────────────────────────────────────────────────────────

def _list_ftp_dir(ftp_url: str) -> List[str]:
    """Lista los archivos en el directorio del ensamblaje usando FTP.

    Acepta ftp:// o https://, internamente convierte a FTP.
    """
    parsed = urlparse(ftp_url.replace("https://", "ftp://"))
    ftp_server, ftp_path = parsed.hostname, parsed.path
    if not ftp_server or not ftp_path:
        raise ValueError(f"URL FTP inválida: {ftp_url}")

    ftp = FTP(ftp_server)
    ftp.login()  # anónimo
    try:
        ftp.cwd(ftp_path)
        files = ftp.nlst()
    finally:
        try:
            ftp.quit()
        except Exception:
            pass
    return files


def _download_file(url: str, out_path: str) -> None:
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    urllib.request.urlretrieve(url, out_path)


def _parse_md5_file(lines: List[str]) -> Dict[str, str]:
    """Parsea el formato habitual de md5checksums.txt (md5  filename)."""
    md5map: Dict[str, str] = {}
    for line in lines:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        # Ejemplo: d41d8cd98f00b204e9800998ecf8427e  GCF_XXXX_genomic.fna.gz
        m = re.match(r"^([0-9a-fA-F]{32})\s+\*?(.+)$", line)
        if m:
            md5, fname = m.group(1).lower(), m.group(2).strip()
            md5map[fname] = md5
    return md5map


def _md5sum(path: str) -> str:
    import hashlib
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
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

    # Buscar IDs
    h = Entrez.esearch(db="assembly", term=q, retmax=max_hits)
    es = Entrez.read(h)
    h.close()
    ids = es.get("IdList", [])

    if not ids:
        # reintentar más laxo
        h = Entrez.esearch(db="assembly", term=organism_query, retmax=max_hits)
        es = Entrez.read(h)
        h.close()
        ids = es.get("IdList", [])
        if not ids:
            raise RuntimeError(f"No se encontraron ensamblajes para: {organism_query}")

    # Resumen completo
    h = Entrez.esummary(db="assembly", id=",".join(ids), report="full")
    recs = Entrez.read(h)
    h.close()
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
        "ftp": ftp.replace("ftp://", "https://"),  # permite descarga vía HTTPS
        "organism": best.get("Organism"),
        "assembly_name": best.get("AssemblyName"),
        "assembly_status": best.get("AssemblyStatus"),
        "refseq_category": best.get("RefSeq_category") or "None",
        "taxid": best.get("Taxid"),
        "bioproject": best.get("BioprojectAccn"),
    }


# ──────────────────────────────────────────────────────────────────────────────
# DESCARGA DE PAQUETE DE REFERENCIA
# ──────────────────────────────────────────────────────────────────────────────

def download_reference_genome(
    ftp_url: str,
    output_dir: str,
    also_download: Tuple[str, ...] = (),  # e.g., ("protein", "cds")
    check_md5: bool = False,
) -> Dict[str, str]:
    """Descarga archivos clave del ensamblaje (FASTA, GFF y opcionales) desde el FTP.

    also_download: tupla con elementos de {"protein", "cds", "gbff"}
    check_md5: si True, intenta verificar checksum con md5checksums.txt

    Retorna un dict con rutas locales descargadas.
    """
    files = _list_ftp_dir(ftp_url)

    required = {
        "fna": next((f for f in files if f.endswith("_genomic.fna.gz")), None),
        "gff": next((f for f in files if f.endswith("_genomic.gff.gz")), None),
    }
    if not required["fna"]:
        raise RuntimeError("No se encontró *_genomic.fna.gz en el directorio del ensamblaje.")
    if not required["gff"]:
        raise RuntimeError("No se encontró *_genomic.gff.gz en el directorio del ensamblaje.")

    optional_map = {
        "protein": next((f for f in files if f.endswith("_protein.faa.gz")), None),
        "cds": next((f for f in files if f.endswith("_cds_from_genomic.fna.gz")), None),
        "gbff": next((f for f in files if f.endswith("_genomic.gbff.gz")), None),
        "md5": next((f for f in files if f.endswith("md5checksums.txt")), None),
    }

    os.makedirs(output_dir, exist_ok=True)

    def build(url_base: str, fname: str) -> Tuple[str, str]:
        return (f"{url_base.rstrip('/')}/{fname}", os.path.join(output_dir, fname))

    downloaded: Dict[str, str] = {}

    # Descargas obligatorias
    fna_url, fna_out = build(ftp_url, required["fna"]) ; _download_file(fna_url, fna_out)
    gff_url, gff_out = build(ftp_url, required["gff"]) ; _download_file(gff_url, gff_out)
    downloaded["fasta"] = fna_out
    downloaded["gff"] = gff_out

    # Descargas opcionales
    if "protein" in also_download and optional_map["protein"]:
        u, o = build(ftp_url, optional_map["protein"]) ; _download_file(u, o)
        downloaded["protein"] = o
    if "cds" in also_download and optional_map["cds"]:
        u, o = build(ftp_url, optional_map["cds"]) ; _download_file(u, o)
        downloaded["cds"] = o
    if "gbff" in also_download and optional_map["gbff"]:
        u, o = build(ftp_url, optional_map["gbff"]) ; _download_file(u, o)
        downloaded["gbff"] = o

    # Verificación de MD5 (opcional)
    if check_md5 and optional_map["md5"]:
        md5_url, md5_out = build(ftp_url, optional_map["md5"]) ; _download_file(md5_url, md5_out)
        with open(md5_out, "r", encoding="utf-8", errors="ignore") as f:
            md5map = _parse_md5_file(f.readlines())
        for label, path in list(downloaded.items()):
            fname = os.path.basename(path)
            if fname in md5map:
                calc = _md5sum(path)
                if calc.lower() != md5map[fname].lower():
                    raise RuntimeError(f"MD5 no coincide para {fname}: {calc} != {md5map[fname]}")

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
    check_md5: bool = False,
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
# MODO STANDALONE
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
