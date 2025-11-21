#!/usr/bin/env python

"""
Prepara insumos para análisis de expresión diferencial (DEA):

1. Lee metadatos de SraRunTable.csv (NCBI SRA).
2. Construye una design matrix con columnas: sample, condition.
3. Reordena la matriz de conteos de featureCounts para que las columnas
   sigan el orden de la design matrix.
4. Escribe:
   - data/metadata/design_matrix.csv
   - data/counts/counts_for_DEA.csv

"""

import argparse
import os
from pathlib import Path

import pandas as pd


def normalize_condition(raw: str) -> str:
    """
    Normaliza el texto de la condición a algo corto y reproducible.

    Ajusta esta función según cómo vengan las condiciones en SraRunTable.
    pensando en nuestro caso:
      - "LB pH 8"           → "control"
      - "LB pH 8/Cr(VI)..." → "treated"
    """
    if pd.isna(raw):
        return "unknown"

    s = str(raw).strip().lower()

    # Ejemplo concreto para nuestro caso (ajustar si cambian cadenas):
    if "cr(vi)" in s or "crvi" in s or "11mm" in s:
        return "treated"
    elif "lb" in s:
        return "control"

    # Fallback: versión limpia del texto original
    return s.replace(" ", "_")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Genera design_matrix.csv y counts_for_DEA.csv a partir de "
            "SraRunTable.csv y counts_final.txt."
        )
    )
    parser.add_argument(
        "--sra-table",
        default="data/metadata/SraRunTable.csv",
        help="Ruta a la tabla SRA (SraRunTable.csv).",
    )
    parser.add_argument(
        "--counts",
        default="data/counts/counts_final.txt",
        help="Matriz de conteos generada por featureCounts.",
    )
    parser.add_argument(
        "--condition-column",
        required=True,
        help=(
            "Nombre de la columna en SraRunTable que contiene la condición "
            "(por ejemplo, 'treatment' o 'LibraryName')."
        ),
    )
    parser.add_argument(
        "--out-design",
        default="data/metadata/design_matrix.csv",
        help="Salida para la design matrix (sample,condition).",
    )
    parser.add_argument(
        "--out-counts",
        default="data/counts/counts_for_DEA.csv",
        help="Salida para la matriz de conteos filtrada/ordenada.",
    )

    args = parser.parse_args()

    sra_path = Path(args.sra_table)
    counts_path = Path(args.counts)
    out_design = Path(args.out_design)
    out_counts = Path(args.out_counts)

    # ─────────────────────────────────────────────────────────────
    # 1) Cargar metadatos de SRA
    # ─────────────────────────────────────────────────────────────
    if not sra_path.exists():
        raise SystemExit(f"[ERROR] No se encontró SraRunTable: {sra_path}")

    meta = pd.read_csv(sra_path)

    if "Run" not in meta.columns:
        raise SystemExit("[ERROR] La tabla SRA no tiene una columna 'Run'.")

    if args.condition_column not in meta.columns:
        raise SystemExit(
            f"[ERROR] La columna de condición '{args.condition_column}' no está en SraRunTable.\n"
            f"Columnas disponibles: {list(meta.columns)}"
        )

    # Nos quedamos con Run y la columna de condición
    meta_sub = meta[["Run", args.condition_column]].copy()
    meta_sub = meta_sub.dropna(subset=["Run"])
    meta_sub["Run"] = meta_sub["Run"].astype(str).str.strip()

    # Normalizar condiciones
    meta_sub["condition"] = meta_sub[args.condition_column].apply(normalize_condition)

    # El nombre de muestra será el Run (SRR1299...)
    meta_sub = meta_sub.rename(columns={"Run": "sample"})

    # Orden estable por nombre de muestra
    meta_sub = meta_sub.sort_values("sample")

    design = meta_sub[["sample", "condition"]].copy()

    out_design.parent.mkdir(parents=True, exist_ok=True)
    design.to_csv(out_design, index=False)
    print(f"Design matrix escrita en: {out_design}")

    # ─────────────────────────────────────────────────────────────
    # 2) Cargar counts_final.txt
    # ─────────────────────────────────────────────────────────────
    if not counts_path.exists():
        raise SystemExit(f"[ERROR] No se encontró la matriz de conteos: {counts_path}")

    # counts_final.txt tiene encabezado con '# Program:...' → usamos comment='#'
    counts = pd.read_csv(counts_path, sep="\t", comment="#")

    # Normalizar nombres de columnas: quitar ruta y sufijo .sorted.bam
    new_cols = []
    for c in counts.columns:
        if c.startswith("data/alignments/") and c.endswith(".sorted.bam"):
            new_cols.append(
                c.replace("data/alignments/", "").replace(".sorted.bam", "")
            )
        else:
            new_cols.append(c)
    counts.columns = new_cols

    # ─────────────────────────────────────────────────────────────
    # 3) Reordenar columnas de conteos según design matrix
    # ─────────────────────────────────────────────────────────────
    ordered_samples = design["sample"].tolist()

    # Comprobar que todas las muestras del diseño están en la matriz de conteos
    missing = [s for s in ordered_samples if s not in counts.columns]
    if missing:
        raise SystemExit(
            "[ERROR] Estas muestras del diseño NO están en la matriz de conteos:\n"
            f"  {missing}\n"
            "Revisa que los nombres de 'sample' coincidan con las columnas "
            "de counts (SRR...)."
        )

    # Seleccionamos solo Geneid + columnas en el orden correcto
    cols_out = ["Geneid"] + ordered_samples
    counts_dea = counts[cols_out].copy()

    out_counts.parent.mkdir(parents=True, exist_ok=True)
    counts_dea.to_csv(out_counts, index=False)
    print(f"Matriz de conteos para DEA escrita en: {out_counts}")

    # Resumen rápido
    print("\nResumen:")
    print(f"  Nº genes: {counts_dea.shape[0]}")
    print(f"  Nº muestras: {len(ordered_samples)}")
    print("  Muestras y condiciones:")
    for s, cond in zip(design["sample"], design["condition"]):
        print(f"    {s}\t{cond}")


if __name__ == "__main__":
    main()
