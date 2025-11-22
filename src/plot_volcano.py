#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Genera un volcano plot a partir de los resultados de pydeseq2.

Requiere:
- data/dea/deseq2_results.csv

Salida:
- figures/volcano_plot.png
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

IN_PATH = Path("data/dea/deseq2_results.csv")
OUT_DIR = Path("figures")

def main():
    # 1) Cargar resultados
    if not IN_PATH.exists():
        raise SystemExit(f"No se encuentra {IN_PATH}")

    df = pd.read_csv(IN_PATH)
    # Aseguramos que Geneid sea índice
    if "Geneid" in df.columns:
        df = df.set_index("Geneid")

    # 2) Limpiar valores NA y padj = 0
    df = df.copy()
    df = df.dropna(subset=["padj", "log2FoldChange"])

    # Para evitar -log10(0): reemplazamos padj=0 por el mínimo >0
    nonzero = df["padj"][df["padj"] > 0]
    if len(nonzero) == 0:
        raise SystemExit("Todos los padj son 0 o NA, no se puede hacer volcano plot.")
    min_nonzero = nonzero.min()
    df.loc[df["padj"] <= 0, "padj"] = min_nonzero

    df["neg_log10_padj"] = -np.log10(df["padj"])

    # 3) Definir categorías (ajusta los umbrales si quieres)
    lfc = df["log2FoldChange"]
    padj = df["padj"]

    sig = padj < 0.05
    up = sig & (lfc > 1)
    down = sig & (lfc < -1)
    ns = ~sig

    # 4) Graficar
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_file = OUT_DIR / "volcano_plot.png"

    plt.figure(figsize=(7, 6))

    # No significativos
    plt.scatter(
        df.loc[ns, "log2FoldChange"],
        df.loc[ns, "neg_log10_padj"],
        s=10,
        alpha=0.4,
        label="No sig."
    )

    # Down-regulated
    plt.scatter(
        df.loc[down, "log2FoldChange"],
        df.loc[down, "neg_log10_padj"],
        s=10,
        alpha=0.7,
        label="Down (padj<0.05, log2FC<-1)"
    )

    # Up-regulated
    plt.scatter(
        df.loc[up, "log2FoldChange"],
        df.loc[up, "neg_log10_padj"],
        s=10,
        alpha=0.7,
        label="Up (padj<0.05, log2FC>1)"
    )

    # Líneas de umbral
    plt.axvline(x=1, linestyle="--", linewidth=1)
    plt.axvline(x=-1, linestyle="--", linewidth=1)
    plt.axhline(y=-np.log10(0.05), linestyle="--", linewidth=1)

    plt.xlabel("log2 Fold Change (treated vs control)")
    plt.ylabel("-log10(padj)")
    plt.title("Volcano plot - Klebsiella sp. AqSCr")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_file, dpi=300)
    plt.close()

    print(f"Volcano plot guardado en: {out_file}")

if __name__ == "__main__":
    main()
