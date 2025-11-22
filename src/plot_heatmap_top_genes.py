#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Genera un heatmap de los genes más significativos (por padj)
usando conteos normalizados.

Requiere:
- data/dea/deseq2_results.csv
- data/dea/normalized_counts.csv
- data/metadata/design_matrix.csv

Salida:
- figures/heatmap_top50_genes.png
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

RES_PATH = Path("data/dea/deseq2_results.csv")
NORM_PATH = Path("data/dea/normalized_counts.csv")
DESIGN_PATH = Path("data/metadata/design_matrix.csv")
OUT_DIR = Path("figures")

N_TOP = 50  # número de genes más significativos para incluir

def main():
    # 1) Cargar resultados de DE
    if not RES_PATH.exists():
        raise SystemExit(f"No se encuentra {RES_PATH}")
    de = pd.read_csv(RES_PATH)
    if "Geneid" in de.columns:
        de = de.set_index("Geneid")

    # Quitar NA en padj
    de = de.dropna(subset=["padj"])

    # Ordenar por padj y quedarnos con los N más significativos
    de_sorted = de.sort_values("padj")
    top_genes = de_sorted.head(N_TOP).index.tolist()
    print(f"Usando top {len(top_genes)} genes más significativos.")

    # 2) Cargar conteos normalizados
    if not NORM_PATH.exists():
        raise SystemExit(f"No se encuentra {NORM_PATH}")
    norm = pd.read_csv(NORM_PATH, index_col=0)  # filas=genes, columnas=muestras

    # Asegurar que todos los top_genes estén en la matriz de conteos
    missing_genes = [g for g in top_genes if g not in norm.index]
    if missing_genes:
        print(f"Advertencia: {len(missing_genes)} genes top no están en normalized_counts.")
        top_genes = [g for g in top_genes if g in norm.index]

    # Extraer submatriz
    sub = norm.loc[top_genes]

    # 3) Cargar design_matrix para ordenar columnas y colorear condiciones
    if not DESIGN_PATH.exists():
        raise SystemExit(f"No se encuentra {DESIGN_PATH}")
    design = pd.read_csv(DESIGN_PATH)
    design = design.set_index("sample")

    # Ordenar columnas por design_matrix
    samples_order = design.index.tolist()
    common_samples = [s for s in samples_order if s in sub.columns]
    sub = sub[common_samples]

    # 4) Transformación log2 y z-score por gen (opcional pero recomendable)
    log_sub = np.log2(sub + 1)

    # z-score por gen (fila)
    zscore = (log_sub - log_sub.mean(axis=1).values.reshape(-1, 1)) / (
        log_sub.std(axis=1).values.reshape(-1, 1) + 1e-9
    )

    # 5) Preparar colores para condiciones
    cond = design.loc[common_samples, "condition"]
    cond_unique = cond.unique()

    # Asignar un color simple por condición
    palette = {
        cond_unique[0]: "#1f77b4",  # azul
        cond_unique[1]: "#ff7f0e",  # naranja
    } if len(cond_unique) == 2 else None

    col_colors = cond.map(palette) if palette is not None else None

    # 6) Graficar heatmap
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_file = OUT_DIR / "heatmap_top50_genes.png"

    plt.figure(figsize=(8, 10))
    sns.clustermap(
        zscore,
        cmap="vlag",
        row_cluster=True,
        col_cluster=True,
        col_colors=col_colors,
        xticklabels=common_samples,
        yticklabels=top_genes,
        figsize=(8, 10)
    )

    # Ajustar leyenda de condiciones
    if col_colors is not None:
        # Construir leyenda manual
        for cond_name, color in palette.items():
            plt.scatter([], [], color=color, label=cond_name)
        plt.legend(title="Condición", bbox_to_anchor=(1.05, 1), loc="upper left")

    plt.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Heatmap guardado en: {out_file}")

if __name__ == "__main__":
    main()
