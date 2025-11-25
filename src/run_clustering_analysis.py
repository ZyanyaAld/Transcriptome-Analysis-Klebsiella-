#!/usr/bin/env python3
"""
Clustering jerárquico y dendrograma - VERSIÓN CORREGIDA
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage, leaves_list
from scipy.spatial.distance import pdist, squareform
import os

def perform_clustering_analysis(counts_file, design_file, output_dir="figures"):
    """Realiza clustering jerárquico - VERSIÓN CORREGIDA"""
    os.makedirs(output_dir, exist_ok=True)
    
    # Cargar datos
    counts_df = pd.read_csv(counts_file, index_col=0)
    design_df = pd.read_csv(design_file)
    
    print(f"Datos cargados: {counts_df.shape}")
    
    # Transformar log2
    log_counts = np.log2(counts_df + 1)
    
    # Calcular matriz de correlación
    corr_matrix = log_counts.corr()
    print("Matriz de correlación calculada")
    
    # CORRECCIÓN: Calcular matriz de distancia correctamente
    distance_matrix = 1 - corr_matrix
    condensed_dist = squareform(distance_matrix, checks=False)
    
    # Clustering jerárquico con matriz de distancia condensada
    linkage_matrix = linkage(condensed_dist, method='average')
    print("Clustering completado")
    
    # Crear figura
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    # Dendrograma
    dendrogram(linkage_matrix, 
               labels=corr_matrix.columns, 
               ax=ax1,
               leaf_rotation=45)
    ax1.set_title('Clustering Jerárquico de Muestras')
    ax1.set_ylabel('Distancia')
    
    # Heatmap de correlación con clustering
    sns.heatmap(
        corr_matrix,
        cmap='coolwarm',
        center=0.8,
        square=True,
        annot=True,
        fmt='.2f',
        ax=ax2,
        cbar_kws={'label': 'Coeficiente de Correlación'}
    )
    ax2.set_title('Matriz de Correlación entre Muestras')
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/clustering_analysis.png", dpi=300, bbox_inches='tight')
    plt.show()
    
    print(f"Gráfica guardada en: {output_dir}/clustering_analysis.png")
    
    return linkage_matrix

if __name__ == "__main__":
    perform_clustering_analysis(
        "data/dea/normalized_counts.csv", 
        "data/metadata/design_matrix.csv"
    )

