#!/usr/bin/env python3
"""
Análisis PCA de los conteos normalizados
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import os

def perform_pca_analysis(counts_file, design_file, output_dir="figures"):
    """Realiza análisis PCA y genera gráficas"""
    os.makedirs(output_dir, exist_ok=True)
    
    # Cargar datos
    counts_df = pd.read_csv(counts_file, index_col=0)
    design_df = pd.read_csv(design_file)
    
    # Transformar log2
    log_counts = np.log2(counts_df + 1)
    
    # Transponer para PCA (muestras x genes)
    X = log_counts.T
    
    # Estandarizar
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    # PCA
    pca = PCA(n_components=3)
    principal_components = pca.fit_transform(X_scaled)
    
    # Crear DataFrame de resultados
    pca_df = pd.DataFrame(
        data=principal_components,
        columns=['PC1', 'PC2', 'PC3'],
        index=X.index
    )
    
    # Unir con metadata
    pca_df = pca_df.merge(design_df.set_index('sample'), left_index=True, right_index=True)
    
    # Graficar
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # PC1 vs PC2
    sns.scatterplot(
        data=pca_df,
        x='PC1',
        y='PC2',
        hue='condition',
        s=100,
        ax=ax1
    )
    ax1.set_title(f'PCA - PC1 vs PC2\n(Varianza: {pca.explained_variance_ratio_[0]:.2%})')
    
    # PC1 vs PC3
    sns.scatterplot(
        data=pca_df,
        x='PC1',
        y='PC3',
        hue='condition',
        s=100,
        ax=ax2
    )
    ax2.set_title(f'PCA - PC1 vs PC3\n(Varianza: {pca.explained_variance_ratio_[2]:.2%})')
    
    # Añadir etiquetas de muestras
    for i, sample in enumerate(pca_df.index):
        ax1.annotate(sample, (pca_df.iloc[i]['PC1'], pca_df.iloc[i]['PC2']), 
                    xytext=(5, 5), textcoords='offset points', fontsize=8)
        ax2.annotate(sample, (pca_df.iloc[i]['PC1'], pca_df.iloc[i]['PC3']), 
                    xytext=(5, 5), textcoords='offset points', fontsize=8)
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/pca_analysis.png", dpi=300, bbox_inches='tight')
    plt.show()
    
    # Guardar resultados
    pca_df.to_csv(f"{output_dir}/pca_results.csv")
    
    print(f"Varianza explicada: {pca.explained_variance_ratio_}")
    print(f"Varianza total: {sum(pca.explained_variance_ratio_):.2%}")
    
    return pca_df

if __name__ == "__main__":
    perform_pca_analysis(
        "data/dea/normalized_counts.csv",
        "data/metadata/design_matrix.csv"
    )

