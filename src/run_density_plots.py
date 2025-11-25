#!/usr/bin/env python3
"""
Density plots (KDE) por muestra y condición
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

def create_density_plots(counts_file, design_file, output_dir="figures"):
    """Crea density plots de expresión"""
    os.makedirs(output_dir, exist_ok=True)
    
    # Cargar datos
    counts_df = pd.read_csv(counts_file, index_col=0)
    design_df = pd.read_csv(design_file)
    
    # Transformar log2
    log_counts = np.log2(counts_df + 1)
    
    # Reorganizar datos
    melted_data = log_counts.reset_index().melt(
        id_vars=['Geneid'],
        var_name='sample',
        value_name='log_count'
    )
    
    # Unir con metadata
    melted_data = melted_data.merge(design_df, on='sample')
    
    # Crear figura
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Density plot por muestra
    for sample in log_counts.columns:
        sns.kdeplot(data=log_counts[sample], label=sample, ax=ax1, alpha=0.7)
    
    ax1.set_title('Distribución de Expresión por Muestra')
    ax1.set_xlabel('log2(Conteos + 1)')
    ax1.set_ylabel('Densidad')
    ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Density plot por condición
    for condition in design_df['condition'].unique():
        condition_samples = design_df[design_df['condition'] == condition]['sample']
        condition_data = melted_data[melted_data['sample'].isin(condition_samples)]
        sns.kdeplot(data=condition_data, x='log_count', label=condition, ax=ax2, alpha=0.7)
    
    ax2.set_title('Distribución de Expresión por Condición')
    ax2.set_xlabel('log2(Conteos + 1)')
    ax2.set_ylabel('Densidad')
    ax2.legend()
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/density_plots.png", dpi=300, bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    create_density_plots(
        "data/dea/normalized_counts.csv",
        "data/metadata/design_matrix.csv"
    )

