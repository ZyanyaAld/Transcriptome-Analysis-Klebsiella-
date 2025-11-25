#!/usr/bin/env python3
"""
Análisis de enriquecimiento de GO (Gene Ontology)
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import gseapy as gp
import os

def perform_go_enrichment(deg_file, output_dir="results/go"):
    """Realiza enriquecimiento de GO"""
    os.makedirs(output_dir, exist_ok=True)
    
    # Cargar DEGs
    deg_df = pd.read_csv(deg_file)
    deg_genes = deg_df[deg_df['padj'] < 0.05]['Geneid'].tolist()
    
    # Aquí iría la anotación GO de tus genes
    # Por ahora es un placeholder
    print(f"DEGs para GO: {len(deg_genes)}")
    
    # Guardar lista de DEGs
    pd.DataFrame({'gene': deg_genes}).to_csv(f"{output_dir}/deg_list.txt", index=False)
    
    return deg_genes

if __name__ == "__main__":
    perform_go_enrichment("data/dea/deseq2_results.csv")

