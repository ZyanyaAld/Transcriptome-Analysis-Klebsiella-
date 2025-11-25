#!/usr/bin/env python3
"""
Diagnóstico de matching entre DEGs y anotaciones KEGG
"""
import pandas as pd
import os

def check_gene_matching():
    """Verifica por qué no hay matching entre DEGs y anotaciones KEGG"""
    
    # Cargar DEGs
    deg_df = pd.read_csv("data/dea/deseq2_results.csv")
    significant_degs = deg_df[
        (deg_df['padj'] < 0.05) & 
        (deg_df['padj'].notna()) &
        (abs(deg_df['log2FoldChange']) > 1.0)
    ]
    
    # Cargar anotaciones KEGG
    kegg_df = pd.read_csv("results/kegg/kegg_annotations.csv")
    
    print("🔍 DIAGNÓSTICO DE MATCHING DE GENES")
    print("=" * 50)
    
    print(f"DEGs significativos: {len(significant_degs)}")
    print(f"Anotaciones KEGG: {len(kegg_df)}")
    
    # Mostrar algunos ejemplos de IDs
    print("\n📋 EJEMPLOS DE IDs EN DEGs:")
    print(significant_degs['Geneid'].head(10).tolist())
    
    print("\n📋 EJEMPLOS DE IDs EN ANOTACIONES KEGG:")
    print(kegg_df['gene_id'].head(10).tolist())
    
    # Verificar matching directo
    deg_ids = set(significant_degs['Geneid'])
    kegg_ids = set(kegg_df['gene_id'])
    
    matches = deg_ids.intersection(kegg_ids)
    
    print(f"\n🎯 MATCHING DIRECTO:")
    print(f"Genes que coinciden: {len(matches)}")
    
    if len(matches) > 0:
        print("Ejemplos de matches:")
        for match in list(matches)[:5]:
            print(f"  - {match}")
    else:
        print("❌ NO HAY COINCIDENCIAS DIRECTAS")
        
        # Verificar si hay patrones comunes
        deg_sample = list(deg_ids)[0] if deg_ids else "N/A"
        kegg_sample = list(kegg_ids)[0] if kegg_ids else "N/A"
        
        print(f"\n🔍 PATRONES DE IDs:")
        print(f"Ejemplo DEG: {deg_sample}")
        print(f"Ejemplo KEGG: {kegg_sample}")
        
        # Verificar si son los mismos pero con diferente formato
        if deg_sample and kegg_sample:
            print(f"¿Coinciden?: {deg_sample == kegg_sample}")

if __name__ == "__main__":
    check_gene_matching()

