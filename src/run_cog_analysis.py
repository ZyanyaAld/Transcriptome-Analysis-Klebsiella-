#!/usr/bin/env python3
"""
Análisis de categorías COG desde anotaciones existentes
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

def assign_cog_categories(functional_analysis_file, output_dir="results/cog"):
    """Asigna categorías COG basado en anotaciones funcionales"""
    
    os.makedirs(output_dir, exist_ok=True)
    
    print("🔬 ANÁLISIS DE CATEGORÍAS COG")
    print("=" * 50)
    
    # Cargar análisis funcional
    functional_df = pd.read_csv(functional_analysis_file)
    
    # Mapeo de categorías funcionales a COG
    functional_to_cog = {
        'Stress_Response': 'O',  # Posttranslational modification, protein turnover, chaperones
        'DNA_Repair': 'L',       # Replication, recombination and repair
        'Transport': 'G',        # Carbohydrate transport and metabolism
        'Metabolism': 'E',       # Amino acid transport and metabolism
        'Transcription': 'K',    # Transcription
        'Translation': 'J',      # Translation, ribosomal structure and biogenesis
        'Membrane': 'M',         # Cell wall/membrane/envelope biogenesis
        'Energy': 'C',           # Energy production and conversion
        'Unknown': 'S'           # Function unknown
    }
    
    # Asignar categorías COG
    functional_df['COG_category'] = functional_df['functional_category'].map(functional_to_cog)
    
    # Descripciones COG
    cog_descriptions = {
        'J': 'Translation, ribosomal structure and biogenesis',
        'K': 'Transcription',
        'L': 'Replication, recombination and repair',
        'D': 'Cell cycle control, cell division, chromosome partitioning',
        'V': 'Defense mechanisms',
        'T': 'Signal transduction mechanisms',
        'M': 'Cell wall/membrane/envelope biogenesis',
        'N': 'Cell motility',
        'U': 'Intracellular trafficking, secretion, and vesicular transport',
        'O': 'Posttranslational modification, protein turnover, chaperones',
        'C': 'Energy production and conversion',
        'G': 'Carbohydrate transport and metabolism',
        'E': 'Amino acid transport and metabolism',
        'F': 'Nucleotide transport and metabolism',
        'H': 'Coenzyme transport and metabolism',
        'I': 'Lipid transport and metabolism',
        'P': 'Inorganic ion transport and metabolism',
        'Q': 'Secondary metabolites biosynthesis, transport and catabolism',
        'R': 'General function prediction only',
        'S': 'Function unknown'
    }
    
    functional_df['COG_description'] = functional_df['COG_category'].map(cog_descriptions)
    
    # Guardar resultados
    functional_df.to_csv(f"{output_dir}/cog_annotated_genes.csv", index=False)
    
    # Análisis por categoría COG
    cog_analysis = functional_df.groupby(['COG_category', 'COG_description']).agg({
        'gene_id': 'count',
        'log2FoldChange': ['mean', 'std']
    }).round(3)
    
    cog_analysis.columns = ['gene_count', 'mean_log2fc', 'std_log2fc']
    cog_analysis = cog_analysis.sort_values('gene_count', ascending=False)
    cog_analysis.to_csv(f"{output_dir}/cog_analysis_summary.csv")
    
    print("🔝 DISTRIBUCIÓN POR CATEGORÍA COG:")
    for (cog, desc), row in cog_analysis.iterrows():
        up_count = len(functional_df[(functional_df['COG_category'] == cog) & 
                                   (functional_df['log2FoldChange'] > 0)])
        down_count = len(functional_df[(functional_df['COG_category'] == cog) & 
                                     (functional_df['log2FoldChange'] < 0)])
        
        print(f"  {cog} - {desc}: {row['gene_count']} genes ({up_count}↑, {down_count}↓)")
    
    # Generar gráfica COG
    plot_cog_analysis(cog_analysis, output_dir)
    
    return functional_df, cog_analysis

def plot_cog_analysis(cog_analysis, output_dir):
    """Genera gráfica del análisis COG"""
    
    # Preparar datos para gráfica
    plot_data = cog_analysis.reset_index()
    
    plt.figure(figsize=(12, 8))
    
    # Colores por categoría
    colors = plt.cm.Set3(np.linspace(0, 1, len(plot_data)))
    
    bars = plt.bar(range(len(plot_data)), plot_data['gene_count'], color=colors, alpha=0.7)
    
    plt.xlabel('Categoría COG')
    plt.ylabel('Número de Genes')
    plt.title('Distribución de DEGs por Categoría COG\nKlebsiella sp. AqSCr - Respuesta a Cr(VI)')
    plt.xticks(range(len(plot_data)), plot_data['COG_category'], rotation=45)
    
    # Añadir valores en las barras
    for i, bar in enumerate(bars):
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                f'{int(height)}', ha='center', va='bottom')
    
    # Añadir leyenda
    legend_labels = [f"{row['COG_category']} - {row['COG_description'][:30]}..." 
                    for _, row in plot_data.iterrows()]
    plt.legend(bars, legend_labels, bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/cog_analysis_plot.png", dpi=300, bbox_inches='tight')
    plt.show()

def main():
    """Función principal"""
    
    functional_file = "results/functional/functional_categories_detailed.csv"
    
    if not os.path.exists(functional_file):
        print(f"❌ Primero ejecuta el análisis funcional: python src/run_eggnog_from_proteins.py")
        return
    
    # Ejecutar análisis COG
    cog_annotated, cog_analysis = assign_cog_categories(functional_file)
    
    print(f"\n✅ ANÁLISIS COG COMPLETADO!")
    print(f"Resultados en: results/cog/")

if __name__ == "__main__":
    main()
