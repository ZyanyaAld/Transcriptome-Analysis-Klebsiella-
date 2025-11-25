#!/usr/bin/env python3
"""
Análisis de enriquecimiento KEGG - VERSIÓN MEJORADA CORREGIDA
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import os

def perform_realistic_kegg_enrichment(deg_kegg_file, output_dir="results/kegg"):
    """Análisis de enriquecimiento más realista"""
    
    os.makedirs(output_dir, exist_ok=True)
    
    print("🔬 ANÁLISIS KEGG MEJORADO - RESPUESTA A Cr(VI)")
    print("=" * 50)
    
    # Cargar datos
    degs_with_kegg = pd.read_csv(deg_kegg_file)
    
    print(f"DEGs analizados: {len(degs_with_kegg)}")
    
    # Análisis por pathway
    pathway_analysis = []
    all_pathways = degs_with_kegg['kegg_pathway'].unique()
    
    for pathway in all_pathways:
        if pd.isna(pathway):
            continue
            
        genes_in_pathway = degs_with_kegg[degs_with_kegg['kegg_pathway'] == pathway]
        total_genes = len(genes_in_pathway)
        
        if total_genes < 5:  # Mínimo para análisis
            continue
            
        up_genes = genes_in_pathway[genes_in_pathway['log2FoldChange'] > 0]
        down_genes = genes_in_pathway[genes_in_pathway['log2FoldChange'] < 0]
        
        # Estadísticas básicas
        up_ratio = len(up_genes) / total_genes
        mean_fc = genes_in_pathway['log2FoldChange'].mean()
        paper_genes = len(genes_in_pathway[genes_in_pathway['paper_reference'] == 'YES'])
        
        pathway_analysis.append({
            'Pathway_ID': pathway,
            'Pathway_Name': get_kegg_pathway_name(pathway),
            'Total_Genes': total_genes,
            'Up_Regulated': len(up_genes),
            'Down_Regulated': len(down_genes),
            'Up_Ratio': up_ratio,
            'Mean_Log2FC': mean_fc,
            'Paper_Genes': paper_genes,
            'Enrichment_Score': up_ratio * abs(mean_fc)  # Score simple
        })
    
    # Crear DataFrame y ordenar por score
    pathway_df = pd.DataFrame(pathway_analysis)
    pathway_df = pathway_df.sort_values('Enrichment_Score', ascending=False)
    
    # Guardar análisis completo
    pathway_df.to_csv(f"{output_dir}/kegg_pathway_analysis_improved.csv", index=False)
    
    # Top pathways
    top_pathways = pathway_df.head(10)
    top_pathways.to_csv(f"{output_dir}/kegg_top_pathways.csv", index=False)
    
    print(f"✅ Pathways analizados: {len(pathway_df)}")
    
    # Mostrar top pathways
    print("\n🔝 TOP PATHWAYS ENRIQUECIDOS:")
    for i, row in pathway_df.head(10).iterrows():
        print(f"  {row['Pathway_Name']}: {row['Total_Genes']} genes ({row['Up_Regulated']}↑, {row['Down_Regulated']}↓), Score: {row['Enrichment_Score']:.2f}")
    
    # Generar gráficas
    generate_kegg_plots(pathway_df, output_dir)
    
    # Análisis específico de respuesta a Cr(VI)
    analyze_cr_response(degs_with_kegg, output_dir)
    
    return pathway_df

def get_kegg_pathway_name(pathway_code):
    """Obtiene nombres de pathways KEGG"""
    pathway_names = {
        'ko02010': 'ABC transporters',
        'ko01200': 'Carbon metabolism', 
        'ko01230': 'Biosynthesis of amino acids',
        'ko00010': 'Glycolysis / Gluconeogenesis',
        'ko00020': 'Citrate cycle (TCA cycle)',
        'ko00190': 'Oxidative phosphorylation',
        'ko03010': 'Ribosome',
        'ko03410': 'Base excision repair',
        'ko02020': 'Two-component system',
        'ko02040': 'Flagellar assembly',
        'ko00520': 'Amino sugar and nucleotide sugar metabolism',
        'ko00680': 'Methane metabolism',
        'ko00910': 'Nitrogen metabolism',
        'ko00920': 'Sulfur metabolism',
        'ko01060': 'Biosynthesis of unsaturated fatty acids'
    }
    return pathway_names.get(pathway_code, pathway_code)

def generate_kegg_plots(pathway_df, output_dir):
    """Genera gráficas del análisis KEGG"""
    
    if pathway_df.empty:
        print("⚠️  No hay datos para generar gráficas")
        return
    
    # Tomar top 15 pathways
    top_plot = pathway_df.head(15).sort_values('Enrichment_Score', ascending=True)
    
    # Gráfica de barras - Score de enriquecimiento
    plt.figure(figsize=(12, 8))
    
    colors = []
    for score in top_plot['Enrichment_Score']:
        if score > 1.0:
            colors.append('red')
        elif score > 0.5:
            colors.append('orange')
        else:
            colors.append('blue')
    
    y_pos = np.arange(len(top_plot))
    bars = plt.barh(y_pos, top_plot['Enrichment_Score'], color=colors, alpha=0.7)
    
    plt.yticks(y_pos, top_plot['Pathway_Name'])
    plt.xlabel('Score de Enriquecimiento (Up-Ratio × |Log2FC|)')
    plt.title('Pathways KEGG Más Enriquecidos en Respuesta a Cr(VI)')
    plt.axvline(x=0.5, color='gray', linestyle='--', alpha=0.5)
    
    # Añadir anotaciones
    for i, (idx, row) in enumerate(top_plot.iterrows()):
        plt.text(row['Enrichment_Score'] + 0.05, i, 
                f"{row['Up_Regulated']}↑/{row['Down_Regulated']}↓", 
                va='center', fontsize=9)
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/kegg_enrichment_improved.png", dpi=300, bbox_inches='tight')
    plt.show()
    
    # Bubble plot - Up vs Down regulation
    plt.figure(figsize=(10, 8))
    
    scatter = plt.scatter(
        top_plot['Up_Regulated'],
        top_plot['Down_Regulated'],
        s=top_plot['Total_Genes'] * 2,
        c=top_plot['Mean_Log2FC'],
        cmap='RdYlBu_r',
        alpha=0.7
    )
    
    plt.xlabel('Genes Up-Regulados')
    plt.ylabel('Genes Down-Regulados')
    plt.title('Regulación por Pathway KEGG\n(Tamaño = total genes, Color = Log2FC promedio)')
    
    # Línea de igual regulación
    max_val = max(top_plot[['Up_Regulated', 'Down_Regulated']].max())
    plt.plot([0, max_val], [0, max_val], 'gray', linestyle='--', alpha=0.5)
    
    plt.colorbar(scatter, label='Log2FC Promedio')
    plt.grid(True, alpha=0.3)
    
    # Añadir etiquetas
    for i, row in top_plot.iterrows():
        plt.annotate(
            row['Pathway_Name'].split()[-1],
            (row['Up_Regulated'], row['Down_Regulated']),
            xytext=(5, 5),
            textcoords='offset points',
            fontsize=8
        )
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/kegg_regulation_bubble.png", dpi=300, bbox_inches='tight')
    plt.show()

def analyze_cr_response(degs_with_kegg, output_dir):
    """Análisis específico de la respuesta a Cr(VI)"""
    
    print("\n🔬 ANÁLISIS ESPECÍFICO - RESPUESTA A Cr(VI)")
    print("=" * 50)
    
    # Categorías funcionales basadas en el paper - CORREGIDO: sin paréntesis
    functional_categories = {
        'Oxidative_Stress': ['ko00190'],
        'DNA_Repair': ['ko03410'],
        'Signal_Transduction': ['ko02020'],
        'Sulfur_Metabolism': ['ko00920'],
        'Iron_Transport': ['ko02010'],
        'Fatty_Acid_Metabolism': ['ko01060', 'ko00061'],
        'Energy_Metabolism': ['ko00020', 'ko01200']
    }
    
    category_analysis = []
    
    # CORREGIDO: usar .items() sin paréntesis
    for category, pathways in functional_categories.items():
        # Filtrar genes en estos pathways
        category_genes = degs_with_kegg[degs_with_kegg['kegg_pathway'].isin(pathways)]
        
        if len(category_genes) == 0:
            continue
            
        up_genes = category_genes[category_genes['log2FoldChange'] > 0]
        down_genes = category_genes[category_genes['log2FoldChange'] < 0]
        paper_genes = category_genes[category_genes['paper_reference'] == 'YES']
        
        category_analysis.append({
            'Functional_Category': category,
            'Total_Genes': len(category_genes),
            'Up_Regulated': len(up_genes),
            'Down_Regulated': len(down_genes),
            'Paper_Genes': len(paper_genes),
            'Mean_Log2FC': category_genes['log2FoldChange'].mean(),
            'Up_Ratio': len(up_genes) / len(category_genes)
        })
        
        print(f"  {category}:")
        print(f"    - Total genes: {len(category_genes)}")
        print(f"    - Up-regulados: {len(up_genes)} ({len(up_genes)/len(category_genes)*100:.1f}%)")
        print(f"    - Down-regulados: {len(down_genes)}")
        print(f"    - Genes del paper: {len(paper_genes)}")
        if len(paper_genes) > 0:
            paper_names = paper_genes['gene_name'].tolist()
            print(f"    - Nombres: {', '.join(paper_names)}")
    
    if not category_analysis:
        print("  No se encontraron genes en las categorías funcionales definidas")
        return
    
    # Guardar análisis funcional
    func_df = pd.DataFrame(category_analysis)
    func_df.to_csv(f"{output_dir}/functional_categories_analysis.csv", index=False)
    
    # Gráfica de categorías funcionales
    plt.figure(figsize=(10, 6))
    
    categories = func_df['Functional_Category']
    up_ratios = func_df['Up_Ratio'] * 100
    
    colors = ['red' if ratio > 60 else 'orange' if ratio > 40 else 'blue' for ratio in up_ratios]
    
    y_pos = np.arange(len(categories))
    plt.barh(y_pos, up_ratios, color=colors, alpha=0.7)
    
    plt.yticks(y_pos, categories)
    plt.xlabel('Porcentaje de Genes Up-Regulados (%)')
    plt.title('Respuesta Funcional a Cr(VI) - Klebsiella sp. AqSCr')
    plt.axvline(x=50, color='gray', linestyle='--', alpha=0.5)
    
    # Añadir valores
    for i, v in enumerate(up_ratios):
        plt.text(v + 1, i, f'{v:.1f}%', va='center', fontsize=9)
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/functional_response_plot.png", dpi=300, bbox_inches='tight')
    plt.show()

def main():
    """Función principal"""
    
    deg_kegg_file = "results/kegg/degs_with_kegg_improved.csv"
    
    if not os.path.exists(deg_kegg_file):
        print(f"❌ Primero ejecuta: python src/get_kegg_annotations_improved.py")
        return
    
    # Ejecutar análisis mejorado
    pathway_results = perform_realistic_kegg_enrichment(deg_kegg_file)
    
    print(f"\n✅ ANÁLISIS KEGG MEJORADO COMPLETADO!")
    print(f"Resultados guardados en: results/kegg/")

if __name__ == "__main__":
    main()


