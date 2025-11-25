#!/usr/bin/env python3
"""
EggNOG Analysis Simple - Usa anotaciones existentes + análisis funcional mejorado
"""
import os
import pandas as pd
import numpy as np
from Bio import SeqIO

def enhanced_functional_analysis(deg_file, gff_file, output_dir="results/eggnog_enhanced"):
    """Análisis funcional mejorado usando anotaciones existentes"""
    
    os.makedirs(output_dir, exist_ok=True)
    
    print("🧬 ANÁLISIS FUNCIONAL MEJORADO")
    print("=" * 50)
    
    # Cargar DEGs
    deg_df = pd.read_csv(deg_file)
    significant_degs = deg_df[
        (deg_df['padj'] < 0.05) & 
        (deg_df['padj'].notna()) & 
        (abs(deg_df['log2FoldChange']) > 1.0)
    ]
    
    print(f"DEGs significativos: {len(significant_degs)}")
    
    # Extraer anotaciones del GFF
    annotations = extract_detailed_annotations(gff_file)
    
    # Integrar con DEGs
    merged_df = significant_degs.merge(
        annotations,
        left_on='Geneid',
        right_on='locus_tag', 
        how='left'
    )
    
    # Análisis por categoría funcional MEJORADO
    functional_analysis = enhanced_category_analysis(merged_df, output_dir)
    
    # Análisis específico para resistencia a metales
    metal_resistance_analysis = analyze_metal_resistance_genes(merged_df, output_dir)
    
    # Generar reportes
    generate_comprehensive_report(merged_df, functional_analysis, metal_resistance_analysis, output_dir)
    
    return merged_df

def extract_detailed_annotations(gff_file):
    """Extrae anotaciones detalladas del GFF"""
    
    annotations = []
    
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('##'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 9 and parts[2] in ['CDS', 'gene']:
                attrs = parts[8]
                gene_info = {}
                for attr in attrs.split(';'):
                    if '=' in attr:
                        key, value = attr.split('=', 1)
                        gene_info[key] = value
                
                annotation = {
                    'locus_tag': gene_info.get('locus_tag', gene_info.get('ID', '')),
                    'product': gene_info.get('product', ''),
                    'protein_id': gene_info.get('protein_id', ''),
                    'gene': gene_info.get('gene', ''),
                    'db_xref': gene_info.get('db_xref', ''),
                    'note': gene_info.get('note', '')
                }
                
                if annotation['locus_tag']:
                    annotations.append(annotation)
    
    return pd.DataFrame(annotations)

def enhanced_category_analysis(df, output_dir):
    """Análisis de categorías funcionales mejorado"""
    
    print("\n📊 ANÁLISIS POR CATEGORÍA FUNCIONAL MEJORADO")
    
    # Categorías expandidas para resistencia a metales
    functional_categories = {
        'Oxidative_Stress': [
            'superoxide', 'catalase', 'peroxidase', 'peroxiredoxin', 'glutathione',
            'thioredoxin', 'redox', 'oxidase', 'peroxidase', 'redoxin'
        ],
        'DNA_Repair': [
            'recA', 'recN', 'recO', 'recR', 'uvrA', 'uvrB', 'uvrC', 'uvrD',
            'mutS', 'mutL', 'nuclease', 'exonuclease', 'endonuclease', 'ligase'
        ],
        'Metal_Resistance': [
            'chromate', 'chromium', 'metal', 'heavy metal', 'efflux', 'transporter',
            'resistance', 'tolerance', 'chrA', 'chrB', 'czc', 'cad', 'cop', 'nik'
        ],
        'Membrane_Transport': [
            'transporter', 'permease', 'channel', 'ABC', 'export', 'import',
            'porin', 'membrane', 'secretion', 'efflux'
        ],
        'Sulfur_Metabolism': [
            'sulfate', 'sulfur', 'cysteine', 'methionine', 'thiosulfate',
            'sulfotransferase', 'sulfatase', 'sulfonate'
        ],
        'Iron_Homeostasis': [
            'iron', 'ferric', 'ferrous', 'siderophore', 'heme', 'hemoglobin',
            'ferritin', 'feo', 'tonB', 'exbB', 'exbD'
        ],
        'Energy_Metabolism': [
            'ATPase', 'cytochrome', 'oxidase', 'reductase', 'dehydrogenase',
            'synthase', 'kinase', 'phosphatase', 'respiratory'
        ],
        'Fatty_Acid_Metabolism': [
            'fatty acid', 'desaturase', 'lipase', 'phospholipase', 'lipid',
            'acyl', 'fabA', 'fabB', 'fabF', 'des'
        ],
        'Signal_Transduction': [
            'sigma', 'rpoE', 'rpoS', 'fecl', 'two-component', 'regulator',
            'transcription factor', 'response regulator'
        ],
        'Unknown': ['hypothetical', 'unknown', 'putative', 'uncharacterized']
    }
    
    categorized_genes = []
    
    for idx, row in df.iterrows():
        product = str(row.get('product', '')).lower()
        note = str(row.get('note', '')).lower()
        combined_text = product + " " + note
        
        assigned_category = 'Unknown'
        
        for category, keywords in functional_categories.items():
            if any(keyword.lower() in combined_text for keyword in keywords):
                assigned_category = category
                break
        
        categorized_genes.append({
            'gene_id': row['Geneid'],
            'locus_tag': row.get('locus_tag', ''),
            'product': row.get('product', ''),
            'log2FoldChange': row['log2FoldChange'],
            'padj': row['padj'],
            'functional_category': assigned_category
        })
    
    category_df = pd.DataFrame(categorized_genes)
    
    # Análisis estadístico
    category_summary = category_df.groupby('functional_category').agg({
        'gene_id': 'count',
        'log2FoldChange': ['mean', 'std'],
        'padj': 'mean'
    }).round(4)
    
    category_summary.columns = ['gene_count', 'mean_log2fc', 'std_log2fc', 'mean_padj']
    category_summary = category_summary.sort_values('gene_count', ascending=False)
    
    # Guardar resultados
    category_df.to_csv(f"{output_dir}/enhanced_functional_categories.csv", index=False)
    category_summary.to_csv(f"{output_dir}/functional_summary_enhanced.csv")
    
    # Imprimir resumen
    print("🔝 DISTRIBUCIÓN MEJORADA:")
    for category, row in category_summary.iterrows():
        up_count = len(category_df[(category_df['functional_category'] == category) & 
                                 (category_df['log2FoldChange'] > 0)])
        down_count = len(category_df[(category_df['functional_category'] == category) & 
                                   (category_df['log2FoldChange'] < 0)])
        
        print(f"  {category}: {row['gene_count']} genes ({up_count}↑, {down_count}↓)")
    
    return category_summary

def analyze_metal_resistance_genes(df, output_dir):
    """Análisis específico para genes de resistencia a metales"""
    
    print("\n🔬 ANÁLISIS ESPECÍFICO - RESISTENCIA A METALES")
    
    # Genes conocidos de resistencia a Cr(VI) y otros metales
    metal_resistance_genes = {
        'Chromate_Resistance': ['chrA', 'chrB', 'chrC', 'chrF', 'yieF'],
        'Oxidative_Stress_Response': ['sodA', 'sodB', 'sodC', 'katG', 'katE', 'ahpC', 'ahpF'],
        'DNA_Repair': ['recA', 'recN', 'uvrA', 'uvrB', 'uvrC', 'uvrD'],
        'Sulfur_Transport': ['cysA', 'cysB', 'cysC', 'cysD', 'cysE', 'cysH', 'cysI', 'cysJ'],
        'Iron_Transport': ['fepA', 'fepB', 'fepC', 'fepD', 'fepG', 'ent', 'tonB', 'exbB', 'exbD'],
        'General_Stress': ['rpoS', 'rpoE', 'rpoH', 'osmC', 'osmY', 'otsA', 'otsB']
    }
    
    metal_analysis = []
    
    for idx, row in df.iterrows():
        gene_id = row['Geneid'].lower()
        product = str(row.get('product', '')).lower()
        
        for category, gene_list in metal_resistance_genes.items():
            if any(gene.lower() in gene_id or gene.lower() in product for gene in gene_list):
                metal_analysis.append({
                    'gene_id': row['Geneid'],
                    'product': row.get('product', ''),
                    'log2FoldChange': row['log2FoldChange'],
                    'padj': row['padj'],
                    'metal_category': category,
                    'known_function': True
                })
                break
        else:
            # Si no es un gen conocido, verificar por función
            if any(keyword in product for keyword in ['chromate', 'chromium', 'heavy metal', 'resistance']):
                metal_analysis.append({
                    'gene_id': row['Geneid'],
                    'product': row.get('product', ''),
                    'log2FoldChange': row['log2FoldChange'],
                    'padj': row['padj'],
                    'metal_category': 'Putative_Metal_Resistance',
                    'known_function': False
                })
    
    metal_df = pd.DataFrame(metal_analysis)
    
    if len(metal_df) > 0:
        metal_df.to_csv(f"{output_dir}/metal_resistance_genes.csv", index=False)
        
        print("🧪 GENES DE RESISTENCIA A METALES IDENTIFICADOS:")
        metal_summary = metal_df.groupby('metal_category').agg({
            'gene_id': 'count',
            'log2FoldChange': 'mean'
        }).round(3)
        
        for category, row in metal_summary.iterrows():
            direction = "↑" if row['log2FoldChange'] > 0 else "↓"
            print(f"  {category}: {row['gene_id']} genes (FC avg: {row['log2FoldChange']:.2f}{direction})")
    
    return metal_df

def generate_comprehensive_report(df, functional_analysis, metal_analysis, output_dir):
    """Genera reporte comprensivo"""
    
    report_file = f"{output_dir}/comprehensive_analysis_report.txt"
    
    with open(report_file, 'w') as f:
        f.write("REPORTE COMPRENSIVO - ANÁLISIS FUNCIONAL KLEBSIELLA AqSCr\n")
        f.write("=" * 70 + "\n\n")
        
        f.write("RESUMEN EJECUTIVO:\n")
        f.write("-" * 50 + "\n")
        f.write(f"Total DEGs analizados: {len(df)}\n")
        f.write(f"DEGs con anotación funcional: {df['product'].notna().sum()}\n")
        f.write(f"DEGs up-regulados: {len(df[df['log2FoldChange'] > 0])}\n")
        f.write(f"DEGs down-regulados: {len(df[df['log2FoldChange'] < 0])}\n\n")
        
        f.write("DISTRIBUCIÓN POR CATEGORÍA FUNCIONAL:\n")
        f.write("-" * 50 + "\n")
        for category, row in functional_analysis.iterrows():
            up_count = len(df[(df.merge(pd.DataFrame([{'functional_category': category}]), 
                                      on='functional_category', how='inner')['log2FoldChange'] > 0)])
            down_count = len(df[(df.merge(pd.DataFrame([{'functional_category': category}]), 
                                        on='functional_category', how='inner')['log2FoldChange'] < 0)])
            f.write(f"{category}: {int(row['gene_count'])} genes ({up_count}↑, {down_count}↓)\n")
        
        f.write("\nGENES CLAVE - RESISTENCIA A Cr(VI):\n")
        f.write("-" * 50 + "\n")
        
        # Genes más significativos en categorías relevantes
        relevant_categories = ['Oxidative_Stress', 'DNA_Repair', 'Metal_Resistance', 
                              'Sulfur_Metabolism', 'Signal_Transduction']
        
        for category in relevant_categories:
            category_genes = df[df.merge(pd.DataFrame([{'functional_category': category}]), 
                                      on='functional_category', how='inner')]
            if len(category_genes) > 0:
                f.write(f"\n{category}:\n")
                top_genes = category_genes.nlargest(5, 'log2FoldChange')
                for idx, row in top_genes.iterrows():
                    f.write(f"  ↑ {row['Geneid']}: {row.get('product', 'N/A')} (FC: {row['log2FoldChange']:.2f})\n")
        
        f.write("\nRECOMENDACIONES:\n")
        f.write("-" * 50 + "\n")
        f.write("1. Validar genes de estrés oxidativo con qPCR\n")
        f.write("2. Estudiar mecanismos de transporte de sulfato\n")
        f.write("3. Analizar modificaciones de membrana (ácidos grasos)\n")
        f.write("4. Investigar sistemas de dos componentes (rpoE, rpoS)\n")
    
    print(f"📄 Reporte comprehensivo generado: {report_file}")

def main():
    """Función principal"""
    
    deg_file = "data/dea/deseq2_results.csv"
    gff_file = "data/genome/GCF_008452795.1_ASM845279v1_genomic.gff"
    
    if not os.path.exists(deg_file):
        print("❌ No se encuentra el archivo de DEGs")
        return
    
    print("🚀 INICIANDO ANÁLISIS FUNCIONAL MEJORADO")
    print("Este análisis usa las anotaciones existentes en tu GFF")
    print("y realiza un análisis comprehensivo de funciones biológicas\n")
    
    enhanced_functional_analysis(deg_file, gff_file)
    
    print(f"\n🎉 ANÁLISIS COMPLETADO!")
    print("Resultados en: results/eggnog_enhanced/")

if __name__ == "__main__":
    main()
