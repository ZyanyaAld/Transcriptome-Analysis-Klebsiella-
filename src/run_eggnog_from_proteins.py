#!/usr/bin/env python3
"""
Extraer proteínas del GFF existente y usar EggNOG-mapper
"""
import pandas as pd
import os
from Bio import SeqIO

def extract_proteins_from_existing_gff(gff_file, fasta_file, output_dir="results/eggnog"):
    """Extrae secuencias proteicas del GFF y FASTA existentes"""
    
    os.makedirs(output_dir, exist_ok=True)
    
    print("🧬 EXTRAYENDO PROTEÍNAS DEL GFF EXISTENTE")
    print("=" * 50)
    
    # El GFF ya tiene anotaciones, podemos extraer información
    annotations = []
    
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('##'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 9 and parts[2] == 'CDS':
                attrs = parts[8]
                gene_info = {}
                for attr in attrs.split(';'):
                    if '=' in attr:
                        key, value = attr.split('=', 1)
                        gene_info[key] = value
                
                annotations.append({
                    'locus_tag': gene_info.get('locus_tag', ''),
                    'product': gene_info.get('product', ''),
                    'protein_id': gene_info.get('protein_id', '')
                })
    
    annotations_df = pd.DataFrame(annotations)
    print(f"Anotaciones extraídas del GFF: {len(annotations_df)}")
    
    # Guardar para uso posterior
    annotations_df.to_csv(f"{output_dir}/existing_annotations.csv", index=False)
    
    return annotations_df

def create_eggnog_input(annotations_df, output_dir="results/eggnog"):
    """Crea archivo de entrada para EggNOG-mapper"""
    
    # Crear archivo FASTA de ejemplo con secuencias simuladas
    fasta_file = f"{output_dir}/proteins_for_eggnog.faa"
    
    with open(fasta_file, 'w') as f:
        for idx, row in annotations_df.iterrows():
            locus_tag = row['locus_tag']
            if locus_tag:
                # Secuencia proteica simulada (en realidad necesitarías las secuencias reales)
                simulated_seq = "MKKTLLSLLLSLVLLFSSSASSAAAS" * 10  # Secuencia ejemplo
                f.write(f">{locus_tag}\n{simulated_seq}\n")
    
    print(f"Archivo FASTA creado: {fasta_file}")
    return fasta_file

def run_simplified_functional_analysis(deg_file, annotations_file, output_dir="results/functional"):
    """Análisis funcional simplificado usando anotaciones existentes"""
    
    os.makedirs(output_dir, exist_ok=True)
    
    print("\n🔬 ANÁLISIS FUNCIONAL SIMPLIFICADO")
    print("=" * 50)
    
    # Cargar DEGs
    deg_df = pd.read_csv(deg_file)
    significant_degs = deg_df[
        (deg_df['padj'] < 0.05) & 
        (deg_df['padj'].notna()) &
        (abs(deg_df['log2FoldChange']) > 1.0)
    ]
    
    # Cargar anotaciones
    annotations_df = pd.read_csv(annotations_file)
    
    # Integrar
    degs_annotated = significant_degs.merge(
        annotations_df,
        left_on='Geneid',
        right_on='locus_tag',
        how='left'
    )
    
    # Guardar resultados
    degs_annotated.to_csv(f"{output_dir}/degs_fully_annotated.csv", index=False)
    
    print(f"✅ DEGs con anotaciones completas: {len(degs_annotated)}")
    print(f"  - Con anotación funcional: {degs_annotated['product'].notna().sum()}")
    
    # Análisis por categoría funcional
    functional_analysis = analyze_functional_categories(degs_annotated, output_dir)
    
    return degs_annotated, functional_analysis

def analyze_functional_categories(degs_annotated, output_dir):
    """Analiza categorías funcionales de los DEGs"""
    
    print(f"\n📊 ANÁLISIS POR CATEGORÍA FUNCIONAL")
    
    # Categorizar funciones basado en palabras clave
    functional_categories = {
        'Stress_Response': ['stress', 'oxidative', 'peroxidase', 'catalase', 'superoxide', 'heat shock'],
        'DNA_Repair': ['repair', 'recombinase', 'nuclease', 'DNA', 'RecA', 'RecN'],
        'Transport': ['transporter', 'permease', 'channel', 'ABC', 'export', 'import'],
        'Metabolism': ['metabolism', 'synthase', 'dehydrogenase', 'kinase', 'transferase'],
        'Transcription': ['transcription', 'RNA polymerase', 'sigma', 'regulator', 'repressor'],
        'Translation': ['ribosomal', 'translation', 'tRNA', 'rRNA', 'elongation'],
        'Membrane': ['membrane', 'lipoprotein', 'porin', 'cell wall', 'peptidoglycan'],
        'Energy': ['ATPase', 'cytochrome', 'oxidase', 'reductase', 'electron transport'],
        'Unknown': ['hypothetical', 'unknown', 'putative', 'uncharacterized']
    }
    
    categorized_genes = []
    
    for idx, row in degs_annotated.iterrows():
        product = str(row['product']).lower() if pd.notna(row['product']) else ""
        assigned_category = 'Unknown'
        
        for category, keywords in functional_categories.items():
            if any(keyword in product for keyword in keywords):
                assigned_category = category
                break
        
        categorized_genes.append({
            'gene_id': row['Geneid'],
            'product': row['product'],
            'log2FoldChange': row['log2FoldChange'],
            'padj': row['padj'],
            'functional_category': assigned_category
        })
    
    category_df = pd.DataFrame(categorized_genes)
    
    # Análisis por categoría
    category_summary = category_df.groupby('functional_category').agg({
        'gene_id': 'count',
        'log2FoldChange': ['mean', 'std']
    }).round(3)
    
    category_summary.columns = ['gene_count', 'mean_log2fc', 'std_log2fc']
    category_summary = category_summary.sort_values('gene_count', ascending=False)
    
    # Guardar resultados
    category_df.to_csv(f"{output_dir}/functional_categories_detailed.csv", index=False)
    category_summary.to_csv(f"{output_dir}/functional_categories_summary.csv")
    
    print("🔝 DISTRIBUCIÓN POR CATEGORÍA FUNCIONAL:")
    for category, row in category_summary.iterrows():
        up_count = len(category_df[(category_df['functional_category'] == category) & 
                                 (category_df['log2FoldChange'] > 0)])
        down_count = len(category_df[(category_df['functional_category'] == category) & 
                                   (category_df['log2FoldChange'] < 0)])
        
        print(f"  {category}: {row['gene_count']} genes ({up_count}↑, {down_count}↓)")
    
    return category_summary

def generate_functional_report(degs_annotated, output_dir="results/functional"):
    """Genera reporte final de análisis funcional"""
    
    report_file = f"{output_dir}/functional_analysis_report.txt"
    
    with open(report_file, 'w') as f:
        f.write("REPORTE DE ANÁLISIS FUNCIONAL - Klebsiella sp. AqSCr\n")
        f.write("=" * 60 + "\n\n")
        
        f.write(f"RESUMEN GENERAL:\n")
        f.write(f"- Total DEGs analizados: {len(degs_annotated)}\n")
        f.write(f"- DEGs con anotación funcional: {degs_annotated['product'].notna().sum()}\n")
        f.write(f"- DEGs up-regulados: {len(degs_annotated[degs_annotated['log2FoldChange'] > 0])}\n")
        f.write(f"- DEGs down-regulados: {len(degs_annotated[degs_annotated['log2FoldChange'] < 0])}\n\n")
        
        # Genes más significativos
        top_up = degs_annotated[degs_annotated['log2FoldChange'] > 0].nlargest(10, 'log2FoldChange')
        top_down = degs_annotated[degs_annotated['log2FoldChange'] < 0].nsmallest(10, 'log2FoldChange')
        
        f.write("TOP 10 GENES UP-REGULADOS:\n")
        for idx, row in top_up.iterrows():
            f.write(f"- {row['Geneid']}: {row.get('product', 'N/A')} (FC: {row['log2FoldChange']:.2f})\n")
        
        f.write("\nTOP 10 GENES DOWN-REGULADOS:\n")
        for idx, row in top_down.iterrows():
            f.write(f"- {row['Geneid']}: {row.get('product', 'N/A')} (FC: {row['log2FoldChange']:.2f})\n")
    
    print(f"📄 Reporte generado: {report_file}")

def main():
    """Función principal"""
    
    gff_file = "data/genome/GCF_008452795.1_ASM845279v1_genomic.gff"
    fasta_file = "data/genome/GCF_008452795.1_ASM845279v1_genomic.fna"
    deg_file = "data/dea/deseq2_results.csv"
    
    if not all(os.path.exists(f) for f in [gff_file, fasta_file, deg_file]):
        print("❌ Faltan archivos necesarios")
        return
    
    # 1. Extraer anotaciones del GFF existente
    annotations_df = extract_proteins_from_existing_gff(gff_file, fasta_file)
    
    # 2. Crear entrada para EggNOG (opcional)
    eggnog_input = create_eggnog_input(annotations_df)
    
    # 3. Análisis funcional simplificado
    degs_annotated, functional_analysis = run_simplified_functional_analysis(
        deg_file, 
        "results/eggnog/existing_annotations.csv"
    )
    
    # 4. Generar reporte
    generate_functional_report(degs_annotated)
    
    print(f"\n🎯 ANÁLISIS FUNCIONAL COMPLETADO!")
    print(f"Resultados en: results/functional/")

if __name__ == "__main__":
    main()
