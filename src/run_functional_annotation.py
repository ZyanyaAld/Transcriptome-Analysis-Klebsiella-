#!/usr/bin/env python3
"""
Anotación funcional usando el archivo GFF del genoma - SIN gffutils
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import re
import os

def parse_gff_simple(gff_file):
    """Extrae anotaciones funcionales del GFF de forma simple"""
    print(f"Parseando {gff_file}...")
    
    annotations = []
    current_sequence = ""
    
    with open(gff_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            if line.startswith('##sequence-region'):
                current_sequence = line.split()[1]
                continue
            if line.startswith('#'):
                continue
                
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
                
            feature_type = parts[2]
            if feature_type not in ['CDS', 'gene']:
                continue
                
            attributes = parts[8]
            
            # Extraer atributos importantes
            locus_tag = extract_attribute(attributes, 'locus_tag')
            gene_name = extract_attribute(attributes, 'gene')
            product = extract_attribute(attributes, 'product')
            protein_id = extract_attribute(attributes, 'protein_id')
            note = extract_attribute(attributes, 'note')
            
            # Si no hay locus_tag, intentar con ID
            if not locus_tag and 'ID=' in attributes:
                locus_tag = extract_attribute(attributes, 'ID')
            
            annotation = {
                'locus_tag': locus_tag,
                'gene': gene_name,
                'product': product,
                'protein_id': protein_id,
                'note': note,
                'feature_type': feature_type,
                'sequence': current_sequence,
                'start': parts[3],
                'end': parts[4],
                'strand': parts[6]
            }
            
            annotations.append(annotation)
    
    df = pd.DataFrame(annotations)
    print(f"Anotaciones extraídas: {len(df)}")
    return df

def extract_attribute(attributes, key):
    """Extrae un atributo específico del string de atributos GFF"""
    patterns = [
        f"{key}=([^;]+)",
        f"{key} \"([^\"]+)\"",
        f"{key}=([^;\\n]+)"
    ]
    
    for pattern in patterns:
        match = re.search(pattern, attributes)
        if match:
            return match.group(1)
    
    return ""

def annotate_degs(deg_file, gff_file, output_dir="results/annotation"):
    """Anota los DEGs con información funcional"""
    os.makedirs(output_dir, exist_ok=True)
    
    print("Cargando DEGs...")
    deg_df = pd.read_csv(deg_file)
    
    # Filtrar DEGs significativos
    significant_degs = deg_df[
        (deg_df['padj'] < 0.05) & 
        (deg_df['padj'].notna()) &
        (np.abs(deg_df['log2FoldChange']) > 1.0)
    ].copy()
    
    print(f"DEGs significativos: {len(significant_degs)}")
    
    # Parsear anotaciones del GFF
    gff_df = parse_gff_simple(gff_file)
    
    if gff_df.empty:
        print("Error: No se pudieron extraer anotaciones del GFF")
        return None
    
    # Unir anotaciones con DEGs
    print("Uniendo anotaciones...")
    deg_annotated = significant_degs.merge(
        gff_df, 
        left_on='Geneid', 
        right_on='locus_tag', 
        how='left'
    )
    
    # Guardar resultados completos
    deg_annotated.to_csv(f"{output_dir}/degs_annotated.csv", index=False)
    
    # Análisis de funciones más comunes
    print("Analizando funciones...")
    functional_summary = deg_annotated['product'].value_counts().head(20)
    functional_summary.to_csv(f"{output_dir}/top_functions.csv")
    
    # Separar up y down regulated
    up_regulated = deg_annotated[deg_annotated['log2FoldChange'] > 0]
    down_regulated = deg_annotated[deg_annotated['log2FoldChange'] < 0]
    
    up_regulated.to_csv(f"{output_dir}/up_regulated_genes.csv", index=False)
    down_regulated.to_csv(f"{output_dir}/down_regulated_genes.csv", index=False)
    
    # Generar reporte
    generate_annotation_report(deg_annotated, output_dir)
    
    return deg_annotated

def generate_annotation_report(deg_annotated, output_dir):
    """Genera un reporte resumen de la anotación"""
    
    report = {
        'total_degs': len(deg_annotated),
        'degs_with_annotation': deg_annotated['product'].notna().sum(),
        'up_regulated': (deg_annotated['log2FoldChange'] > 0).sum(),
        'down_regulated': (deg_annotated['log2FoldChange'] < 0).sum(),
        'top_functions': deg_annotated['product'].value_counts().head(10).to_dict()
    }
    
    # Guardar reporte
    with open(f"{output_dir}/annotation_report.txt", 'w') as f:
        f.write("REPORTE DE ANOTACIÓN FUNCIONAL\n")
        f.write("=" * 50 + "\n")
        f.write(f"Total DEGs: {report['total_degs']}\n")
        f.write(f"DEGs con anotación: {report['degs_with_annotation']}\n")
        f.write(f"Genes sobre-regulados: {report['up_regulated']}\n")
        f.write(f"Genes sub-regulados: {report['down_regulated']}\n")
        f.write("\nTOP 10 FUNCIONES MÁS COMUNES:\n")
        for func, count in report['top_functions'].items():
            if pd.notna(func) and func != "":
                f.write(f"  {func}: {count} genes\n")
    
    # Gráfica de funciones más comunes
    top_funcs = deg_annotated['product'].value_counts().head(15)
    if not top_funcs.empty:
        plt.figure(figsize=(12, 8))
        top_funcs.plot(kind='barh')
        plt.title('Funciones Más Comunes en DEGs - Klebsiella sp. AQSCr')
        plt.xlabel('Número de Genes')
        plt.tight_layout()
        plt.savefig(f"{output_dir}/top_functions_plot.png", dpi=300, bbox_inches='tight')
        plt.show()
    
    print(f"Reporte generado en: {output_dir}/annotation_report.txt")

def search_cr_related_genes(deg_annotated, output_dir):
    """Busca genes potencialmente relacionados con respuesta a Cr(VI)"""
    
    # Términos de búsqueda para genes de respuesta a metales pesados
    cr_keywords = [
        'chromate', 'chromium', 'Cr(VI)', 'heavy metal', 
        'metal resistance', 'transporter', 'efflux',
        'reductase', 'oxidase', 'peroxidase', 'detoxification',
        'stress', 'oxidative', 'antioxidant', 'glutathione',
        'catalase', 'superoxide', 'dismutase'
    ]
    
    print("Buscando genes relacionados con respuesta a Cr(VI)...")
    
    cr_related = []
    for idx, row in deg_annotated.iterrows():
        product = str(row['product']).lower() if pd.notna(row['product']) else ""
        note = str(row['note']).lower() if pd.notna(row['note']) else ""
        
        for keyword in cr_keywords:
            if (keyword in product or keyword in note):
                cr_related.append({
                    'gene_id': row['Geneid'],
                    'log2FoldChange': row['log2FoldChange'],
                    'padj': row['padj'],
                    'product': row['product'],
                    'keyword_found': keyword,
                    'regulation': 'UP' if row['log2FoldChange'] > 0 else 'DOWN'
                })
                break
    
    if cr_related:
        cr_df = pd.DataFrame(cr_related)
        cr_df.to_csv(f"{output_dir}/cr_related_genes.csv", index=False)
        print(f"Genes potencialmente relacionados con Cr(VI): {len(cr_df)}")
        print(cr_df[['gene_id', 'product', 'regulation', 'log2FoldChange']])
    else:
        print("No se encontraron genes directamente relacionados con Cr(VI)")
    
    return cr_related

if __name__ == "__main__":
    # Configurar rutas
    deg_file = "data/dea/deseq2_results.csv"
    gff_file = "data/genome/GCF_008452795.1_ASM845279v1_genomic.gff"
    
    # Ejecutar anotación
    annotated_degs = annotate_degs(deg_file, gff_file)
    
    if annotated_degs is not None:
        # Buscar genes relacionados con Cr(VI)
        search_cr_related_genes(annotated_degs, "results/annotation")
        
        print("\n✅ ANOTACIÓN FUNCIONAL COMPLETADA!")
        print("Archivos generados:")
        print("  - results/annotation/degs_annotated.csv")
        print("  - results/annotation/up_regulated_genes.csv") 
        print("  - results/annotation/down_regulated_genes.csv")
        print("  - results/annotation/top_functions.csv")
        print("  - results/annotation/annotation_report.txt")
        print("  - results/annotation/top_functions_plot.png")
        print("  - results/annotation/cr_related_genes.csv (si se encontraron)")
