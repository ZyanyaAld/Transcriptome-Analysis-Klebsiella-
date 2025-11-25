#!/usr/bin/env python3
"""
EggNOG-mapper REAL - Análisis funcional completo
"""
import os
import subprocess
import pandas as pd
from Bio import SeqIO

def extract_protein_sequences(gff_file, fasta_file, output_dir="results/eggnog_real"):
    """Extrae secuencias proteicas REALES del genoma"""
    
    os.makedirs(output_dir, exist_ok=True)
    
    print("🧬 EXTRAYENDO SECUENCIAS PROTEICAS REALES")
    print("=" * 50)
    
    # Leer el genoma
    genome = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    
    protein_sequences = []
    
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
                
                # Extraer coordenadas
                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6]
                seq_id = parts[0]
                
                if seq_id in genome:
                    # Extraer secuencia nucleotídica
                    nucleotide_seq = genome[seq_id].seq[start-1:end]
                    
                    if strand == '-':
                        nucleotide_seq = nucleotide_seq.reverse_complement()
                    
                    # Traducir a proteína (simplificado)
                    protein_seq = nucleotide_seq.translate()
                    
                    # Remover stop codon al final
                    if protein_seq.endswith('*'):
                        protein_seq = protein_seq[:-1]
                    
                    locus_tag = gene_info.get('locus_tag', f"gene_{len(protein_sequences)}")
                    
                    protein_sequences.append({
                        'id': locus_tag,
                        'sequence': str(protein_seq),
                        'product': gene_info.get('product', ''),
                        'location': f"{seq_id}:{start}-{end}({strand})"
                    })
    
    print(f"Secuencias proteicas extraídas: {len(protein_sequences)}")
    return protein_sequences

def create_protein_fasta(protein_sequences, output_file):
    """Crea archivo FASTA con secuencias reales"""
    
    with open(output_file, 'w') as f:
        for protein in protein_sequences:
            f.write(f">{protein['id']} {protein['product']}\n")
            f.write(f"{protein['sequence']}\n")
    
    print(f"Archivo FASTA creado: {output_file}")
    return output_file

def run_eggnog_mapper(fasta_file, output_dir="results/eggnog_real"):
    """Ejecuta EggNOG-mapper real"""
    
    output_prefix = os.path.join(output_dir, "eggnog_results")
    
    print("\n🚀 EJECUTANDO EGGNOG-MAPPER REAL")
    print("=" * 50)
    
    cmd = [
        "emapper.py",
        "-i", fasta_file,
        "-o", output_prefix,
        "--output_dir", output_dir,
        "-m", "diamond",  # Más rápido que hmmer
        "--cpu", "4",
        "--decorate_gff", "none"
    ]
    
    print(f"Comando: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        print("✅ EggNOG-mapper ejecutado exitosamente")
        return True
    except subprocess.CalledProcessError as e:
        print(f"❌ Error ejecutando EggNOG-mapper: {e}")
        print(f"Stderr: {e.stderr}")
        return False
    except FileNotFoundError:
        print("❌ EggNOG-mapper no está instalado")
        return False

def parse_eggnog_results(output_dir):
    """Parsea los resultados de EggNOG"""
    
    annotations_file = os.path.join(output_dir, "eggnog_results.emapper.annotations")
    
    if not os.path.exists(annotations_file):
        print("❌ No se encontraron resultados de EggNOG")
        return None
    
    print(f"\n📊 RESULTADOS DE EGGNOG-MAPPER")
    print("=" * 50)
    
    # Leer resultados
    eggnog_df = pd.read_csv(annotations_file, sep='\t', comment='#')
    
    print(f"Genes anotados: {len(eggnog_df)}")
    print(f"Columnas disponibles: {list(eggnog_df.columns)}")
    
    # Estadísticas básicas
    if 'COG_category' in eggnog_df.columns:
        cog_stats = eggnog_df['COG_category'].value_counts()
        print("\n📈 DISTRIBUCIÓN COG:")
        for cog, count in cog_stats.head(10).items():
            print(f"  {cog}: {count} genes")
    
    if 'Preferred_name' in eggnog_df.columns:
        top_functions = eggnog_df['Preferred_name'].value_counts().head(10)
        print("\n🔝 FUNCIONES MÁS COMUNES:")
        for func, count in top_functions.items():
            print(f"  {func}: {count}")
    
    return eggnog_df

def integrate_with_degs(eggnog_df, deg_file, output_dir):
    """Integra resultados EggNOG con DEGs"""
    
    deg_df = pd.read_csv(deg_file)
    significant_degs = deg_df[
        (deg_df['padj'] < 0.05) & 
        (deg_df['padj'].notna()) & 
        (abs(deg_df['log2FoldChange']) > 1.0)
    ]
    
    # Integrar anotaciones
    merged_df = significant_degs.merge(
        eggnog_df,
        left_on='Geneid', 
        right_on='#query',
        how='left'
    )
    
    # Guardar resultados integrados
    merged_df.to_csv(os.path.join(output_dir, "degs_with_eggnog.csv"), index=False)
    
    print(f"\n🎯 DEGs CON ANOTACIONES EGGNOG: {len(merged_df)}")
    print(f"  - Con anotación completa: {merged_df['Preferred_name'].notna().sum()}")
    
    return merged_df

def main():
    """Función principal"""
    
    gff_file = "data/genome/GCF_008452795.1_ASM845279v1_genomic.gff"
    fasta_file = "data/genome/GCF_008452795.1_ASM845279v1_genomic.fna" 
    deg_file = "data/dea/deseq2_results.csv"
    
    output_dir = "results/eggnog_real"
    
    if not all(os.path.exists(f) for f in [gff_file, fasta_file, deg_file]):
        print("❌ Faltan archivos necesarios")
        return
    
    # 1. Extraer secuencias proteicas reales
    protein_sequences = extract_protein_sequences(gff_file, fasta_file, output_dir)
    
    # 2. Crear FASTA para EggNOG
    fasta_output = os.path.join(output_dir, "proteins_real.faa")
    create_protein_fasta(protein_sequences, fasta_output)
    
    # 3. Ejecutar EggNOG-mapper
    success = run_eggnog_mapper(fasta_output, output_dir)
    
    if success:
        # 4. Parsear resultados
        eggnog_df = parse_eggnog_results(output_dir)
        
        if eggnog_df is not None:
            # 5. Integrar con DEGs
            integrated_df = integrate_with_degs(eggnog_df, deg_file, output_dir)
            
            print(f"\n🎉 ANÁLISIS EGGNOG COMPLETADO!")
            print(f"Resultados en: {output_dir}/")
            
            # Resumen final
            up_regulated = integrated_df[integrated_df['log2FoldChange'] > 0]
            down_regulated = integrated_df[integrated_df['log2FoldChange'] < 0]
            
            print(f"\n📊 RESUMEN FINAL:")
            print(f"  - DEGs up-regulados con EggNOG: {len(up_regulated)}")
            print(f"  - DEGs down-regulados con EggNOG: {len(down_regulated)}")
            
            if 'COG_category' in integrated_df.columns:
                print(f"\n🧬 COGs EN DEGs:")
                cog_summary = integrated_df['COG_category'].value_counts()
                for cog, count in cog_summary.head(5).items():
                    print(f"  {cog}: {count}")

if __name__ == "__main__":
    main()
