#!/usr/bin/env python3
"""
PROKKA - Versión mejorada con manejo de errores
"""
import os
import subprocess
import pandas as pd
import shutil
import sys

def check_prokka_installation():
    """Verifica que PROKKA esté instalado correctamente"""
    print("🔍 VERIFICANDO INSTALACIÓN DE PROKKA...")
    
    try:
        # Verificar comando prokka
        result = subprocess.run(["prokka", "--version"], 
                              capture_output=True, text=True, check=True)
        print(f"✅ PROKKA encontrado: {result.stdout.strip()}")
        return True
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        print(f"❌ PROKKA no está instalado o no funciona:")
        print(f"   Error: {e}")
        print("\n💡 SOLUCIÓN: Instala PROKKA con:")
        print("   conda install -c bioconda prokka")
        return False

def setup_environment():
    """Configura el environment para pandas"""
    try:
        import pandas as pd
        import matplotlib.pyplot as plt
        print("✅ Pandas y matplotlib disponibles")
        return True
    except ImportError as e:
        print(f"❌ Error de importación: {e}")
        print("\n💡 SOLUCIÓN: Instala las dependencias:")
        print("   conda install pandas matplotlib seaborn biopython")
        return False

def run_prokka_safe(genome_file, output_dir="results/prokka"):
    """Ejecuta PROKKA de forma segura con manejo de errores"""
    
    # Limpiar directorio existente
    if os.path.exists(output_dir):
        print(f"🗑️  Limpiando directorio existente: {output_dir}")
        shutil.rmtree(output_dir)
    
    os.makedirs(output_dir, exist_ok=True)
    
    print("🧬 EJECUTANDO PROKKA - ANOTACIÓN GENÓMICA")
    print("=" * 50)
    print(f"Genoma: {genome_file}")
    print(f"Salida: {output_dir}")
    
    # Comando PROKKA simplificado y robusto
    cmd = [
        "prokka",
        "--outdir", output_dir,
        "--prefix", "klebsiella", 
        "--cpus", "4",
        "--force",
        "--compliant",  # Hacer salida compliant con NCBI
        genome_file
    ]
    
    print(f"Comando: {' '.join(cmd)}")
    print("⏳ Esto puede tomar 10-30 minutos...")
    
    try:
        # Ejecutar PROKKA con timeout (2 horas)
        result = subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=7200)
        print("✅ PROKKA ejecutado exitosamente!")
        
        # Mostrar output de PROKKA
        if result.stdout:
            print("📋 OUTPUT DE PROKKA:")
            for line in result.stdout.split('\n')[-20:]:  # Últimas 20 líneas
                if line.strip():
                    print(f"   {line}")
        
        return output_dir
        
    except subprocess.TimeoutExpired:
        print("❌ PROKKA excedió el tiempo límite (2 horas)")
        return None
    except subprocess.CalledProcessError as e:
        print(f"❌ Error en PROKKA (código {e.returncode}):")
        if e.stderr:
            print("📋 ERROR DETALLADO:")
            for line in e.stderr.split('\n')[-10:]:  # Últimas 10 líneas de error
                if line.strip():
                    print(f"   {line}")
        return None
    except Exception as e:
        print(f"❌ Error inesperado: {e}")
        return None

def check_prokka_results(prokka_dir):
    """Verifica y procesa los resultados de PROKKA"""
    
    expected_files = [
        "klebsiella.gff",      # Anotaciones
        "klebsiella.faa",      # Proteínas
        "klebsiella.fna",      # Nucleótidos
        "klebsiella.tsv",      # Tabla de anotaciones
    ]
    
    print(f"\n🔍 VERIFICANDO RESULTADOS EN: {prokka_dir}")
    
    missing_files = []
    for file in expected_files:
        file_path = f"{prokka_dir}/{file}"
        if os.path.exists(file_path):
            file_size = os.path.getsize(file_path) / 1024 / 1024  # MB
            print(f"✅ {file}: {file_size:.2f} MB")
        else:
            print(f"❌ {file}: NO ENCONTRADO")
            missing_files.append(file)
    
    if missing_files:
        print(f"⚠️  Faltan {len(missing_files)} archivos importantes")
        return False
    
    return True

def parse_prokka_annotations(prokka_dir):
    """Parsea las anotaciones de PROKKA"""
    
    tsv_file = f"{prokka_dir}/klebsiella.tsv"
    gff_file = f"{prokka_dir}/klebsiella.gff"
    
    print(f"\n📊 PROCESANDO ANOTACIONES PROKKA...")
    
    # Método 1: Usar archivo TSV (más fácil)
    if os.path.exists(tsv_file):
        try:
            annotations_df = pd.read_csv(tsv_file, sep='\t')
            print(f"✅ Anotaciones del TSV: {len(annotations_df)} genes")
            
            # Guardar como CSV
            annotations_df.to_csv(f"{prokka_dir}/prokka_annotations.csv", index=False)
            return annotations_df
            
        except Exception as e:
            print(f"❌ Error leyendo TSV: {e}")
    
    # Método 2: Parsear GFF manualmente
    if os.path.exists(gff_file):
        try:
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
                            'gene': gene_info.get('gene', ''),
                            'protein_id': gene_info.get('protein_id', ''),
                            'contig': parts[0],
                            'start': parts[3],
                            'end': parts[4],
                            'strand': parts[6]
                        })
            
            annotations_df = pd.DataFrame(annotations)
            print(f"✅ Anotaciones del GFF: {len(annotations_df)} genes")
            annotations_df.to_csv(f"{prokka_dir}/prokka_annotations.csv", index=False)
            return annotations_df
            
        except Exception as e:
            print(f"❌ Error parseando GFF: {e}")
    
    return None

def integrate_with_degs(prokka_annotations, deg_file, output_dir="results/prokka"):
    """Integra anotaciones PROKKA con DEGs"""
    
    if prokka_annotations is None:
        print("❌ No hay anotaciones para integrar")
        return None
    
    # Cargar DEGs
    deg_df = pd.read_csv(deg_file)
    significant_degs = deg_df[
        (deg_df['padj'] < 0.05) & 
        (deg_df['padj'].notna()) &
        (abs(deg_df['log2FoldChange']) > 1.0)
    ]
    
    print(f"\n🔗 INTEGRANDO CON {len(significant_degs)} DEGs...")
    
    # Integrar
    degs_with_prokka = significant_degs.merge(
        prokka_annotations,
        left_on='Geneid',
        right_on='locus_tag',
        how='left'
    )
    
    # Guardar resultados
    output_file = f"{output_dir}/degs_with_prokka_annotations.csv"
    degs_with_prokka.to_csv(output_file, index=False)
    
    # Estadísticas
    annotated_count = degs_with_prokka['product'].notna().sum()
    print(f"✅ DEGs con anotaciones PROKKA: {annotated_count}/{len(degs_with_prokka)}")
    
    # Análisis de funciones
    if annotated_count > 0:
        top_functions = degs_with_prokka['product'].value_counts().head(15)
        print(f"\n🔝 TOP 10 FUNCIONES EN DEGs:")
        for func, count in top_functions.head(10).items():
            if pd.notna(func):
                print(f"  {func}: {count} genes")
    
    return degs_with_prokka

def main():
    """Función principal"""
    
    print("🚀 PROKKA - ANOTACIÓN GENÓMICA COMPLETA")
    print("=" * 60)
    
    # 1. Verificar instalación
    if not check_prokka_installation():
        sys.exit(1)
    
    # 2. Verificar dependencias
    if not setup_environment():
        sys.exit(1)
    
    # 3. Archivo de genoma
    genome_file = "data/genome/GCF_008452795.1_ASM845279v1_genomic.fna"
    if not os.path.exists(genome_file):
        print(f"❌ No se encuentra el genoma: {genome_file}")
        sys.exit(1)
    
    # 4. Ejecutar PROKKA
    prokka_dir = run_prokka_safe(genome_file)
    
    if prokka_dir is None:
        print("❌ PROKKA falló. Revisa los errores arriba.")
        sys.exit(1)
    
    # 5. Verificar resultados
    if not check_prokka_results(prokka_dir):
        print("⚠️  Resultados incompletos, pero continuando...")
    
    # 6. Procesar anotaciones
    annotations_df = parse_prokka_annotations(prokka_dir)
    
    # 7. Integrar con DEGs
    deg_file = "data/dea/deseq2_results.csv"
    if os.path.exists(deg_file) and annotations_df is not None:
        integrate_with_degs(annotations_df, deg_file, prokka_dir)
    
    print(f"\n🎯 PROKKA COMPLETADO!")
    print(f"📁 Resultados en: {prokka_dir}")
    print(f"📊 Archivos generados:")
    print(f"   - klebsiella.gff : Anotaciones en formato GFF")
    print(f"   - klebsiella.faa : Secuencias proteicas")
    print(f"   - klebsiella.tsv : Tabla de anotaciones")
    print(f"   - prokka_annotations.csv : Anotaciones procesadas")

if __name__ == "__main__":
    main()
