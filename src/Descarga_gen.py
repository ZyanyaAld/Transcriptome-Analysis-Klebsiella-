from Bio import SeqIO
import os
import glob

def identificar_genoma_referencia(bioproject_dir):
    """
    Analiza las secuencias descargadas y identifica el genoma de referencia
    """
    
    print(f"\n{'='*70}")
    print("IDENTIFICANDO GENOMA DE REFERENCIA")
    print(f"{'='*70}")
    
    secuencias_dir = os.path.join(bioproject_dir, "secuencias")
    
    if not os.path.exists(secuencias_dir):
        print(" No se encontró la carpeta de secuencias")
        return None
    
    # Encontrar todos los archivos FASTA
    fasta_files = glob.glob(os.path.join(secuencias_dir, "*.fasta"))
    
    if not fasta_files:
        print(" No se encontraron archivos FASTA")
        return None
    
    print(f" Analizando {len(fasta_files)} archivos FASTA...")
    
    candidatos_genoma = []
    todas_secuencias = []
    
    for fasta_file in fasta_files:
        print(f"\n Analizando: {os.path.basename(fasta_file)}")
        
        try:
            # Leer todas las secuencias del archivo
            for record in SeqIO.parse(fasta_file, "fasta"):
                seq_info = {
                    'id': record.id,
                    'descripcion': record.description,
                    'longitud': len(record.seq),
                    'archivo': os.path.basename(fasta_file),
                    'secuencia': record
                }
                todas_secuencias.append(seq_info)
                
                # Mostrar información básica
                print(f"   ├─ {record.id}")
                print(f"   ├─ Descripción: {record.description[:80]}...")
                print(f"   └─ Longitud: {len(record.seq):,} bp")
                
                # CRITERIOS PARA IDENTIFICAR GENOMA DE REFERENCIA
                es_candidato = False
                criterios = []
                
                # Criterio 1: Longitud típica de genoma completo (> 1 Mb)
                if len(record.seq) > 1000000:
                    criterios.append(f"Genoma completo ({len(record.seq):,} bp)")
                    es_candidato = True
                
                # Criterio 2: Palabras clave en descripción
                keywords_genoma = [
                    'complete genome', 'whole genome', 'chromosome', 
                    'complete sequence', 'genomic DNA'
                ]
                desc_lower = record.description.lower()
                for keyword in keywords_genoma:
                    if keyword in desc_lower:
                        criterios.append(f"Keyword: '{keyword}'")
                        es_candidato = True
                        break
                
                # Criterio 3: Patrones de accession number
                # Los genomas completos suelen tener accession numbers específicos
                accession = record.id
                if any(prefix in accession for prefix in ['NC_', 'NZ_', 'AC_', 'AE', 'CP']):
                    criterios.append(f"Accession de genoma: {accession}")
                    es_candidato = True
                
                if es_candidato:
                    candidatos_genoma.append({
                        'info': seq_info,
                        'criterios': criterios,
                        'puntuacion': len(record.seq)  # Usar longitud como puntuación
                    })
                    print(f"    POSIBLE GENOMA - Criterios: {', '.join(criterios)}")
        
        except Exception as e:
            print(f"    Error leyendo archivo: {e}")
    
    # ANALIZAR RESULTADOS
    print(f"\n{'='*70}")
    print("RESULTADOS DEL ANÁLISIS")
    print(f"{'='*70}")
    
    print(f"Total de secuencias analizadas: {len(todas_secuencias)}")
    print(f"Candidatos a genoma de referencia: {len(candidatos_genoma)}")
    
    if not candidatos_genoma:
        print("\n No se identificaron candidatos claros para genoma de referencia")
        print("Mostrando las secuencias más largas:")
        
        # Ordenar por longitud
        todas_secuencias.sort(key=lambda x: x['longitud'], reverse=True)
        
        for i, seq in enumerate(todas_secuencias[:10]):  # Top 10
            print(f"{i+1}. {seq['id']} - {seq['longitud']:,} bp")
            print(f"   Descripción: {seq['descripcion'][:100]}...")
            print(f"   Archivo: {seq['archivo']}")
        
        return None
    
    # Ordenar candidatos por puntuación (longitud)
    candidatos_genoma.sort(key=lambda x: x['puntuacion'], reverse=True)
    
    print(f"\n  CANDIDATOS A GENOMA DE REFERENCIA (ordenados por probabilidad):")
    
    for i, candidato in enumerate(candidatos_genoma[:5]):  # Top 5
        info = candidato['info']
        print(f"\n{'#'*60}")
        print(f"#{i+1}. {info['id']}")
        print(f"{'#'*60}")
        print(f" Longitud: {info['longitud']:,} bp")
        print(f" Descripción: {info['descripcion']}")
        print(f"Archivo: {info['archivo']}")
        print(f" Criterios: {', '.join(candidato['criterios'])}")
        
        # Calcular GC content
        seq_str = str(info['secuencia'].seq)
        gc_content = (seq_str.count('G') + seq_str.count('C')) / len(seq_str) * 100
        print(f"🧬 GC Content: {gc_content:.2f}%")
    
    # EL MEJOR CANDIDATO
    mejor_candidato = candidatos_genoma[0]
    print(f"\n{'='*70}")
    print(f"  GENOMA DE REFERENCIA IDENTIFICADO:")
    print(f"{'='*70}")
    print(f"Accession: {mejor_candidato['info']['id']}")
    print(f"Longitud: {mejor_candidato['info']['longitud']:,} bp")
    print(f"Descripción: {mejor_candidato['info']['descripcion']}")
    print(f"Archivo: {mejor_candidato['info']['archivo']}")
    
    # Guardar el genoma de referencia por separado
    genoma_ref_dir = os.path.join(bioproject_dir, "genoma_referencia")
    os.makedirs(genoma_ref_dir, exist_ok=True)
    
    genoma_file = os.path.join(genoma_ref_dir, "genoma_referencia.fasta")
    SeqIO.write(mejor_candidato['info']['secuencia'], genoma_file, "fasta")
    
    print(f"\n Genoma de referencia guardado en: {genoma_file}")
    
    return mejor_candidato['info']

def analizar_estructura_genoma(genoma_info):
    """
    Analiza la estructura del genoma identificado
    """
    print(f"\n{'='*70}")
    print("ANÁLISIS DE ESTRUCTURA DEL GENOMA")
    print(f"{'='*70}")
    
    secuencia = genoma_info['secuencia']
    
    print(f" Estadísticas del genoma:")
    print(f"   • Longitud total: {len(secuencia.seq):,} bp")
    print(f"   • Longitud de la secuencia: {len(secuencia.seq):,} bp")
    
    # Calcular contenido GC
    seq_str = str(secuencia.seq).upper()
    a_count = seq_str.count('A')
    t_count = seq_str.count('T')
    g_count = seq_str.count('G')
    c_count = seq_str.count('C')
    n_count = seq_str.count('N')  # Bases desconocidas
    
    total_bases = len(seq_str)
    
    print(f"   • Contenido GC: {((g_count + c_count) / total_bases * 100):.2f}%")
    print(f"   • Contenido AT: {((a_count + t_count) / total_bases * 100):.2f}%")
    print(f"   • Bases indeterminadas (N): {n_count} ({n_count/total_bases*100:.4f}%)")
    
    # Buscar características en la descripción
    desc_lower = genoma_info['descripcion'].lower()
    
    if 'plasmid' in desc_lower:
        print("   • Tipo: Plásmido")
    elif 'chromosome' in desc_lower:
        print("   • Tipo: Cromosoma")
    elif 'complete genome' in desc_lower:
        print("   • Tipo: Genoma completo")
    else:
        print("   • Tipo: No especificado")
    
    # Verificar si es circular o lineal
    if 'circular' in desc_lower:
        print("   • Topología: Circular")
    elif 'linear' in desc_lower:
        print("   • Topología: Lineal")
    else:
        print("   • Topología: No especificada")

# EJECUTAR EL ANÁLISIS
if __name__ == "__main__":
    bioproject_id = "PRJNA341863"
    bioproject_dir = os.path.join("results", f"bioproject_{bioproject_id}")
    
    if os.path.exists(bioproject_dir):
        print(f" Analizando BioProject: {bioproject_id}")
        genoma_referencia = identificar_genoma_referencia(bioproject_dir)
        
        if genoma_referencia:
            analizar_estructura_genoma(genoma_referencia)
            
            print(f"\n{'='*70}")
            print(" PRÓXIMOS PASOS RECOMENDADOS:")
            print(f"{'='*70}")
            print("1. Usa el genoma identificado como referencia para alineamientos")
            print("2. Para análisis de RNA-seq, necesitarás descargar los datos SRA por separado")
            print("3. Las proteínas descargadas pueden usarse para anotación funcional")
            print(f"4. El genoma está guardado en: results/bioproject_{bioproject_id}/genoma_referencia/")
    else:
        print(f"No se encontró el directorio: {bioproject_dir}")
        print("Ejecuta primero el script de descarga del BioProject")