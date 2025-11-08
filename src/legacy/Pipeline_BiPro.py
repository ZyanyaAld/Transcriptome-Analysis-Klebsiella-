from Bio import Entrez
from Bio import SeqIO
import time
import os

# Configuración
Entrez.email = "zyanyava@lcg.unam.mx"
results_dir = "results"
os.makedirs(results_dir, exist_ok=True)


# =============================================================================
# DOWNLOAD REFERENCE GENOME AND RNA-SEQ DATA FROM BIOPROJECT
# =============================================================================

def descargar_datos_bioproject(bioproject_id):
    """Download reference genome and associated RNA-Seq data from a BioProject ID."""

    print("\nDESCARGANDO DATOS DEL BIOPROJECT")
    print("-" * 50)
    print(f"BioProject ID: {bioproject_id}")

    try:
        # Search for BioProject record
        handle = Entrez.esearch(db="bioproject", term=bioproject_id)
        record = Entrez.read(handle)
        handle.close()

        if not record["IdList"]:
            print("  ✗ No BioProject records found.")
            return None, None

        bioproject_uid = record["IdList"][0]
        print(f" Found BioProject UID: {bioproject_uid}")

        # Link to the reference genome
        handle = Entrez.elink(dbfrom="bioproject", db="nuccore", id=bioproject_uid)
        links = Entrez.read(handle)
        handle.close()

        genome_ids = []
        for linkset in links[0]["LinkSetDb"]:
            if linkset["LinkName"] == "bioproject_nuccore":
                genome_ids = [l["Id"] for l in linkset["Link"]]
        print(f" Found {len(genome_ids)} linked genome(s).")

        # Download the first reference genome
        if genome_ids:
            genome_handle = Entrez.efetch(db="nuccore", id=genome_ids[0], rettype="fasta", retmode="text")
            genome_seq = genome_handle.read()
            genome_handle.close()

            ref_file = os.path.join(results_dir, f"reference_genome_{bioproject_id}.fasta")
            with open(ref_file, "w") as f:
                f.write(genome_seq)
            print(f"  Reference genome saved as {ref_file}")
        else:
            print("   No linked genome found for this BioProject.")

        # Link to SRA (RNA-Seq data)
        handle = Entrez.elink(dbfrom="bioproject", db="sra", id=bioproject_uid)
        sra_links = Entrez.read(handle)
        handle.close()

        sra_ids = []
        for linkset in sra_links[0]["LinkSetDb"]:
            if linkset["LinkName"] == "bioproject_sra":
                sra_ids = [l["Id"] for l in linkset["Link"]]

        print(f"   Found {len(sra_ids)} linked RNA-Seq runs.")

        # Save metadata for the first few RNA-Seq datasets
        if sra_ids:
            sra_summary = Entrez.esummary(db="sra", id=",".join(sra_ids[:10]))
            sra_records = Entrez.read(sra_summary)

            sra_file = os.path.join(results_dir, f"sra_metadata_{bioproject_id}.txt")
            with open(sra_file, "w", encoding="utf-8") as f:
                for sra in sra_records:
                    f.write(f"Title: {sra.get('Title', 'N/A')}\n")
                    f.write(f"Study: {sra.get('Study', 'N/A')}\n")
                    f.write(f"Organism: {sra.get('Organism', 'N/A')}\n")
                    f.write(f"Sample: {sra.get('Sample', 'N/A')}\n")
                    f.write(f"Accession: {sra.get('Accession', 'N/A')}\n")
                    f.write("-" * 60 + "\n")
            print(f"  SRA metadata saved as {sra_file}")
        else:
            print("   No RNA-Seq links found.")

        return genome_ids, sra_ids

    except Exception as e:
        print(f"  Error downloading BioProject data: {e}")
        return None, None


# Run the BioProject data download
genome_ids, sra_ids = descargar_datos_bioproject("341863")

print("=" * 70)
print("DESCARGANDO SECUENCIAS ESPECIFICAS DE KLEBSIELLA")
print("=" * 70)

# =============================================================================
# BUSCAR Y DESCARGAR SECUENCIAS REALES DE KLEBSIELLA
# =============================================================================

def buscar_secuencias_klebsiella():
    """Buscar secuencias específicas de Klebsiella relacionadas con nuestro proyecto"""
    
    print("\n1. BUSCANDO SECUENCIAS DE KLEBSIELLA RELEVANTES")
    print("-" * 50)
    
    # Términos de búsqueda MÁS AMPLIOS para obtener resultados
    terminos_busqueda = [
        "Klebsiella pneumoniae chromium resistance",  # Más general
        "Klebsiella sp. chromium", 
        "Klebsiella heavy metal resistance",
        "Klebsiella oxidoreductase",
        "Klebsiella chromate",
        "Klebsiella superoxide dismutase",
        "Klebsiella transporter chrA"
    ]
    
    secuencias_encontradas = {}
    todos_ids = set()  # Para evitar duplicados
    
    for termino in terminos_busqueda:
        print(f"\nBuscando: {termino}")
        
        try:
            # Buscar en nucleotide
            handle = Entrez.esearch(db="nucleotide", term=termino, retmax=10)
            result_nuc = Entrez.read(handle)
            handle.close()
            
            count = int(result_nuc['Count'])
            print(f"   Secuencias de DNA encontradas: {count}")
            
            if count > 0:
                secuencias_encontradas[termino] = {
                    'count': count,
                    'ids': result_nuc['IdList'],
                    'db': 'nucleotide'
                }
                
                # Agregar IDs únicos
                for seq_id in result_nuc['IdList']:
                    todos_ids.add(seq_id)
                
                print(f"   IDs encontrados: {len(result_nuc['IdList'])}")
                if result_nuc['IdList']:
                    print(f"   Ejemplo: {result_nuc['IdList'][0]}")
            
            time.sleep(1)  # Respetar límites de NCBI
            
        except Exception as e:
            print(f"   Error en búsqueda: {e}")
            continue
    
    print(f"\nTotal de IDs únicos encontrados: {len(todos_ids)}")
    return secuencias_encontradas, list(todos_ids)

# Buscar secuencias
secuencias_encontradas, todos_ids = buscar_secuencias_klebsiella()

# =============================================================================
# DESCARGAR SECUENCIAS ESPECÍFICAS DE KLEBSIELLA
# =============================================================================

print("\n\n2. DESCARGANDO SECUENCIAS DE KLEBSIELLA")
print("-" * 50)

# Crear subcarpeta para secuencias
secuencias_dir = os.path.join(results_dir, "secuencias_klebsiella")
os.makedirs(secuencias_dir, exist_ok=True)

def descargar_secuencia_klebsiella(seq_id, descripcion="", formato="gb"):
    """Descargar una secuencia específica de Klebsiella"""
    
    try:
        print(f"   Descargando secuencia {seq_id}...")
        
        if formato == "gb":
            filename = os.path.join(secuencias_dir, f"klebsiella_{seq_id}.gb")
            handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="gb", retmode="text")
            
            # Leer y guardar el contenido crudo primero
            raw_data = handle.read()
            handle.close()
            
            if not raw_data or len(raw_data) < 100:
                print(f"   Datos insuficientes para {seq_id}")
                return None
            
            # Guardar datos crudos
            with open(filename, "w", encoding='utf-8') as f:
                f.write(raw_data)
            
            # Intentar parsear para obtener información
            try:
                handle2 = Entrez.efetch(db="nucleotide", id=seq_id, rettype="gb", retmode="text")
                record = SeqIO.read(handle2, "genbank")
                handle2.close()
                
                print(f"   GUARDADO: {filename}")
                print(f"   Organismo: {record.annotations.get('organism', 'N/A')}")
                print(f"   Longitud: {len(record.seq)} bp")
                print(f"   Descripción: {record.description[:100]}...")
                
                # Verificar si es Klebsiella
                organism = record.annotations.get('organism', '')
                if "Klebsiella" in organism:
                    print("  SECUENCIA DE KLEBSIELLA CONFIRMADA")
                    return record
                else:
                    print(f"  Organismo diferente: {organism}")
                    return record
                    
            except Exception as parse_error:
                print(f"   GUARDADO (pero error al parsear): {filename}")
                print(f"   Error parsing: {parse_error}")
                return None
                
    except Exception as e:
        print(f"   ✗ ERROR descargando {seq_id}: {e}")
        return None

# Descargar secuencias encontradas
secuencias_descargadas = []
if todos_ids:
    print(f"\nDescargando {min(5, len(todos_ids))} secuencias de {len(todos_ids)} encontradas...")
    
    for i, seq_id in enumerate(todos_ids[:5]):  # Máximo 5 para prueba
        print(f"\n--- Descargando secuencia {i+1}/{min(5, len(todos_ids))}: {seq_id} ---")
        record = descargar_secuencia_klebsiella(seq_id, "Secuencia de Klebsiella", "gb")
        
        if record is not None:
            secuencias_descargadas.append(record)
            
            # Guardar también en FASTA si es Klebsiella
            try:
                organism = record.annotations.get('organism', '')
                if "Klebsiella" in organism:
                    fasta_filename = os.path.join(secuencias_dir, f"klebsiella_{seq_id}.fasta")
                    with open(fasta_filename, "w") as f:
                        SeqIO.write(record, f, "fasta")
                    print(f"   También guardado como FASTA: {fasta_filename}")
            except Exception as e:
                print(f"   Error guardando FASTA: {e}")
        
        time.sleep(1)  # Respetar límites de NCBI

else:
    print("\nNo se encontraron secuencias con los términos de búsqueda.")
    print("Usando secuencias de referencia...")
    
    # Secuencias de referencia de Klebsiella (IDs conocidos que SÍ existen)
    secuencias_referencia = [
        ("CP000647", "Klebsiella pneumoniae subsp. pneumoniae MGH 78578"),
        ("NC_012731", "Klebsiella pneumoniae NTUH-K2044"),
        ("NZ_CP008827", "Klebsiella pneumoniae strain KPNIH1"),
    ]
    
    for seq_id, descripcion in secuencias_referencia:
        print(f"\n--- Intentando secuencia de referencia {seq_id} ---")
        print(f"Descripción: {descripcion}")
        
        try:
            record = descargar_secuencia_klebsiella(seq_id, descripcion, "gb")
            if record is not None:
                secuencias_descargadas.append(record)
            time.sleep(1)
        except Exception as e:
            print(f"   No se pudo descargar {seq_id}: {e}")

# =============================================================================
# BUSCAR Y DESCARGAR PROTEÍNAS RELACIONADAS
# =============================================================================

print("\n\n3. BUSCANDO PROTEINAS DE KLEBSIELLA RELACIONADAS CON Cr(VI)")
print("-" * 50)

def buscar_proteinas_klebsiella():
    """Buscar proteínas específicas mencionadas en el artículo"""
    
    proteinas_interes = [
        "sodA Klebsiella",  # Más general
        "sodB Klebsiella",  
        "chrA Klebsiella", 
        "chromate transporter Klebsiella",
        "oxidoreductase Klebsiella",
    ]
    
    proteinas_encontradas = {}
    
    for termino in proteinas_interes:
        print(f"\nBuscando: {termino}")
        
        try:
            handle = Entrez.esearch(db="protein", term=termino, retmax=3)
            result = Entrez.read(handle)
            handle.close()
            
            count = int(result['Count'])
            print(f"   Proteínas encontradas: {count}")
            
            if count > 0:
                proteinas_encontradas[termino] = {
                    'count': count,
                    'ids': result['IdList']
                }
                print(f"   IDs: {result['IdList']}")
            else:
                print(f"   No se encontraron proteínas para: {termino}")
            
            time.sleep(1)
            
        except Exception as e:
            print(f"   Error: {e}")
    
    return proteinas_encontradas

# Buscar proteínas
proteinas_encontradas = buscar_proteinas_klebsiella()

# Descargar algunas proteínas si se encontraron
if proteinas_encontradas:
    print("\nDescargando información de proteínas...")
    
    # Crear archivo de resumen de proteínas
    proteinas_file = os.path.join(results_dir, "proteinas_klebsiella.txt")
    
    with open(proteinas_file, "w", encoding="utf-8") as f:
        f.write("PROTEINAS DE KLEBSIELLA RELACIONADAS CON Cr(VI)\n")
        f.write("=" * 50 + "\n\n")
        
        for termino, info in proteinas_encontradas.items():
            f.write(f"PROTEINA: {termino}\n")
            f.write(f"Encontradas: {info['count']}\n")
            f.write(f"IDs: {', '.join(info['ids'])}\n")
            
            # Descargar detalles de la primera proteína de cada categoría
            if info['ids']:
                prot_id = info['ids'][0]
                try:
                    handle = Entrez.efetch(db="protein", id=prot_id, rettype="gb", retmode="text")
                    prot_data = handle.read()
                    handle.close()
                    
                    # Guardar proteína
                    prot_filename = os.path.join(secuencias_dir, f"proteina_{prot_id}.gb")
                    with open(prot_filename, "w", encoding='utf-8') as pf:
                        pf.write(prot_data)
                    
                    print(f" Proteína guardada: {prot_filename}")
                    f.write(f"Archivo: proteina_{prot_id}.gb\n")
                    
                except Exception as e:
                    print(f" Error descargando proteína {prot_id}: {e}")
                    f.write(f"Error: {e}\n")
            
            f.write("-" * 30 + "\n\n")
    
    print(f"Resumen de proteínas guardado en: {proteinas_file}")
else:
    print("\nNo se encontraron proteínas específicas.")

# =============================================================================
# CREAR ARCHIVO CONSOLIDADO DE SECUENCIAS (CORREGIDO COMPLETAMENTE)
# =============================================================================

print("\n\n4. CREANDO ARCHIVOS CONSOLIDADOS")
print("-" * 50)

# Verificar si hay archivos descargados antes de consolidar
if os.path.exists(secuencias_dir):
    archivos_gb = [f for f in os.listdir(secuencias_dir) if f.endswith(".gb")]
else:
    archivos_gb = []

if archivos_gb:
    print(f"Encontrados {len(archivos_gb)} archivos GenBank para procesar")
    
    # SEPARAR archivos de nucleótidos vs proteínas
    archivos_nucleotidos = [f for f in archivos_gb if f.startswith('klebsiella_') or f.startswith('genoma_') or f.startswith('bioproyecto_') or f.startswith('cepa_')]
    archivos_proteinas = [f for f in archivos_gb if f.startswith('proteina_')]
    
    print(f"  - Secuencias de nucleótidos: {len(archivos_nucleotidos)}")
    print(f"  - Secuencias de proteínas: {len(archivos_proteinas)}")
    
    # =========================================================================
    # PROCESAR SECUENCIAS DE NUCLEÓTIDOS PARA FASTA
    # =========================================================================
    
    if archivos_nucleotidos:
        print(f"\nProcesando {len(archivos_nucleotidos)} secuencias de nucleótidos...")
        
        todas_secuencias_adn = []
        secuencias_con_problemas = []
        
        for archivo in archivos_nucleotidos:
            try:
                filepath = os.path.join(secuencias_dir, archivo)
                record = SeqIO.read(filepath, "genbank")
                
                # VERIFICAR CRÍTICAMENTE que la secuencia esté definida y sea de ADN
                if (record.seq is not None and 
                    len(record.seq) > 0 and 
                    hasattr(record, 'seq') and
                    str(record.seq) != ''):
                    
                    # Verificar que no sea una secuencia de proteína disfrazada
                    if not any(keyword in archivo.lower() for keyword in ['proteina', 'protein']):
                        todas_secuencias_adn.append(record)
                        print(f" Nucleótido: {archivo} ({len(record.seq)} bp)")
                    else:
                        print(f" Omitido (es proteína): {archivo}")
                else:
                    secuencias_con_problemas.append(archivo)
                    print(f"  Secuencia indefinida/vacía: {archivo}")
                    
            except Exception as e:
                secuencias_con_problemas.append(archivo)
                print(f" Error procesando {archivo}: {e}")
        
        # CREAR ARCHIVO FASTA SOLO PARA NUCLEÓTIDOS
        if todas_secuencias_adn:
            fasta_consolidado = os.path.join(results_dir, "klebsiella_secuencias_consolidadas.fasta")
            
            try:
                with open(fasta_consolidado, "w") as f:
                    # Usar write directamente para mayor control
                    for record in todas_secuencias_adn:
                        # Formatear manualmente el FASTA para evitar problemas
                        header = f">{record.id} {record.description}"
                        sequence = str(record.seq)
                        
                        # Escribir header y secuencia en líneas de 80 caracteres
                        f.write(header + "\n")
                        for i in range(0, len(sequence), 80):
                            f.write(sequence[i:i+80] + "\n")
                
                print(f" Archivo consolidado FASTA creado: {fasta_consolidado}")
                print(f"  Contiene {len(todas_secuencias_adn)} secuencias de nucleótidos")
                print(f"  Longitud total: {sum(len(rec.seq) for rec in todas_secuencias_adn):,} bp")
                
            except Exception as e:
                print(f" Error creando archivo FASTA: {e}")
        else:
            print(" No hay secuencias de nucleótidos válidas para el archivo FASTA")
    
    # =========================================================================
    # PROCESAR SECUENCIAS DE PROTEÍNAS POR SEPARADO
    # =========================================================================
    
    if archivos_proteinas:
        print(f"\nProcesando {len(archivos_proteinas)} secuencias de proteínas...")
        
        todas_proteinas = []
        
        for archivo in archivos_proteinas:
            try:
                filepath = os.path.join(secuencias_dir, archivo)
                record = SeqIO.read(filepath, "genbank")
                
                # Para proteínas, verificar que tenga secuencia de aminoácidos
                if (record.seq is not None and 
                    len(record.seq) > 0 and
                    hasattr(record, 'seq')):
                    
                    todas_proteinas.append(record)
                    print(f"  Proteína: {archivo} ({len(record.seq)} aa)")
                    
            except Exception as e:
                print(f"  Error procesando proteína {archivo}: {e}")
        
        # CREAR ARCHIVO FASTA SEPARADO PARA PROTEÍNAS
        if todas_proteinas:
            fasta_proteinas = os.path.join(results_dir, "klebsiella_proteinas_consolidadas.fasta")
            
            try:
                with open(fasta_proteinas, "w") as f:
                    for record in todas_proteinas:
                        header = f">{record.id} {record.description}"
                        sequence = str(record.seq)
                        
                        f.write(header + "\n")
                        for i in range(0, len(sequence), 80):
                            f.write(sequence[i:i+80] + "\n")
                
                print(f"Archivo de proteínas FASTA creado: {fasta_proteinas}")
                print(f"  Contiene {len(todas_proteinas)} secuencias de proteínas")
                
            except Exception as e:
                print(f"Error creando archivo de proteínas: {e}")
    
    # =========================================================================
    # CREAR ARCHIVO DE METADATOS MEJORADO
    # =========================================================================
    
    metadatos_file = os.path.join(results_dir, "metadatos_secuencias.txt")
    
    with open(metadatos_file, "w", encoding="utf-8") as f:
        f.write("METADATOS DE SECUENCIAS DESCARGADAS\n")
        f.write("=" * 50 + "\n\n")
        
        # Procesar todos los archivos para metadatos
        for archivo in archivos_gb:
            try:
                filepath = os.path.join(secuencias_dir, archivo)
                record = SeqIO.read(filepath, "genbank")
                
                f.write(f"ARCHIVO: {archivo}\n")
                f.write(f"ID: {record.id}\n")
                f.write(f"Tipo: {'PROTEÍNA' if archivo.startswith('proteina_') else 'NUCLEÓTIDO'}\n")
                f.write(f"Organismo: {record.annotations.get('organism', 'N/A')}\n")
                
                # Longitud según el tipo
                if record.seq is not None:
                    if archivo.startswith('proteina_'):
                        f.write(f"Longitud: {len(record.seq)} aminoácidos\n")
                    else:
                        f.write(f"Longitud: {len(record.seq)} bp\n")
                else:
                    f.write(f"Longitud: SECUENCIA INDEFINIDA\n")
                
                f.write(f"Fecha: {record.annotations.get('date', 'N/A')}\n")
                f.write(f"Definición: {record.description}\n")
                
                # Información específica por tipo
                if archivo.startswith('proteina_'):
                    # Información de proteínas
                    product = record.annotations.get('product', [''])[0] if isinstance(record.annotations.get('product'), list) else record.annotations.get('product', 'N/A')
                    f.write(f"Producto: {product}\n")
                else:
                    # Información de genes para nucleótidos
                    genes = []
                    for feature in record.features:
                        if feature.type == "CDS" and 'gene' in feature.qualifiers:
                            genes.append(feature.qualifiers['gene'][0])
                    
                    if genes:
                        f.write(f"Genes encontrados: {', '.join(set(genes[:8]))}")
                        if len(genes) > 8:
                            f.write(f"... y {len(genes)-8} más")
                        f.write("\n")
                
                f.write("-" * 40 + "\n\n")
                
            except Exception as e:
                f.write(f"ARCHIVO: {archivo} - ERROR: {e}\n\n")
    
    print(f"Metadatos guardados en: {metadatos_file}")
    
    # RESUMEN FINAL DE PROCESAMIENTO
    print(f"\n" + "=" * 50)
    print("RESUMEN DE PROCESAMIENTO:")
    if 'todas_secuencias_adn' in locals() and todas_secuencias_adn:
        print(f"Secuencias de nucleótidos: {len(todas_secuencias_adn)}")
    if 'todas_proteinas' in locals() and todas_proteinas:
        print(f"Secuencias de proteínas: {len(todas_proteinas)}")
    if 'secuencias_con_problemas' in locals() and secuencias_con_problemas:
        print(f"Archivos con problemas: {len(secuencias_con_problemas)}")
    
else:
    print("No hay archivos GenBank para consolidar")
    print("  La carpeta de secuencias está vacía o no existe")

# =============================================================================
# RESUMEN FINAL MEJORADO
# =============================================================================

print("\n" + "=" * 70)
print("RESUMEN DE DESCARGAS - KLEBSIELLA")
print("=" * 70)

# Mostrar lo que se descargó
if os.path.exists(secuencias_dir):
    archivos = os.listdir(secuencias_dir)
    print(f"\nARCHIVOS EN {secuencias_dir}/:")
    if archivos:
        for archivo in archivos:
            tamaño = os.path.getsize(os.path.join(secuencias_dir, archivo))
            print(f"  - {archivo} ({tamaño} bytes)")
    else:
        print("  (ningún archivo descargado)")
else:
    print(f"\nLA CARPETA {secuencias_dir} NO EXISTE")

print(f"\nRESUMEN:")
print(f"  - Secuencias encontradas en búsqueda: {len(todos_ids)}")
print(f"  - Secuencias descargadas exitosamente: {len(secuencias_descargadas)}")
print(f"  - Proteínas encontradas: {len(proteinas_encontradas)}")

print(f"\nPRÓXIMOS PASOS:")
if secuencias_descargadas:
    print("1. Secuencias descargadas en results/secuencias_klebsiella/")
    print("2. Revisa que sean realmente de Klebsiella")
    print("3. Usa estas secuencias para tus análisis bioinformáticos")
else:
    print("1. No se pudieron descargar secuencias automáticamente")
    print("2. Intenta con IDs específicos manualmente")
    print("3. Revisa la conexión a internet y los términos de búsqueda")

print("\n" + "=" * 70)