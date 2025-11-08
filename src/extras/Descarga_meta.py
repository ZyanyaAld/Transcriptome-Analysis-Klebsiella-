from Bio import Entrez
from Bio import SeqIO
import time
import os
import xml.etree.ElementTree as ET
from urllib.error import HTTPError

# Configuración
Entrez.email = "zyanyava@lcg.unam.mx"
Entrez.api_key = None  # Opcional: agregar API key si tienes
results_dir = "results"
os.makedirs(results_dir, exist_ok=True)

# =============================================================================
# FUNCIÓN PRINCIPAL
# =============================================================================

def descargar_bioproject_completo(bioproject_id):
    """
    Descarga todos los datos disponibles de un BioProject:
    - Metadatos del BioProject
    - Genomas de referencia
    - Secuencias relacionadas (nucleotide)
    - Datos de SRA (RNA-Seq, etc.)
    - Proteínas relacionadas
    """
    
    print(f"\n{'='*70}")
    print(f"DESCARGANDO DATOS COMPLETOS DEL BIOPROJECT: {bioproject_id}")
    print(f"{'='*70}")
    
    # Crear carpeta específica para este BioProject
    bioproject_dir = os.path.join(results_dir, f"bioproject_{bioproject_id}")
    os.makedirs(bioproject_dir, exist_ok=True)
    
    try:
        # =====================================================================
        # 1. OBTENER METADATOS DEL BIOPROJECT
        # =====================================================================
        print("\n1. OBTENIENDO METADATOS DEL BIOPROJECT...")
        
        # Buscar BioProject
        handle = Entrez.esearch(db="bioproject", term=bioproject_id)
        record = Entrez.read(handle)
        handle.close()
        
        if not record["IdList"]:
            print(f"  ✗ No se encontró el BioProject: {bioproject_id}")
            return None
        
        bioproject_uid = record["IdList"][0]
        print(f"  ✓ BioProject UID: {bioproject_uid}")
        
        # Obtener detalles completos del BioProject - CORREGIDO
        try:
            handle = Entrez.efetch(db="bioproject", id=bioproject_uid, retmode="xml")
            bioproject_data = handle.read()
            
            # Decodificar si es bytes
            if isinstance(bioproject_data, bytes):
                bioproject_xml = bioproject_data.decode('utf-8')
            else:
                bioproject_xml = bioproject_data
                
            handle.close()
            
            # Guardar metadatos XML
            metadata_file = os.path.join(bioproject_dir, f"bioproject_{bioproject_id}_metadata.xml")
            with open(metadata_file, "w", encoding="utf-8") as f:
                f.write(bioproject_xml)
            print(f"  ✓ Metadatos guardados: {metadata_file}")
            
            # Intentar parsear información básica
            try:
                root = ET.fromstring(bioproject_xml)
                # Buscar título en diferentes ubicaciones posibles
                title_elem = root.find(".//Project/ProjectDescr/Title")
                if title_elem is None:
                    title_elem = root.find(".//Project/Project/ProjectDescr/Title")
                if title_elem is None:
                    title_elem = root.find(".//Title")
                
                title = title_elem.text if title_elem is not None else "N/A"
                print(f"  ✓ Título del proyecto: {title}")
            except Exception as parse_error:
                print(f"  ℹ No se pudo parsear XML del proyecto: {parse_error}")
                
        except Exception as e:
            print(f"  ✗ Error obteniendo metadatos XML: {e}")
            # Continuar con otros datos aunque falle este paso
        
        # =====================================================================
        # 2. DESCARGAR GENOMAS DE REFERENCIA Y SECUENCIAS RELACIONADAS
        # =====================================================================
        print("\n2. BUSCANDO GENOMAS Y SECUENCIAS RELACIONADAS...")
        
        # Buscar en nucleotide database - término de búsqueda corregido
        search_term = f"{bioproject_id}[BioProject]"
        print(f"  Término de búsqueda: {search_term}")
        
        try:
            handle = Entrez.esearch(db="nuccore", term=search_term, retmax=100)
            nuc_records = Entrez.read(handle)
            handle.close()
            
            nuc_count = int(nuc_records["Count"])
            nuc_ids = nuc_records["IdList"]
            
            print(f"  ✓ Secuencias encontradas en Nucleotide: {nuc_count}")
            
            # Crear subcarpeta para secuencias
            secuencias_dir = os.path.join(bioproject_dir, "secuencias")
            os.makedirs(secuencias_dir, exist_ok=True)
            
            # Descargar secuencias en lotes
            if nuc_ids:
                print(f"  Descargando {len(nuc_ids)} secuencias...")
                descargar_secuencias_lotes(nuc_ids, secuencias_dir, "nuccore")
            else:
                print("  ℹ No se encontraron secuencias de nucleotide")
                
        except Exception as e:
            print(f"  ✗ Error buscando secuencias: {e}")
            nuc_count = 0
            nuc_ids = []
        
        # =====================================================================
        # 3. BUSCAR Y DESCARGAR DATOS SRA (RNA-SEQ, ETC.)
        # =====================================================================
        print("\n3. BUSCANDO DATOS DE EXPRESIÓN (SRA)...")
        
        try:
            handle = Entrez.esearch(db="sra", term=search_term, retmax=50)
            sra_records = Entrez.read(handle)
            handle.close()
            
            sra_count = int(sra_records["Count"])
            sra_ids = sra_records["IdList"]
            
            print(f"  ✓ Experimentos SRA encontrados: {sra_count}")
            
            if sra_ids:
                # Descargar metadatos SRA
                descargar_metadatos_sra(sra_ids, bioproject_dir, bioproject_id)
            else:
                print("  ℹ No se encontraron experimentos SRA")
                
        except Exception as e:
            print(f"  ✗ Error buscando datos SRA: {e}")
            sra_count = 0
            sra_ids = []
        
        # =====================================================================
        # 4. BUSCAR PROTEÍNAS RELACIONADAS
        # =====================================================================
        print("\n4. BUSCANDO PROTEÍNAS RELACIONADAS...")
        
        try:
            handle = Entrez.esearch(db="protein", term=search_term, retmax=50)
            prot_records = Entrez.read(handle)
            handle.close()
            
            prot_count = int(prot_records["Count"])
            prot_ids = prot_records["IdList"]
            
            print(f"  ✓ Proteínas encontradas: {prot_count}")
            
            if prot_ids:
                # Crear subcarpeta para proteínas
                proteinas_dir = os.path.join(bioproject_dir, "proteinas")
                os.makedirs(proteinas_dir, exist_ok=True)
                
                print(f"  Descargando {len(prot_ids)} proteínas...")
                descargar_secuencias_lotes(prot_ids, proteinas_dir, "protein")
            else:
                print("  ℹ No se encontraron proteínas")
                
        except Exception as e:
            print(f"   Error buscando proteínas: {e}")
            prot_count = 0
            prot_ids = []
        
        # =====================================================================
        # 5. BUSCAR EN PUBMED (ARTÍCULOS RELACIONADOS)
        # =====================================================================
        print("\n5. BUSCANDO PUBLICACIONES RELACIONADAS...")
        
        try:
            handle = Entrez.esearch(db="pubmed", term=search_term, retmax=20)
            pubmed_records = Entrez.read(handle)
            handle.close()
            
            pubmed_count = int(pubmed_records["Count"])
            pubmed_ids = pubmed_records["IdList"]
            
            print(f"   Publicaciones encontradas: {pubmed_count}")
            
            if pubmed_ids:
                descargar_metadatos_pubmed(pubmed_ids, bioproject_dir, bioproject_id)
            else:
                print("  ℹ No se encontraron publicaciones")
                
        except Exception as e:
            print(f"   Error buscando publicaciones: {e}")
            pubmed_count = 0
            pubmed_ids = []
        
        # =====================================================================
        # 6. CREAR ARCHIVOS CONSOLIDADOS Y RESUMEN
        # =====================================================================
        print("\n6. CREANDO ARCHIVOS CONSOLIDADOS Y RESUMEN...")
        crear_resumen_bioproject(bioproject_dir, bioproject_id, {
            'nucleotide': nuc_count,
            'sra': sra_count,
            'protein': prot_count,
            'pubmed': pubmed_count
        })
        
        # =====================================================================
        # 7. RESUMEN FINAL
        # =====================================================================
        print(f"\n{'='*70}")
        print(f"DESCARGA COMPLETADA: {bioproject_id}")
        print(f"{'='*70}")
        print(f"✓ Metadatos del BioProject")
        print(f"✓ Secuencias de nucleotide: {nuc_count}")
        print(f"✓ Experimentos SRA: {sra_count}")
        print(f"✓ Secuencias de proteínas: {prot_count}")
        print(f"✓ Publicaciones: {pubmed_count}")
        print(f"\nTodos los datos guardados en: {bioproject_dir}")
        
        return bioproject_dir
        
    except Exception as e:
        print(f"  ✗ ERROR general procesando BioProject: {e}")
        import traceback
        traceback.print_exc()
        return None

# =============================================================================
# FUNCIONES AUXILIARES 
# =============================================================================

def descargar_secuencias_lotes(ids_list, output_dir, db_type):
    """Descargar secuencias en lotes para evitar timeouts"""
    
    batch_size = 5  # Reducido para mayor estabilidad
    formato = "fasta"  # Formato más simple y confiable
    
    if not ids_list:
        print("    ℹ Lista de IDs vacía")
        return
    
    for i in range(0, len(ids_list), batch_size):
        batch_ids = ids_list[i:i + batch_size]
        batch_num = (i // batch_size) + 1
        total_batches = (len(ids_list) + batch_size - 1) // batch_size
        
        print(f"    Procesando lote {batch_num}/{total_batches} ({len(batch_ids)} secuencias)...")
        
        try:
            # Descargar en formato FASTA
            handle = Entrez.efetch(db=db_type, id=",".join(batch_ids), 
                                 rettype=formato, retmode="text")
            batch_data = handle.read()
            
            # Verificar y decodificar si es necesario
            if isinstance(batch_data, bytes):
                batch_data = batch_data.decode('utf-8')
                
            handle.close()
            
            # Verificar que los datos no estén vacíos
            if batch_data and len(batch_data.strip()) > 0:
                # Guardar lote
                batch_file = os.path.join(output_dir, f"{db_type}_batch_{batch_num}.{formato}")
                with open(batch_file, "w", encoding="utf-8") as f:
                    f.write(batch_data)
                
                print(f"    Lote {batch_num} guardado: {batch_file}")
            else:
                print(f"    Lote {batch_num} vacío o sin datos")
            
            # Pausa para ser amable con el servidor
            time.sleep(2)
            
        except HTTPError as e:
            print(f"    Error HTTP en lote {batch_num}: {e}")
            continue
        except Exception as e:
            print(f"    Error inesperado en lote {batch_num}: {e}")
            continue

def descargar_metadatos_sra(sra_ids, output_dir, bioproject_id):
    """Descargar metadatos de SRA"""
    
    try:
        # Obtener resúmenes SRA (limitar a 10 para evitar problemas)
        ids_to_fetch = sra_ids[:10]
        handle = Entrez.esummary(db="sra", id=",".join(ids_to_fetch))
        sra_summary = Entrez.read(handle)
        handle.close()
        
        # Guardar metadatos SRA
        sra_file = os.path.join(output_dir, f"sra_metadata_{bioproject_id}.txt")
        with open(sra_file, "w", encoding="utf-8") as f:
            f.write("METADATOS SRA - EXPERIMENTOS DE SECUENCIACIÓN\n")
            f.write("=" * 60 + "\n\n")
            
            for i, sra in enumerate(sra_summary):
                f.write(f"Experimento {i+1}:\n")
                f.write(f"  Accession: {sra.get('Accession', 'N/A')}\n")
                f.write(f"  Title: {sra.get('Title', 'N/A')}\n")
                f.write(f"  Platform: {sra.get('Platform', 'N/A')}\n")
                f.write(f"  Bases: {sra.get('Bases', 'N/A')}\n")
                f.write(f"  Spots: {sra.get('Spots', 'N/A')}\n")
                f.write(f"  Organism: {sra.get('Organism', 'N/A')}\n")
                f.write("-" * 40 + "\n\n")
        
        print(f"  Metadatos SRA guardados: {sra_file}")
        
    except Exception as e:
        print(f"  Error descargando metadatos SRA: {e}")

def descargar_metadatos_pubmed(pubmed_ids, output_dir, bioproject_id):
    """Descargar metadatos de publicaciones"""
    
    try:
        # Obtener resúmenes de PubMed (limitar a 10)
        ids_to_fetch = pubmed_ids[:10]
        handle = Entrez.esummary(db="pubmed", id=",".join(ids_to_fetch))
        pubmed_summary = Entrez.read(handle)
        handle.close()
        
        # Guardar metadatos PubMed
        pubmed_file = os.path.join(output_dir, f"pubmed_metadata_{bioproject_id}.txt")
        with open(pubmed_file, "w", encoding="utf-8") as f:
            f.write("PUBLICACIONES RELACIONADAS - PubMed\n")
            f.write("=" * 50 + "\n\n")
            
            for i, pub in enumerate(pubmed_summary):
                f.write(f"Publicación {i+1}:\n")
                f.write(f"  PMID: {pub.get('Id', 'N/A')}\n")
                f.write(f"  Title: {pub.get('Title', 'N/A')}\n")
                authors = pub.get('AuthorList', ['N/A'])
                if isinstance(authors, list):
                    f.write(f"  Authors: {', '.join(authors[:5])}")  # Primeros 5 autores
                    if len(authors) > 5:
                        f.write(f"... y {len(authors)-5} más")
                    f.write("\n")
                else:
                    f.write(f"  Authors: {authors}\n")
                f.write(f"  Journal: {pub.get('Source', 'N/A')}\n")
                f.write(f"  PubDate: {pub.get('PubDate', 'N/A')}\n")
                f.write(f"  DOI: {pub.get('DOI', 'N/A')}\n")
                f.write("-" * 40 + "\n\n")
        
        print(f"  Metadatos PubMed guardados: {pubmed_file}")
        
    except Exception as e:
        print(f"  Error descargando metadatos PubMed: {e}")

def crear_resumen_bioproject(bioproject_dir, bioproject_id, counts):
    """Crear archivo de resumen del BioProject"""
    
    resumen_file = os.path.join(bioproject_dir, f"RESUMEN_BIOPROJECT_{bioproject_id}.txt")
    
    with open(resumen_file, "w", encoding="utf-8") as f:
        f.write(f"RESUMEN DE DATOS DESCARGADOS - BIOPROJECT {bioproject_id}\n")
        f.write("=" * 70 + "\n\n")
        
        f.write("ESTADÍSTICAS DE DESCARGA:\n")
        f.write(f"- Secuencias de nucleotide: {counts['nucleotide']}\n")
        f.write(f"- Experimentos SRA: {counts['sra']}\n")
        f.write(f"- Secuencias de proteínas: {counts['protein']}\n")
        f.write(f"- Publicaciones relacionadas: {counts['pubmed']}\n")
        
        f.write("\nARCHIVOS DESCARGADOS:\n")
        
        # Listar archivos en cada subcarpeta
        try:
            for root, dirs, files in os.walk(bioproject_dir):
                level = root.replace(bioproject_dir, '').count(os.sep)
                indent = ' ' * 2 * level
                f.write(f"{indent}{os.path.basename(root)}/\n")
                subindent = ' ' * 2 * (level + 1)
                for file in files:
                    if file.endswith(('.fasta', '.txt', '.xml')):
                        file_path = os.path.join(root, file)
                        file_size = os.path.getsize(file_path)
                        f.write(f"{subindent}{file} ({file_size} bytes)\n")
        except Exception as e:
            f.write(f"Error listando archivos: {e}\n")
    
    print(f"   Resumen creado: {resumen_file}")

# =============================================================================
# EJECUCIÓN PRINCIPAL CON MANEJO DE ERRORES 
# =============================================================================

if __name__ == "__main__":
    
    print("SCRIPT DE DESCARGA COMPLETA DE BIOPROJECTS")
    print("Este script descargará todos los datos disponibles de cualquier BioProject")
    
    try:
        # Solicitar BioProject ID al usuario
        bioproject_id = input("\nIngresa el ID del BioProject (ej: PRJNA341863 o 341863): ").strip()
        
        # Asegurar formato correcto
        if bioproject_id and not bioproject_id.startswith('PRJ'):
            # Si solo dan el número, agregar prefijo
            if bioproject_id.isdigit():
                bioproject_id = f"PRJNA{bioproject_id}"
        
        if bioproject_id:
            print(f"\nIniciando descarga para: {bioproject_id}")
            resultado = descargar_bioproject_completo(bioproject_id)
            
            if resultado:
                print(f"\n ¡DESCARGA EXITOSA!")
                print(f"Revisa la carpeta: {resultado}")
            else:
                print(f"\n Hubo problemas con la descarga")
                print("Verifica que el BioProject ID sea correcto")
        else:
            print(" No se proporcionó un ID de BioProject válido")
            
    except KeyboardInterrupt:
        print("\n  Proceso interrumpido por el usuario")
    except Exception as e:
        print(f"\n ERROR inesperado: {e}")
        import traceback
        traceback.print_exc()
    
    print(f"\n{'='*70}")
    print("PROCESO COMPLETADO")
    print(f"{'='*70}")
