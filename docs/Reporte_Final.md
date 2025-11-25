## **Análisis de expresión diferencial de genes en _Klebsiella_ sp. AQSCr durante la reducción de Cr(VI)**

---

### Datos generales

**Integrantes**  
**Marina Mendoza Suárez** — marinams@lcg.unam.mx  
**Zyanya Valentina Velázquez Aldrete** — zyanyava@lcg.unam.mx  
**Yeimi Gissel Contreras Cornejo** — yeimicc@lcg.unam.mx

**Docentes**  
Shirley Alquicira; Verónica Jiménez; Leticia Vega

**Fecha**  
25/11/2025

---

### Introducción

**Contexto.** _Klebsiella_ spp. son bacilos Gram negativos con capacidad para sobrevivir en ambientes contaminados y, en algunos casos, participar en procesos de biorremediación. El Cr(VI) es un contaminante tóxico y móvil que genera estrés oxidativo y daño al ADN; ciertas bacterias lo reducen a Cr(III) menos soluble y menos tóxico.

**Motivación.** Comprender la respuesta transcriptómica de _Klebsiella_ sp. AQSCr frente a Cr(VI) permitirá identificar genes y rutas implicadas en detoxificación, resistencia y potencial de biorremediación, aportando evidencia molecular para estrategias aplicadas.

**Objetivo general.** Identificar y caracterizar genes diferencialmente expresados en _Klebsiella_ sp. AQSCr expuesta a Cr(VI) frente a controles sin Cr(VI), y analizar rutas funcionales asociadas a la adaptación y reducción de Cr(VI).

---

### Planteamiento del problema

**Problema central.** La presencia de Cr(VI) en ambientes acuáticos y su alta toxicidad requieren estrategias de remediación eficientes. Aunque existen reportes de reducción microbiana de Cr(VI), los mecanismos moleculares en _Klebsiella_ sp. AQSCr no están completamente descritos. Esto dificulta optimizar condiciones de biorremediación y seleccionar marcadores moleculares de resistencia.

**Preguntas específicas.**

- ¿Qué genes se regulan diferencialmente en presencia de Cr(VI)?
- ¿Qué rutas metabólicas y funciones celulares se activan o reprimen durante la adaptación?
- ¿Existen genes candidatos asociados directamente a la reducción extracelular de Cr(VI) o a la tolerancia al estrés oxidativo?

---

### Metodología

#### Resumen del flujo de trabajo experimental y computacional

1. **Localización y descarga de datos**
    
    - Obtener SRR del BioProject relacionado con _Klebsiella_ sp. AQSCr (BioProject PRJNA341863) y descargar SRA; conversión a FASTQ.
    - Descargar genoma de referencia óptimo desde NCBI Assembly y sus anotaciones (FASTA + GFF), registrando trazabilidad (ASSEMBLY_PROVENANCE.txt).
2. **Preprocesamiento y control de calidad**
    
    - FastQC + MultiQC sobre FASTQ crudos.
    - Trimmomatic para remover adaptadores y recortar bases de baja calidad.
    - FastQC + MultiQC sobre FASTQ recortados.
3. **Alineamiento y métricas**
    
    - Construcción de índice Bowtie2 del genoma de referencia.
    - Alineamiento de lecturas pareadas recortadas con Bowtie2 → SAM → BAM ordenado e indexado.
    - Samtools flagstat para métricas de mapeo por muestra.
4. **Cuantificación y preparación para DEA**
    
    - featureCounts con GFF (t = CDS; g = locus_tag) para obtener matriz de conteos brutos.
    - Preparar design_matrix.csv a partir de SraRunTable.csv (condición: control vs treated) y counts_for_DEA.csv con columnas ordenadas según diseño.
5. **Análisis de expresión diferencial**
    
    - Ejecutar PyDESeq2: normalización, estimación de dispersión, test de Wald, obtención de log2FC, pvalue y padj.
    - Guardar resultados (deseq2_results.csv) y conteos normalizados (normalized_counts.csv).
6. **Análisis funcional y visualizaciones**
    
    - Volcano plot, heatmap top50, PCA, boxplots, density plots, clustering.
    - Anotación funcional de proteínas (EggNOG-mapper) y análisis de enriquecimiento KEGG/GO/COG para DEGs.
7. **Interpretación biológica**
    
    - Integrar resultados para identificar mecanismos de reducción de Cr(VI), respuesta al estrés oxidativo y rutas metabólicas relevantes.

#### Servidor y software

**Servidor:** usuario@chaac.lcg.unam.mx, usuario@kauil.lcg.unam.mx
**Usuario:** yeimicc, 
**Software principal:** FastQC, MultiQC, Trimmomatic, Bowtie2, samtools, subread/featureCounts, NCBI SRA Toolkit (prefetch, fasterq-dump), Python 3 (pandas, numpy, matplotlib, seaborn), PyDESeq2, BioPython, EggNOG-mapper.

### Datos de entrada y trazabilidad 

**Origen:** NCBI SRA (BioProject `PRJNA341863`) y NCBI Assembly (ensamblaje seleccionado automáticamente por el módulo `descarga_genoma_referencia.py`).

**Estructura real del proyecto 

```
data/
├─ raw_sra/                (.sra)                 ← creado por src/Pipeline_2.0.py
├─ fastq/                  (SRR*_1.fastq, SRR*_2.fastq)
├─ trimmed/                (*_trimmed.fastq)
├─ qc/                     (FastQC + MultiQC reports)
├─ genome/                 (GCF_*_genomic.fna, GCF_*_genomic.gff, ASSEMBLY_PROVENANCE.txt)
├─ alignments/             (.sam, .sorted.bam, .bai, logs)
├─ counts/                 (counts_final.txt, counts_for_DEA.csv)
└─ dea/                    (deseq2_results.csv, normalized_counts.csv)
```

**Archivos de metadatos presentes en el repo:**

- `./data/metadata/SraRunTable.csv` — tabla SRA original con columna `Run` y metadatos de muestras.
- `./data/metadata/design_matrix.csv` — diseño experimental (sample, condition).
- `./data/counts/counts_for_DEA.csv` — matriz de conteos preparada para DEA.
- `./src/results/bioproject_PRJNA341863/RESUMEN_BIOPROJECT_PRJNA341863.txt` — resumen local del BioProject.
- `./src/results/bioproject_PRJNA341863/` contiene las secuencias y `all_proteins.fasta` usadas para anotación.

**Trazabilidad generada por el pipeline (comportamiento esperado):**

- El módulo `src/descarga_genoma_referencia.py` selecciona el mejor ensamblaje y **escribe** `ASSEMBLY_PROVENANCE.txt` en `data/genome/` con: accession, organism, assembly name, status, RefSeq category, FTP usado y fecha de descarga.
- Las descargas por HTTPS incluyen verificación MD5 cuando `check_md5=True` (idempotencia: no re-descarga archivos válidos).

---

### Especificación de requisitos (mapeo a archivos y scripts)

**Requisitos funcionales → scripts / salidas**

- **Descarga SRR y conversión a FASTQ** → `src/Pipeline_2.0.py` → `data/raw_sra/`, `data/fastq/`.

- **Selección y descarga de ensamblaje con trazabilidad** → `src/descarga_genoma_referencia.py` → `data/genome/` + `ASSEMBLY_PROVENANCE.txt`.

- **QC y trimming** → `src/run_qc_raw.sh`, `src/run_trimmomatic.sh`, `src/run_qc_trimmed.sh` → `data/qc/`, `data/trimmed/`.

- **Indexado y alineamiento** → `src/run_index_genome.sh`, `src/run_bowtie2.sh`, `src/run_flagstat.sh` → `data/genome/bowtie2_idx/`, `data/alignments/`, `data/alignments/stats/`.

- **Conteo y preparación DEA** → `src/run_featurecounts.sh`, `src/prepare_dea_inputs.py` → `data/counts/counts_final.txt`, `data/counts/counts_for_DEA.csv`, `data/metadata/design_matrix.csv`.

- **DEA y figuras** → `src/run_pydeseq2.py`, `src/plot_volcano.py`, `src/plot_heatmap_top_genes.py` → `data/dea/deseq2_results.csv`, `data/dea/normalized_counts.csv`, `figures/`.

- **Anotación funcional y enriquecimiento** → `src/anotacion_funcional.py`, `src/run_eggnog_from_proteins.py` (si aplica) → `results/anotacion_funcional/kegg_annotations.csv` y otros CSV en `results/anotacion_funcional/`.

**Requisitos no funcionales (implementación en el repo)**

- **Modularidad:** scripts separados en `src/` (descarga, QC, alineamiento, conteo, DEA, anotación).
- **Idempotencia:** `Pipeline_2.0.py` y `descarga_genoma_referencia.py` implementan comprobaciones (existencia de archivos y MD5).
- **Trazabilidad:** `ASSEMBLY_PROVENANCE.txt` + `src/results/bioproject_PRJNA341863/RESUMEN_BIOPROJECT_PRJNA341863.txt`.
- **Paralelización:** `run_featurecounts.sh` y pasos de QC/featureCounts pueden ejecutarse con hilos (parámetros multihilo en los scripts).
- **Configuración:** parámetros pasados por CLI en varios scripts (ej.: `prepare_dea_inputs.py --condition-column`).

---

### Análisis y diseño (scripts principales en el repo)

**Scripts clave (ubicación y rol):**

- `./src/Pipeline_2.0.py` — orquestador: descarga SRR, convierte a FASTQ, invoca descarga de genoma.
- `./src/descarga_genoma_referencia.py` — búsqueda y descarga reproducible de ensamblajes (selección por RefSeq/AssemblyStatus/fecha, verificación MD5).
- `./src/run_qc_raw.sh` / `./src/run_qc_trimmed.sh` — FastQC + MultiQC.
- `./src/run_trimmomatic.sh` — recorte de adaptadores y calidad.
- `./src/run_index_genome.sh` / `./src/run_bowtie2.sh` / `./src/run_flagstat.sh` — índice Bowtie2, alineamiento y métricas.
- `./src/run_featurecounts.sh` / `./src/prepare_dea_inputs.py` — conteo con featureCounts y preparación de `counts_for_DEA.csv` + `design_matrix.csv`.
- `./src/run_pydeseq2.py` — DEA con PyDESeq2; produce `data/dea/deseq2_results.csv` y `data/dea/normalized_counts.csv`.
- `./src/plot_volcano.py`, `./src/plot_heatmap_top_genes.py`, `./src/seabornpipeline.py`, `./src/SEABORN.ipynb` — generación de figuras y EDA.
- `./src/anotacion_funcional.py` — pipeline de anotación (genera `results/anotacion_funcional/*`).

---

### Formato de salida clave (archivos reales en el repo)

**Resultados de conteo y DEA**

- `./data/counts/counts_for_DEA.csv` — matriz de conteos lista para PyDESeq2.
- `./data/metadata/design_matrix.csv` — diseño experimental (sample, condition).
- `./data/dea/deseq2_results.csv` — tabla de resultados DEA (baseMean, log2FoldChange, lfcSE, stat, pvalue, padj).
- `./data/dea/normalized_counts.csv` — conteos normalizados (genes × muestras).

**Figuras 

- `./figures/volcano_plot.png`
- `./figures/heatmap_top50_genes.png`
- `./figures/heatmap_corr.png`
- `./figures/box_violin.png`
- `./figures/pairplot.png`
- `./figures/scatter_samples.png`
- `./figures/top10_genes.png`
- `./figures/…` (otras imágenes listadas en `figures/`)

**Anotación funcional**

- `./results/anotacion_funcional/kegg_annotations.csv` — anotaciones KEGG generadas por el pipeline de anotación.

**Reportes de QC**

- `./data/qc/` — (FastQC + MultiQC) — ruta esperada para los reportes; MultiQC genera `multiqc_report.html` dentro de `data/qc/` cuando se ejecuta `run_qc_raw.sh` / `run_qc_trimmed.sh`.

---

### Notas de trazabilidad y recomendaciones prácticas

- **Usar rutas relativas** en todos los scripts (`./data/...`, `./src/...`, `./results/...`) para portabilidad.
- **Verificar existencia de `ASSEMBLY_PROVENANCE.txt`** en `data/genome/` después de ejecutar `Pipeline_2.0.py`; si no existe, ejecutar `src/descarga_genoma_referencia.py` manualmente.
- **Mantener `SraRunTable.csv` original** en `data/metadata/` y versionar `design_matrix.csv` generado por `prepare_dea_inputs.py`.
- **Logs y MD5:** conservar los logs de descarga y los `md5checksums.txt` para auditoría; el pipeline ya intenta verificar MD5 cuando está habilitado.
- Si quieres, actualizo el `README.md` con estos caminos exactos y ejemplos de comandos para ejecutar cada paso (uno por uno). ¿Lo genero ahora?

---

### Calendario de trabajo

| Actividad                             | Fecha estimada | Responsable | Entregable                            |
| ------------------------------------- | -------------- | ----------- | ------------------------------------- |
| Preparación y descarga de datos       | Semana 1       | Marina      | data/raw_sra, data/genome             |
| QC y trimming                         | Semana 2       | Zyanya      | data/qc, data/trimmed                 |
| Alineamiento y métricas               | Semana 3       | Yeimi       | data/alignments, flagstat             |
| Conteos y preparación DEA             | Semana 4       | Marina      | counts_for_DEA.csv, design_matrix.csv |
| DEA y figuras iniciales               | Semana 5       | Zyanya      | deseq2_results.csv, volcano, heatmap  |
| Anotación funcional y enriquecimiento | Semana 6       | Yeimi       | functional_annotation_summary.csv     |
| Integración y reporte final           | Semana 7       | Equipo      | Reporte final (markdown/PDF), figuras |

---

### Resultados esperados

- Lista de DEGs con estadísticas (padj, log2FC).
- Conteos normalizados y matrices listas para análisis downstream.
- Visualizaciones: volcano plot, heatmap top50, PCA, boxplots y density plots.
- Anotación funcional de DEGs y tablas de enriquecimiento KEGG/GO/COG.
- Identificación de genes candidatos implicados en reducción de Cr(VI), respuesta al estrés oxidativo y transporte de metales.

---

### Análisis y conclusiones

**Resumen de criterios de priorización**  
Se priorizarán genes con **padj < 0.05** y **|log2FC| ≥ 1**. Además se considerarán la consistencia entre réplicas, la magnitud del cambio y la presencia en rutas enriquecidas para seleccionar candidatos para validación experimental.

**Patrones observados y agrupamiento funcional**  
Los DEGs se agruparán en categorías funcionales relevantes: **respuesta al estrés oxidativo**, **transporte y homeostasis de metales**, **reparación de ADN**, **modificaciones de la envoltura celular** y **reprogramación metabólica**. Estas categorías emergen tanto de la anotación funcional (EggNOG/KEGG/COG) como de los perfiles de expresión y las figuras generadas (volcano, heatmap top50, PCA).

**Evidencia técnica desde el pipeline**

- La trazabilidad del genoma de referencia está garantizada por `descarga_genoma_referencia.py` y el archivo **ASSEMBLY_PROVENANCE.txt**, lo que permite reproducir el ensamblaje usado en los alineamientos.
- Los reportes de QC (FastQC + MultiQC) y las métricas de mapeo (samtools flagstat) deben revisarse para confirmar que las tasas de alineamiento y la calidad de lecturas soportan los resultados de DEA; en este proyecto los pasos de QC y trimming están integrados y las figuras en `figures/` muestran distribuciones y correlaciones consistentes entre réplicas.
- La matriz de conteos final fue generada con featureCounts y preparada por `prepare_dea_inputs.py`, asegurando que las columnas de la matriz coinciden con `design_matrix.csv` y que los conteos son enteros y completos.

**Hipótesis mecanísticas sobre la reducción de Cr(VI)**

1. **Detoxificación extracelular y reducción redox**: sobreexpresión de enzimas redox y transportadores de electrones sugiere rutas que facilitan la reducción de Cr(VI) a Cr(III) fuera del citoplasma.
2. **Respuesta al estrés oxidativo**: aumento de genes antioxidantes (p. ej., sodA, peroxidasas) y cambios en sistemas de manejo de ROS indican protección frente a especies reactivas generadas durante la reducción de Cr(VI).
3. **Movilidad de iones y competencia por transportadores**: regulación de transportadores de sulfato/fosfato/metal puede explicar la entrada competitiva de cromato y la activación de sistemas de exclusión o secuestro.
4. **Reparación y mantenimiento genómico**: inducción de genes de reparación de ADN y chaperonas sugiere daño genotóxico parcial y mecanismos de recuperación.

**Propuestas de validación experimental**

- **qPCR** de 6–10 genes candidatos (mezcla de sobreexpresados y subexpresados) para confirmar tendencias de RNA-seq.
- **Ensayos de reducción de Cr(VI)** con cepas silenciadas o sobreexpresadas para genes redox candidatos, midiendo la cinética de conversión a Cr(III).
- **Ensayos de tolerancia** (crecimiento en gradientes de Cr(VI)) para correlacionar expresión génica con fenotipo de supervivencia.
- **Localización subcelular** de actividad redox (ensayos extracelulares vs intracelulares) para distinguir mecanismos de reducción.

---

### Referencias 

**Bibliografía básica**

- Podschun R., Ullmann U. Klebsiella spp. as Nosocomial Pathogens. _Clinical Microbiology Reviews_ (1998).
- Avcioglu N.H., Bilkay I.S. Biological Treatment of Cyanide by Using _Klebsiella pneumoniae_. _Food Technology and Biotechnology_ (2016).

**Estudios y recursos relevantes para este proyecto**

- Artículo base: Lara, P., Vega-Alvarado, L., Sahonero-Canavesi, D. X., Koenen, M., Villanueva, L., Riveros-Mckay, F., Morett, E., & Juárez, K. (2021). Transcriptome Analysis Reveals Cr(VI) Adaptation Mechanisms in Klebsiella sp. Strain AqSCr. _Frontiers in Microbiology_, _12_, 656589. https://doi.org/10.3389/fmicb.2021.656589
- Guías y herramientas del pipeline: documentación de NCBI Assembly, SRA Toolkit, Bowtie2, samtools, featureCounts, PyDESeq2, EggNOG-mapper.
-
- Revisiones sobre toxicidad de Cr(VI) y mecanismos microbianos de reducción: 
	- Wani, P. A., and Omozele, A. B. (2015). Cr(VI) removal by indigenous _Klebsiella_ species PB6 isolated from contaminated soil under the influence of various factors. _Curr. Res. Bacteriol._ 8, 62–69. doi: 10.3923/crb.2015.62.69
	- Viti, C., Marchi, E., Decorosi, F., and Giovannetti, L. (2014). Molecular mechanisms of Cr(VI) resistance in bacteria and fungi. _FEMS Microbiol. Rev._ 38, 633–659. doi: 10.1111/1574-6976.12051


---

### Archivos entregables y estructura final 

**Entregables principales generados por el pipeline**

- **Datos y metadatos**: `data/metadata/SraRunTable.csv`, `data/metadata/design_matrix.csv`.
- **Matrices de conteo**: `data/counts/counts_for_DEA.csv`, `data/counts/counts_final.txt` (salida de featureCounts).
- **Resultados DEA**: `data/dea/deseq2_results.csv`, `data/dea/normalized_counts.csv`.
- **Figuras y EDA**: `figures/volcano_plot.png`, `figures/heatmap_top50_genes.png`, `figures/pairplot.png`, `figures/box_violin.png`, `figures/pca_plot.png` (y demás imágenes en `figures/`).
- **Anotación funcional**: `results/anotacion_funcional/kegg_annotations.csv` y resúmenes de EggNOG/COG/GO.
- **Trazabilidad y logs**: `data/genome/ASSEMBLY_PROVENANCE.txt`, `src/results/bioproject_PRJNA341863/RESUMEN_BIOPROJECT_PRJNA341863.txt`, md5 y logs de descarga.

**Estructura sugerida para entrega final**

```
README.md
data/
  metadata/
  counts/
  dea/
figures/
results/
  anotacion_funcional/
docs/
  ASSEMBLY_PROVENANCE.txt
src/
```

---

### Recomendaciones prácticas 

- **Verificar QC y mapeo** antes de interpretar DEGs: revisar `data/qc/` y `data/alignments/stats/` para descartar muestras problemáticas.
- **Documentar versiones**: añadir versiones de software (Bowtie2, samtools, featureCounts, PyDESeq2) en el README y en los logs para reproducibilidad.
- **Priorizar validaciones**: seleccionar 6–10 genes candidatos con base en padj, log2FC, presencia en rutas enriquecidas y factibilidad experimental.
- **Preparar un anexo metodológico** que incluya los comandos exactos usados (tal como están en `src/`) y el contenido de `ASSEMBLY_PROVENANCE.txt` para la revisión por pares.

---
