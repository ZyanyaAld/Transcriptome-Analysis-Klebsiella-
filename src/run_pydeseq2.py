#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Script para análisis de expresión diferencial usando pydeseq2
Análisis de transcriptoma de Klebsiella sp. AQSCr - Condición treated vs control
"""

from pathlib import Path
import pandas as pd
import numpy as np
import warnings
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# Configuración de rutas de archivos
COUNTS_PATH = Path("data/counts/counts_for_DEA.csv")
DESIGN_PATH = Path("data/metadata/design_matrix.csv")
OUT_DIR = Path("data/dea")


def main():
    """
    Función principal que ejecuta el análisis completo de expresión diferencial
    """

    # Suprimir warnings de deprecación específicos de pydeseq2
    warnings.filterwarnings("ignore", category=DeprecationWarning, module="pydeseq2")

    print("=" * 60)
    print("ANÁLISIS DE EXPRESIÓN DIFERENCIAL CON PYDESEQ2")
    print("=" * 60)

    # =========================================================================
    # 1. CARGA Y PREPARACIÓN DE DATOS
    # =========================================================================
    print("\n1. CARGANDO Y PREPARANDO DATOS...")

    # Cargar matriz de conteos
    print(f"   • Cargando conteos desde: {COUNTS_PATH}")
    counts_df = pd.read_csv(COUNTS_PATH)

    # Cargar matriz de diseño experimental
    print(f"   • Cargando diseño experimental desde: {DESIGN_PATH}")
    design_df = pd.read_csv(DESIGN_PATH)

    # Verificar que las muestras coincidan entre conteos y diseño experimental
    muestras_diseno = design_df["sample"].tolist()
    print(f"\n   • Muestras en el diseño experimental ({len(muestras_diseno)}):")
    for muestra in muestras_diseno:
        print(f"     - {muestra}")

    # Preparar matriz de conteos para DESeq2
    # Filtrar y ordenar conteos según el diseño experimental
    counts_filtrados = counts_df[["Geneid"] + muestras_diseno].set_index("Geneid")

    # Transponer: DESeq2 espera muestras como filas y genes como columnas
    counts_para_deseq = counts_filtrados.T

    # Preparar metadata ordenada según las muestras
    design_df = design_df.set_index("sample").loc[muestras_diseno]

    # Mostrar distribución de condiciones experimentales
    print(f"\n   • Distribución de condiciones:")
    conteo_condiciones = design_df["condition"].value_counts()
    for condicion, conteo in conteo_condiciones.items():
        print(f"     - {condicion}: {conteo} muestras")

    # =========================================================================
    # 2. ANÁLISIS DE EXPRESIÓN DIFERENCIAL CON DESEQ2
    # =========================================================================
    print("\n2. INICIANDO ANÁLISIS DE EXPRESIÓN DIFERENCIAL...")

    # Crear objeto DeseqDataSet con sintaxis moderna si se puede
    try:
        dds = DeseqDataSet(
            counts=counts_para_deseq.astype(int),  # Conteos deben ser enteros
            metadata=design_df,
            design="~ condition",  # Fórmula moderna para el diseño
        )
        print("   • Usando sintaxis moderna de pydeseq2")
    except (TypeError, ValueError):
        # Usar sintaxis antigua si falla la moderna
        print("   • Usando sintaxis antigua (compatibilidad)")
        dds = DeseqDataSet(
            counts=counts_para_deseq.astype(int),
            metadata=design_df,
            design_factors=["condition"],  # <- LISTA, lo que sabemos que funciona
        )

    # Ejecutar pipeline completo de DESeq2
    print("\n   • Ejecutando pipeline DESeq2...")
    print("     - Ajustando factores de tamaño...")
    print("     - Estimando dispersiones...")
    print("     - Ajustando modelo lineal...")

    dds.deseq2()

    # =========================================================================
    # 3. PRUEBAS ESTADÍSTICAS Y RESULTADOS
    # =========================================================================
    print("\n3. REALIZANDO PRUEBAS ESTADÍSTICAS...")

    # Configurar contraste: treated vs control
    estadisticas = DeseqStats(dds, contrast=("condition", "treated", "control"))

    # Ejecutar prueba de Wald
    estadisticas.run_wald_test()

    # Generar resumen estadístico
    print("   • Resumen de resultados de expresión diferencial:")
    estadisticas.summary()

    # Obtener dataframe de resultados
    resultados = estadisticas.results_df
    resultados.index.name = "Geneid"

    # Crear directorio de salida si no existe
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # Guardar resultados completos
    ruta_resultados = OUT_DIR / "deseq2_results.csv"
    resultados.to_csv(ruta_resultados)
    print(f"\n   • Resultados guardados en: {ruta_resultados}")

    # =========================================================================
    # 4. OBTENER Y GUARDAR CONTEOS NORMALIZADOS
    # =========================================================================
    print("\n4. PROCESANDO CONTEOS NORMALIZADOS...")
    conteos_normalizados_df = None

    # Estrategia 1: Buscar en capas (layers)
    if hasattr(dds, "layers"):
        for nombre_capa in ["normed_counts", "normalized_counts", "counts"]:
            if nombre_capa in dds.layers:
                array_normalizado = dds.layers[nombre_capa]
                print(f"   • Conteos normalizados encontrados en capa: '{nombre_capa}'")
                conteos_normalizados_df = pd.DataFrame(
                    array_normalizado.T,
                    index=counts_filtrados.index,
                    columns=muestras_diseno,
                )
                break

    # Estrategia 2: Buscar en atributos directos
    if conteos_normalizados_df is None:
        for nombre_atributo in ["normed_counts", "normalized_counts"]:
            if hasattr(dds, nombre_atributo):
                objeto_conteos = getattr(dds, nombre_atributo)
                if isinstance(objeto_conteos, pd.DataFrame):
                    conteos_normalizados_df = objeto_conteos
                    print(
                        f"   • Conteos normalizados encontrados en atributo: '{nombre_atributo}'"
                    )
                    break
                elif isinstance(objeto_conteos, np.ndarray):
                    conteos_normalizados_df = pd.DataFrame(
                        objeto_conteos.T,
                        index=counts_filtrados.index,
                        columns=muestras_diseno,
                    )
                    print(
                        f"   • Conteos normalizados convertidos desde: '{nombre_atributo}'"
                    )
                    break

    # Estrategia 3: Cálculo manual usando factores de tamaño
    if (
        conteos_normalizados_df is None
        and hasattr(dds, "counts")
        and hasattr(dds, "size_factors")
    ):
        print("   • Calculando conteos normalizados manualmente...")
        try:
            conteos_crudos = dds.counts.T  # genes x muestras
            if dds.size_factors is not None and len(dds.size_factors) == conteos_crudos.shape[1]:
                conteos_normalizados_df = conteos_crudos.div(dds.size_factors, axis=1)
                print("   ✓ Conteos normalizados calculados exitosamente")
            else:
                print("   ✗ Factores de tamaño no disponibles o con longitud incorrecta")
        except Exception as error:
            print(f"   ✗ Error en cálculo manual: {error}")

    # Guardar conteos normalizados si se encontraron
    if conteos_normalizados_df is not None:
        ruta_conteos_normalizados = OUT_DIR / "normalized_counts.csv"
        conteos_normalizados_df.to_csv(ruta_conteos_normalizados)
        print(f"   • Conteos normalizados guardados en: {ruta_conteos_normalizados}")
        print(f"   • Estadísticas de conteos normalizados:")
        print(f"     - Forma: {conteos_normalizados_df.shape} (genes × muestras)")
        print(
            f"     - Rango: {conteos_normalizados_df.values.min():.2f} - "
            f"{conteos_normalizados_df.values.max():.2f}"
        )
        print(
            f"     - Mediana: {np.median(conteos_normalizados_df.values):.2f}"
        )
    else:
        print("   • No se pudieron obtener conteos normalizados")
        print("   • Los resultados DE están completos, pero faltan conteos normalizados")

    # =========================================================================
    # 5. ANÁLISIS DE RESULTADOS Y ESTADÍSTICAS
    # =========================================================================
    print("\n5. RESUMEN DE RESULTADOS...")
    print("=" * 50)

    # Top 10 genes más significativos
    top_genes = resultados.dropna(subset=["padj"]).sort_values("padj").head(10)
    print("\n   • TOP 10 GENES MÁS SIGNIFICATIVOS (menor padj):")
    print("     " + "-" * 70)
    for i, (gen_id, fila) in enumerate(top_genes.iterrows(), 1):
        print(f"     {i:2d}. {gen_id}")
        print(f"         Expresión media: {fila['baseMean']:8.1f}")
        print(f"         Fold-change (log2): {fila['log2FoldChange']:7.3f}")
        print(f"         p-valor: {fila['pvalue']:10.2e}")
        print(f"         p-ajustado: {fila['padj']:10.2e}")

    # Estadísticas generales del análisis
    print("\n   • ESTADÍSTICAS GENERALES DEL ANÁLISIS:")
    print("     " + "-" * 40)

    total_genes = len(resultados)
    genes_significativos = len(resultados[resultados["padj"] < 0.05])
    genes_sobrerregulados = len(
        resultados[(resultados["padj"] < 0.05) & (resultados["log2FoldChange"] > 0)]
    )
    genes_subregulados = len(
        resultados[(resultados["padj"] < 0.05) & (resultados["log2FoldChange"] < 0)]
    )

    print(f"     • Total de genes analizados: {total_genes}")
    print(
        f"     • Genes diferencialmente expresados (padj < 0.05): {genes_significativos}"
    )
    if genes_significativos > 0:
        print(f"        - Sobreexpresados: {genes_sobrerregulados}")
        print(f"        - Subexpresados: {genes_subregulados}")
        print(f"        - Porcentaje: {(genes_significativos/total_genes)*100:.1f}%")

    genes_fc_alto = len(
        resultados[
            (resultados["padj"] < 0.05)
            & (abs(resultados["log2FoldChange"]) > 2)
        ]
    )
    if genes_fc_alto > 0:
        print(f"     • Genes con |log2FC| > 2: {genes_fc_alto}")

    print("\n" + "=" * 60)
    print("ANÁLISIS COMPLETADO EXITOSAMENTE")
    print("=" * 60)


if __name__ == "__main__":
    main()
