import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os


import os
print("Directorio actual:", os.getcwd())
# Crear carpeta figures si no existe
os.makedirs("../figures", exist_ok=True)

# 1. Leer archivos
counts_df = pd.read_csv("../data/counts/counts_for_DEA.csv")      # Conteos de genes
design_df = pd.read_csv("../data/metadata/design_matrix.csv")     # Condiciones experimentales

counts_df.columns = counts_df.columns.str.strip()   # quitar espacios

# 3. Reorganizar datos: melt
datos_melted = counts_df.melt(
    id_vars=['Geneid'],
    var_name='sample',
    value_name='count'
)

# 4. Unir con diseño experimental
# design_df = design_df.rename(columns={'SampleID': 'sample'})  # si fuera necesario
datos_melted = datos_melted.merge(design_df, on='sample')

# 5. Calcular log2(count+1)
datos_melted['log_count'] = np.log2(datos_melted['count'] + 1)

# 6. Graficar Boxplot y Violinplot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
sns.boxplot(data=datos_melted, x='condition', y='log_count',
            palette={'control': 'lightblue', 'treated': 'lightcoral'}, ax=ax1)
ax1.set_title('Boxplot - Distribución por condición')
ax1.set_ylabel('log2(Conteos + 1)')

sns.violinplot(data=datos_melted, x='condition', y='log_count',
               palette={'control': 'lightblue', 'treated': 'lightcoral'}, ax=ax2)
ax2.set_title('Violinplot - Densidad por condición')
ax2.set_ylabel('log2(Conteos + 1)')
plt.tight_layout()
plt.savefig("../figures/box_violin.png", dpi=300)
plt.close()

# Scatterplot entre dos muestras
muestra1 = counts_df.columns[1]     # Primera muestra, puedes cambiar para que analice otras muestras.
muestra2 = counts_df.columns[2]    # Segunda muestra, puedes cambiar para que analice otras muestras.
plt.figure(figsize=(8, 6))
sns.scatterplot(data=counts_df, x=muestra1, y=muestra2, alpha=0.6, s=30)
plt.title(f'Relación entre {muestra1} vs {muestra2}')
plt.xlabel(f'Conteos {muestra1}')
plt.ylabel(f'Conteos {muestra2}')
plt.xscale('log')
plt.yscale('log')
plt.savefig("../figures/scatter_samples.png", dpi=300)
plt.close()

# Histograma de una muestra
muestra = counts_df.columns[3] # Primera muestra, puedes cambiar para que analice otras muestras.
datos = np.log2(counts_df[muestra] + 1)
plt.figure(figsize=(10, 6))
sns.histplot(data=datos, kde=True, bins=30, alpha=0.7, edgecolor='black')
plt.title(f'Distribución de expresión - {muestra}')
plt.xlabel('log2(Conteos + 1)')
plt.ylabel('Frecuencia')
plt.grid(True, alpha=0.3)
plt.savefig("../figures/hist_muestra.png", dpi=300)
plt.close()

# Pairplot
muestras_subset = counts_df.columns[1:5]   #Pueden ser mas o menos muestras
datos_subset = counts_df[muestras_subset]
datos_log = np.log2(datos_subset + 1)
sns.pairplot(datos_log, diag_kind='kde', plot_kws={'alpha': 0.6})
plt.suptitle('Relaciones entre Muestras - Pairplot', y=1.02)
plt.savefig("../figures/pairplot.png", dpi=300)
plt.close()

# Heatmap de correlación
counts_data = counts_df.set_index('Geneid')
corr_matrix = counts_data.corr()
plt.figure(figsize=(10, 8))
sns.heatmap(corr_matrix, annot=True, fmt='.2f', cmap='coolwarm', center=0.8, square=True)
plt.title('Correlación entre Muestras - Klebsiella sp. AQSCr')
plt.tight_layout()
plt.savefig("../figures/heatmap_corr.png", dpi=300)
plt.close()

# Top 10 genes más expresados
mean_expression = counts_data.mean(axis=1)
top_10 = mean_expression.nlargest(10).reset_index()
top_10.columns = ['Gene', 'Mean_Expression']
plt.figure(figsize=(12, 6))
sns.barplot(data=top_10, x='Mean_Expression', y='Gene', palette='viridis')
plt.title('Top 10 Genes Más Expresados - Klebsiella sp. AQSCr')
plt.xlabel('Expresión Promedio')
plt.ylabel('Gen')
plt.tight_layout()
plt.savefig("../figures/top10_genes.png", dpi=300)
plt.close()

# Lineplot de genes 
genes_ejemplo = counts_df['Geneid'].head(10)  #Pueden ser mas o menos genes
line_data = counts_df[counts_df['Geneid'].isin(genes_ejemplo)]
line_data = line_data.melt(id_vars=['Geneid'], var_name='sample', value_name='count')
line_data = line_data.merge(design_df, on='sample')
plt.figure(figsize=(12, 6))
sns.lineplot(data=line_data, x='sample', y='count', hue='Geneid', marker='o', markersize=8)
plt.title('Expresión de Genes por Muestra')
plt.xlabel('Muestra')
plt.ylabel('Conteos')
plt.xticks(rotation=45)
plt.legend(title='Gen')
plt.tight_layout()
plt.savefig("../figures/lineplot_genes.png", dpi=300)
plt.close()

# Distribución de genes por categoría
bins = [0, 10, 100, 1000, float('inf')]
labels = ['Muy bajo', 'Bajo', 'Medio', 'Alto']
counts_data['categoria'] = pd.cut(mean_expression, bins=bins, labels=labels)
categoria_counts = []
for condition in design_df['condition'].unique():
    samples = design_df.loc[design_df['condition'] == condition, 'sample']
    samples = [s for s in samples if s in counts_data.columns]
    for categoria in labels:
        mask_categoria = counts_data['categoria'] == categoria
        mask_expresados = counts_data[samples].mean(axis=1) > 0
        count = (mask_categoria & mask_expresados).sum()
        categoria_counts.append({'condition': condition, 'categoria': categoria, 'count': count})
stack_data = pd.DataFrame(categoria_counts)
plt.figure(figsize=(10, 6))
sns.barplot(data=stack_data, x='condition', y='count', hue='categoria', palette='Set2')
plt.title('Distribución de Genes por Categoría de Expresión')
plt.ylabel('Número de Genes')
plt.xlabel('Condición')
plt.legend(title='Categoría de Expresión', bbox_to_anchor=(1.05, 1))
plt.tight_layout()
plt.savefig("../figures/categorias_expr.png", dpi=300)
plt.close()

# Expresión acumulada por percentil
percentiles = np.arange(0, 101, 10)
accum_data = []
for condition in design_df['condition'].unique():
    samples = design_df.loc[design_df['condition'] == condition, 'sample']
    samples = [s for s in samples if s in counts_data.columns]
    condition_data = counts_data[samples].mean(axis=1)
    for p in percentiles:
        threshold = np.percentile(condition_data, p)
        accum_expression = condition_data[condition_data <= threshold].sum()
        accum_data.append({'condition': condition, 'percentile': p, 'accum_expression': accum_expression})
accum_df = pd.DataFrame(accum_data)
plt.figure(figsize=(10, 6))
sns.lineplot(data=accum_df, x='percentile', y='accum_expression',
             hue='condition', palette={'control': 'blue', 'treated': 'red'}, marker='o')
plt.title('Expresión Acumulada por Percentil')
plt.xlabel('Percentil')
plt.ylabel('Expresión Acumulada')
plt.grid(True, alpha=0.3)
plt.legend(title='Condición')
plt.tight_layout()
plt.savefig("../figures/percentiles_expr.png", dpi=300)
plt.close()