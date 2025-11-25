import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# 1. Leer archivos
counts_df = pd.read_csv("data/counts/counts_for_DEA.csv")
design_df = pd.read_csv("data/metadata/design_matrix.csv")

# 2. Reorganizar datos: melt
datos_melted = counts_df.melt(
    id_vars=['Geneid'],
    var_name='sample',
    value_name='count'
)

# 3. Unir con diseño experimental
datos_melted = datos_melted.merge(design_df, on='sample')

# 4. Calcular log2(count+1)
datos_melted['log_count'] = np.log2(datos_melted['count'] + 1)

# 5. Graficar
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

# Boxplot
sns.boxplot(
    data=datos_melted,
    x='condition',
    y='log_count',
    hue='condition',
    palette={'control': 'lightblue', 'treated': 'lightcoral'},
    ax=ax1,
    legend=False
)
ax1.set_title('Boxplot - Distribución por condición')
ax1.set_ylabel('log2(Conteos + 1)')

# Violinplot
sns.violinplot(
    data=datos_melted,
    x='condition',
    y='log_count',
    hue='condition',
    palette={'control': 'lightblue', 'treated': 'lightcoral'},
    ax=ax2,
    legend=False
)
ax2.set_title('Violinplot - Densidad por condición')
ax2.set_ylabel('log2(Conteos + 1)')

plt.tight_layout()
plt.savefig('figures/distribution_plots.png', dpi=300, bbox_inches='tight')
plt.show()

