import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Leer datos
counts_df = pd.read_csv("data/counts/counts_for_DEA.csv")

# Tomar una muestra y aplicar log-transform
muestra = counts_df.columns[1]  # Primera muestra
datos = np.log2(counts_df[muestra] + 1)

plt.figure(figsize=(10, 6))
sns.histplot(
    data=datos, 
    kde=True,
    bins=30,
    alpha=0.7,
    edgecolor='black'
)
plt.title(f'Distribución de expresión - {muestra}')
plt.xlabel('log2(Conteos + 1)')
plt.ylabel('Frecuencia')
plt.grid(True, alpha=0.3)
plt.savefig('figures/expression_histogram.png', dpi=300, bbox_inches='tight')
plt.show()

