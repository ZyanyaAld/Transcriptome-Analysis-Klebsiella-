import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

# Leer datos
counts_df = pd.read_csv("data/counts/counts_for_DEA.csv")

# Preparar datos: expresión de 2 muestras para scatter
muestra1 = counts_df.columns[1]  # Primera muestra
muestra2 = counts_df.columns[2]  # Segunda muestra

plt.figure(figsize=(8, 6))
sns.scatterplot(
    data=counts_df,
    x=muestra1,
    y=muestra2,
    alpha=0.6,
    s=30
)
plt.title(f'Relación entre {muestra1} vs {muestra2}')
plt.xlabel(f'Conteos {muestra1}')
plt.ylabel(f'Conteos {muestra2}')
plt.xscale('log')
plt.yscale('log')
plt.savefig('figures/scatter_plot.png', dpi=300, bbox_inches='tight')
plt.show()
