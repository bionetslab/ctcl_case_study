import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
sns.set_theme(style="whitegrid")

# ========== Load and pre-process data required for generating plots: ==========
sample_results=pd.read_csv(os.path.join('../results', 'sample_results.csv'))
fig, axes = plt.subplots(figsize=(20,10))
sns.scatterplot(data=sample_results, x='number_of_cells', y='SHouT_execution_time', ax=axes)

# ========== Generate, save and show final plot (fig1): ==========
# plt.subplots_adjust(wspace=0.02, hspace=0.5)
plt.subplots_adjust(wspace=0.45, hspace=1.8)
plt.savefig('fig_scalability_varying_across_no_of_cells.pdf', format='pdf', bbox_inches='tight')
# fig.tight_layout(pad=100.0)
plt.show()