import pandas as pd
import os
import seaborn as sns
sns.set_theme(style="whitegrid")
import matplotlib.pyplot as plt
from _generate_boxplot_ import _generate_boxplot_
sns.set_theme(style="whitegrid")
radii=[1, 5, 10, 20, 40, 80, 100]

layout = [
    ["all_samples", "scatter"]
]
fig, axes = plt.subplot_mosaic(layout, figsize=(25,12.5))

# ===== I. "all_cells": =====
data_all_cells=pd.read_csv('../../../results/data_all_cells.csv')
var_name='variable'
val_name='value'
   
# Generate cell count plots for all celltypes:
title_prefix='B'
plot_title=None
# Generate violinplot:
subplot_axis_id=axes["scatter"]
legend_loc='upper right'
title_loc='left'
_generate_boxplot_(data_all_cells, var_name, val_name, subplot_axis_id, radii, legend_loc=legend_loc, title_prefix=title_prefix, plot_title=plot_title, ylabel='radius', title_loc=title_loc)

# ===== II. Generate scatter plot for "scatter": =====
title_prefix="A"
plot_title=None
subplot_axis_id=axes['all_samples']
sample_results=pd.read_csv(os.path.join('../../../results', 'sample_results.csv'))
ax=sns.regplot(data=sample_results, x='number_of_cells', y='SHouT_execution_time', ax=subplot_axis_id, scatter_kws={'s':60})
yticks=[]
for j in ax.get_yticks():
    yticks.append(round(j,1))
ax.set_yticklabels(yticks, size = 20) # size = 20
ax.tick_params(axis='x', which='major', labelsize=20)
ax.grid(False)
if title_prefix:
    if plot_title:
        plot_title=r"$\bf{" + title_prefix + "}$" + plot_title
    else:
        plot_title=r"$\bf{" + title_prefix + "}$"
ax.set_title(plot_title, loc=title_loc, fontsize=25) # fontsize=25, pad=20
ax.set_ylabel('Runtime, radius $r=5$ (sec)', fontsize=25, labelpad=5) # fontsize=25, labelpad=20
ax.set_xlabel('Cell count', fontsize=25, labelpad=5)

# ========== Generate, save and show final plot (fig1): ==========
plt.savefig('fig6_scalability_plot.pdf', format='pdf', bbox_inches='tight')
plt.show()

