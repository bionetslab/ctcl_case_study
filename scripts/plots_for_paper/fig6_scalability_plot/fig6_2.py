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
# fig, axes = plt.subplot_mosaic(layout, figsize=(25,12.5), sharey=True)
fig, axes = plt.subplot_mosaic(layout, figsize=(25,6.5))

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
samples=[]
for filename in os.listdir('../../../data'):
    filename = os.fsdecode(filename)
    if filename.endswith('.h5ad'):
        samples.append(filename.split('.')[0])

title_prefix="A"
plot_title=None
subplot_axis_id=axes['all_samples']
sample_results=pd.read_csv(os.path.join('../../../results', 'sample_results.csv'))
no_of_cells_dict=dict(zip(list(sample_results['sample_id']), list(sample_results['number_of_cells'])))
data_all_cells_radius_equals_5=data_all_cells[data_all_cells['variable']==5]
sample_ids=list(data_all_cells_radius_equals_5['sample'])
runtimes=list(data_all_cells_radius_equals_5['value'])
no_of_cells_list=[]
for i in samples:
    no_of_cells_list.append(no_of_cells_dict[int(i)])
df=pd.DataFrame()
df['number_of_cells']=no_of_cells_list
df['SHouT_execution_time']=runtimes
ax=sns.regplot(data=df, x='number_of_cells', y='SHouT_execution_time', ax=subplot_axis_id, scatter_kws={'s':60}, color=(0.8666666666666667, 0.5176470588235295, 0.3215686274509804))
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
ax.set_title(plot_title, loc=title_loc, fontsize=30) # fontsize=25, pad=20
ax.set_ylabel('Runtime, radius $r=5$ (sec)', fontsize=25, labelpad=5) # fontsize=25, labelpad=20
ax.set_xlabel('Cell count', fontsize=25, labelpad=5)
ax.yaxis.set_tick_params(which='both', labelleft=True)
# ========== Generate, save and show final plot (fig1): ==========
plt.savefig('fig6_scalability_plot.pdf', format='pdf', bbox_inches='tight')
plt.show()

