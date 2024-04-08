import pandas as pd
import os
import seaborn as sns
sns.set_theme(style="whitegrid")
import matplotlib.pyplot as plt
import scanpy as sc
import time
from _generate_boxplot_ import _generate_boxplot_
import shout
sns.set_theme(style="whitegrid")
radii=[1, 5, 10, 20, 40, 80, 100]

layout = [
    ["all_samples", "scatter", ],
    ["315_few_cells", "309_many_cells"],
]
fig, axes = plt.subplot_mosaic(layout, figsize=(25,25))

# ===== I. "all_cells": =====
# times=[]
# for filename in os.listdir('test/CTCL/data'):
#     filename = os.fsdecode(filename)
#     if filename.endswith('.h5ad'):
#         print(filename)
#         adata=sc.read_h5ad('test/CTCL/data/'+filename)
#         times_=pd.DataFrame()
#         for radius in radii:
#             start = time.time()
#             shout.all_scores(adata, cluster_key='celltype', radii=[radius])
#             end = time.time()
#             times_[radius]=[end-start]
#         times.append(times_)
# times=pd.concat(times, axis=0)
# data_all_cells=pd.melt(times)
# data_all_cells.to_csv('data_all_cells.csv')
data_all_cells=pd.read_csv('data_all_cells.csv')
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
sample_results=pd.read_csv(os.path.join('test/CTCL/results', 'sample_results.csv'))
# ax=sns.scatterplot(data=sample_results, x='number_of_cells', y='SHouT_execution_time', ax=subplot_axis_id, s=10)
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
# ax.set_xlabel(xlabel, fontsize=25) # fontsize=25, labelpad=20
ax.set_ylabel('Runtime, radius $r=5$ (sec)', fontsize=25, labelpad=5) # fontsize=25, labelpad=20
ax.set_xlabel('Cell count', fontsize=25, labelpad=5)

# ===== III. "315_few_cells": =====
title_prefix="C"
plot_title=None
subplot_axis_id=axes["315_few_cells"]
# runtimes={}
# n=[1, 5, 10, 20, 40, 80, 100]
# adata=sc.read_h5ad('test/CTCL/data/315.h5ad')

# for _count_ in range (100):
#     for i in n:
#         print('Iteration '+str(_count_+1)+', i= '+str(i))
#         # ---
#         start = time.time()
#         scores_mibitof = shout.all_scores(adata, cluster_key='celltype', radii=[i])
#         end = time.time()
#         try:
#             runtimes[i].append(end-start)
#         except:
#             runtimes[i]=[]
#             runtimes[i].append(end-start)
# runtimes_df=pd.DataFrame.from_dict(runtimes)
# sns.set_theme(style="whitegrid", palette="muted")
# df_melted_315=pd.melt(runtimes_df)
# df_melted_315.to_csv('df_melted_315.csv')
df_melted_315=pd.read_csv('df_melted_315.csv')
ax=sns.boxplot(x = 'variable', y = 'value', data = df_melted_315, ax=subplot_axis_id, orient='v')
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
# ax.set_xlabel(xlabel, fontsize=25) # fontsize=25, labelpad=20
ax.set_ylabel('Runtime, samples with few cells (sec)', fontsize=25, labelpad=5) # fontsize=25, labelpad=20
ax.set_xlabel('Radius ($r$)', fontsize=25, labelpad=5)

# ===== I. "309_many_cells": =====
title_prefix="D"
plot_title=None
subplot_axis_id=axes["309_many_cells"]
# runtimes={}
# n=[1, 5, 10, 20, 40, 80, 100]
# adata=sc.read_h5ad('test/CTCL/data/309.h5ad')

# for _count_ in range (100):
#     for i in n:
#         print('Iteration '+str(_count_+1)+', i= '+str(i))
#         start = time.time()
#         scores_mibitof = shout.all_scores(adata, cluster_key='celltype', radii=[i])
#         end = time.time()
#         try:
#             runtimes[i].append(end-start)
#         except:
#             runtimes[i]=[]
#             runtimes[i].append(end-start)
# runtimes_df=pd.DataFrame.from_dict(runtimes)
# sns.set_theme(style="whitegrid", palette="muted")
# df_melted_309=pd.melt(runtimes_df)
# df_melted_309.to_csv('df_melted_309.csv')
df_melted_309=pd.read_csv('df_melted_309.csv')
ax=sns.boxplot(x = 'variable', y = 'value', data = df_melted_309, ax=subplot_axis_id, orient='v')
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
# ax.set_xlabel(xlabel, fontsize=25) # fontsize=25, labelpad=20
ax.set_ylabel('Runtime, samples with many cells (sec)', fontsize=25, labelpad=5) # fontsize=25, labelpad=20
ax.set_xlabel('Radius ($r$)', fontsize=25, labelpad=5)

# ========== Generate, save and show final plot (fig1): ==========
plt.savefig('scalability_plot.pdf', format='pdf', bbox_inches='tight')
plt.show()

