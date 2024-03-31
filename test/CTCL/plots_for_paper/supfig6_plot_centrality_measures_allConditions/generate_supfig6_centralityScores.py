import pandas as pd
import os
import seaborn as sns
sns.set_theme(style="whitegrid")
import matplotlib.pyplot as plt
from _plot_centrality_measures_supfig6_ import _plot_centrality_measures_

# ========== Load and pre-process data required for generating plots: ==========
squidpy_centralityScores_results=pd.read_csv(os.path.join('../../results', 'squidpy_centralityScores_results.csv'))
p_values_cell_type_squidpy_centralityScores=pd.read_csv(os.path.join('../../results', 'p_values_cell_type_squidpy_centralityScores.csv'))
conditions_=['Eczema', 'T-Cell Lymphoma', 'Psoriasis']
conditions=['AD', 'PSO', 'CTCL']
conditions_abbreviations_dict=dict(zip(conditions_, conditions))
squidpy_centralityScores_results=squidpy_centralityScores_results.replace({"condition": conditions_abbreviations_dict})
p_values_cell_type_squidpy_centralityScores=p_values_cell_type_squidpy_centralityScores.replace({"condition_1": conditions_abbreviations_dict})
p_values_cell_type_squidpy_centralityScores=p_values_cell_type_squidpy_centralityScores.replace({"condition_2": conditions_abbreviations_dict})
celltypes=sorted(list(set(squidpy_centralityScores_results.celltypes)))

# ========== Create mosaic for plot: ==========
palette={"AD":(0.12156862745098039, 0.4666666666666667, 0.7058823529411765),
         "PSO": (0.17254901960784313, 0.6274509803921569, 0.17254901960784313),
         "CTCL": (1.0, 0.4980392156862745, 0.054901960784313725)}
layout = [
    ["B-cells", "B-cells", "Basal keratinocytes", "Basal keratinocytes", "Endothelial cells", "Endothelial cells"],
    ["Fibroblasts", "Fibroblasts", "Langerhans cells", "Langerhans cells", "Macrophages", "Macrophages"],
    ["Melanocytes", "Melanocytes", "Smooth muscle cells", "Smooth muscle cells", "Suprabasal keratinocytes", "Suprabasal keratinocytes"],
    [".", "T-cells", "T-cells", "Unknown", "Unknown", "."]
]

# ========== Data for plot supfig6: ==========
fig, axes = plt.subplot_mosaic(layout, figsize=(25,15))
celltype_count=0
for celltype in celltypes:
    celltype_count+=1
    if celltype_count==11:
        legend=True
    else:
        legend=None
    scores = ['degree_centrality', 'closeness_centrality']
    data=squidpy_centralityScores_results[squidpy_centralityScores_results['celltypes']==celltype]
    data_=[]
    for condition in conditions:
        data_.append(data[data['condition']==condition])
    data=pd.concat(data_, axis=0)
    columns=scores+['condition']
    data=data[columns]
    subplot_axis_id=axes[celltype]
    _plot_centrality_measures_(p_values_cell_type_squidpy_centralityScores, data, celltype, conditions, scores, subplot_axis_id, legend=legend, palette=palette)
# suptitle='Centrality scores for AD, CTCL and PSO samples'
# plt.suptitle(suptitle)
plt.subplots_adjust(wspace=0.60, hspace=0.50)
plt.savefig('supfig6.pdf', format='pdf', bbox_inches='tight')
# fig.tight_layout()
plt.show()







