import pandas as pd
import os
import seaborn as sns
sns.set_theme(style="whitegrid")
import matplotlib.pyplot as plt
from _plot_tcells_centrality_measures_supfig6_ import _plot_tcells_centrality_measures_

# ========== Load and pre-process data required for generating plots: ==========
squidpy_centralityScores_results=pd.read_csv(os.path.join('../results', 'squidpy_centralityScores_results.csv'))
p_values_cell_type_squidpy_centralityScores=pd.read_csv(os.path.join('../results', 'p_values_cell_type_squidpy_centralityScores.csv'))
conditions_=['Eczema', 'T-Cell Lymphoma', 'Psoriasis']
conditions=['AD', 'CTCL', 'PSO']
conditions_abbreviations_dict=dict(zip(conditions_, conditions))
squidpy_centralityScores_results=squidpy_centralityScores_results.replace({"condition": conditions_abbreviations_dict})
p_values_cell_type_squidpy_centralityScores=p_values_cell_type_squidpy_centralityScores.replace({"condition_1": conditions_abbreviations_dict})
p_values_cell_type_squidpy_centralityScores=p_values_cell_type_squidpy_centralityScores.replace({"condition_2": conditions_abbreviations_dict})
celltypes=squidpy_centralityScores_results.celltypes

# ========== Create mosaic for plot: ==========

celltypes=["Suprabasal keratinocytes", "Smooth muscle cells",
"Macrophages", "Fibroblasts", 
"T-cells", "Endothelial cells"]

layout = [
    ["A", "B", "C"],
    ["D", "E", "F"],
]

# ========== Data for plot supfig6: ==========
fig, axes = plt.subplot_mosaic(layout, layout='constrained', figsize=(20,10))
subplot_str='A'
for celltype in celltypes:
    # conditions=['AD', 'PSO', 'CTCL']
    conditions=['AD', 'CTCL']
    scores = ['degree_centrality', 'closeness_centrality']
    data=squidpy_centralityScores_results[squidpy_centralityScores_results['celltypes']==celltype]
    data_=[]
    for condition in conditions:
        data_.append(data[data['condition']==condition])
    data=pd.concat(data_, axis=0)
    columns=scores+['condition']
    data=data[columns]
    subplot_axis_id=axes[subplot_str]
    _plot_tcells_centrality_measures_(p_values_cell_type_squidpy_centralityScores, data, celltype, conditions, scores, subplot_axis_id)
    subplot_str=chr(ord(subplot_str)+1)
suptitle='Centrality scores for AD and CTCL samples'
plt.suptitle(suptitle)
plt.savefig('supfig6.pdf', format='pdf', bbox_inches='tight')
fig.tight_layout()
plt.show()







