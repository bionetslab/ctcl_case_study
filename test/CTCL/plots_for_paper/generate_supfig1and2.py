import pandas as pd
import os
import seaborn as sns
sns.set_theme(style="whitegrid")
import matplotlib.pyplot as plt
from _plot_distributions_per_individual_celltypes_ import _plot_distributions_per_individual_celltypes_
from _plot_tcells_nhood_enrichment_ import _plot_tcells_nhood_enrichment_
from _plot_tcells_centrality_measures_ import _plot_tcells_centrality_measures_
from _generate_scatterplot_ import _generate_scatterplot_
from _generate_histplot_supfig1and2_ import _generate_histplot_

# ========== Load and pre-process data required for generating plots: ==========
cell_results=pd.read_csv(os.path.join('../results', 'cell_results.csv'))
p_values_cell_type=pd.read_csv(os.path.join('../results', 'p_values_cell_type.csv'))
conditions_=['Eczema', 'T-Cell Lymphoma', 'Psoriasis']
conditions=['AD', 'CTCL', 'PSO']
conditions_abbreviations_dict=dict(zip(conditions_, conditions))
cell_results=cell_results.replace({"condition": conditions_abbreviations_dict})
p_values_cell_type=p_values_cell_type.replace({"condition_1": conditions_abbreviations_dict})
p_values_cell_type=p_values_cell_type.replace({"condition_2": conditions_abbreviations_dict})
celltypes=list(set(cell_results.cell_type))
celltypes=[i for i in celltypes if i != 'T-cells']
celltypes.append('T-cells')
# ========== Create mosaic for plot: ==========


layout = [
    ["Basal keratinocytes", "Smooth muscle cells", "Melanocytes"],
    ["Unknown", "B-cells", "Fibroblasts"],
    [ "Suprabasal keratinocytes", "Endothelial cells", "Macrophages"],
    ["Langerhans cells", "T-cells", ""]
]





# ========== Data for plot supfig1: ==========
fig, axes = plt.subplot_mosaic(layout, layout='constrained', figsize=(20,10))
subplot_str='A'
for celltype in celltypes:
    conditions=['AD', 'CTCL']
    score='local_entropy_3'
    hue='condition'
    data=cell_results[cell_results['cell_type']==celltype]
    title=subplot_str+f'. {celltype}'
    subplot_axis_id=axes[celltype]
    data_=[]
    for condition in conditions:
        data_.append(data[data['condition']==condition])
    data=pd.concat(data_, axis=0)
    legend_loc="upper left"
    _generate_histplot_(score, data, subplot_axis_id, hue, legend_loc, title)
    subplot_str=chr(ord(subplot_str)+1)
# ========== Generate, save and show final plot (supfig1): ==========
suptitle=f'Distributions of {score} scores for AD and CTCL samples'
plt.suptitle(suptitle)
plt.savefig('supfig1.pdf', format='pdf', bbox_inches='tight')
fig.tight_layout()
plt.show()





# ========== Data for plot supfig2: ==========
fig, axes = plt.subplot_mosaic(layout, layout='constrained', figsize=(20,10))
subplot_str='A'
for celltype in celltypes:
    conditions=['AD', 'CTCL']
    score='egophily_3'
    hue='condition'
    data=cell_results[cell_results['cell_type']==celltype]
    title=subplot_str+f'. {celltype}'
    subplot_axis_id=axes[celltype]
    data_=[]
    for condition in conditions:
        data_.append(data[data['condition']==condition])
    data=pd.concat(data_, axis=0)
    legend_loc="upper right"
    _generate_histplot_(score, data, subplot_axis_id, hue, legend_loc, title)
    subplot_str=chr(ord(subplot_str)+1)
# ========== Generate, save and show final plot (supfig2): ==========
suptitle=f'Distributions of {score} scores for AD and CTCL samples'
plt.suptitle(suptitle)
plt.savefig('supfig2.pdf', format='pdf', bbox_inches='tight')
fig.tight_layout()
plt.show()

