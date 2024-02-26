import pandas as pd
import os
import seaborn as sns
sns.set_theme(style="whitegrid")
import matplotlib.pyplot as plt
from _plot_distributions_per_individual_celltypes_ import _plot_distributions_per_individual_celltypes_
from _plot_tcells_nhood_enrichment_ import _plot_tcells_nhood_enrichment_
from _plot_tcells_centrality_measures_ import _plot_tcells_centrality_measures_
from _generate_scatterplot_ import _generate_scatterplot_
from _generate_histplot_ import _generate_histplot_

# ========== Load and pre-process data required for generating plots: ==========
# if __name__ == '__main__':
cell_results=pd.read_csv(os.path.join('../results', 'cell_results.csv'))
p_values_cell_type=pd.read_csv(os.path.join('../results', 'p_values_cell_type.csv'))
squidpy_nhoodEnrichment_results=pd.read_csv(os.path.join('../results', 'squidpy_nhoodEnrichment_results.csv'))
p_values_cell_type_squidpy_nhoodEnrichment=pd.read_csv(os.path.join('../results', 'p_values_cell_type_squidpy_nhoodEnrichment.csv'))
squidpy_centralityScores_results=pd.read_csv(os.path.join('../results', 'squidpy_centralityScores_results.csv'))
p_values_cell_type_squidpy_centralityScores=pd.read_csv(os.path.join('../results', 'p_values_cell_type_squidpy_centralityScores.csv'))
conditions_=['Eczema', 'T-Cell Lymphoma', 'Psoriasis']
conditions=['AD', 'CTCL', 'PSO']
conditions_abbreviations_dict=dict(zip(conditions_, conditions))
cell_results=cell_results.replace({"condition": conditions_abbreviations_dict})
p_values_cell_type=p_values_cell_type.replace({"condition_1": conditions_abbreviations_dict})
p_values_cell_type=p_values_cell_type.replace({"condition_2": conditions_abbreviations_dict})
squidpy_nhoodEnrichment_results=squidpy_nhoodEnrichment_results.replace({"condition": conditions_abbreviations_dict})
p_values_cell_type_squidpy_nhoodEnrichment=p_values_cell_type_squidpy_nhoodEnrichment.replace({"condition_1": conditions_abbreviations_dict})
p_values_cell_type_squidpy_nhoodEnrichment=p_values_cell_type_squidpy_nhoodEnrichment.replace({"condition_2": conditions_abbreviations_dict})
squidpy_centralityScores_results=squidpy_centralityScores_results.replace({"condition": conditions_abbreviations_dict})
p_values_cell_type_squidpy_centralityScores=p_values_cell_type_squidpy_centralityScores.replace({"condition_1": conditions_abbreviations_dict})
p_values_cell_type_squidpy_centralityScores=p_values_cell_type_squidpy_centralityScores.replace({"condition_2": conditions_abbreviations_dict})

# ========== Create mosaic for plot: ==========

layout = [
    ["A", "A", "A", "A", "A", "A", "A", "A",    "C1", "C1", "C2", "C2", "C3", "C3"],
    ["A", "A", "A", "A", "A", "A", "A", "A",    "C1", "C1", "C2", "C2", "C3", "C3"],
    ["A", "A", "A", "A", "A", "A", "A", "A",    "C4", "C4", "C5", "C5", "C6", "C6"],
    ["A", "A", "A", "A", "A", "A", "A", "A",    "C4", "C4", "C5", "C5", "C6", "C6"],
    # -------------------------------------------------------------------------------
    ["B", "B", "B", "B", "B", "B", "B", "B",    "D", "D", "D",             "E", "E", "E"],
    ["B", "B", "B", "B", "B", "B", "B", "B",    "D", "D", "D",             "E", "E", "E"],
    ["B", "B", "B", "B", "B", "B", "B", "B",    "F1", "F1", "F1",       "F2", "F2", "F2"],
    ["B", "B", "B", "B", "B", "B", "B", "B",    "F1", "F1", "F1",       "F2", "F2", "F2"]
    # -------------------------------------------------------------------------------
]

fig, axes = plt.subplot_mosaic(layout, layout='constrained', figsize=(20,10))

# ========== Data for plot A: ==========
heterogeneity_measure='entropy'
celltype='T-cells'
conditions=['AD', 'CTCL']
radii=[1, 2, 3, 4, 5]
subplot_axis_id=axes["A"]
legend_loc='lower right'
title_prefix='A. '
# Generate plot-1:
_plot_distributions_per_individual_celltypes_(p_values_cell_type, cell_results, celltype, heterogeneity_measure, conditions, radii, subplot_axis_id, legend_loc, title_prefix)

# ========== Data for plot B: ==========
heterogeneity_measure='egophily'
celltype='T-cells'
conditions=['AD', 'CTCL']
radii=[1, 2, 3, 4, 5]
subplot_axis_id=axes["B"]
legend_loc='upper right'
title_prefix='C. '
# Generate plot-B:
_plot_distributions_per_individual_celltypes_(p_values_cell_type, cell_results, celltype, heterogeneity_measure, conditions, radii, subplot_axis_id, legend_loc, title_prefix)


# ========== Data for plot C: ==========
celltypes_=['Smooth muscle cells', 'Macrophages', 'Endothelial cells', 'B-cells', 'Fibroblasts', 'Basal keratinocytes']
conditions=['AD', 'CTCL']
radii=[1, 2, 3, 4, 5]
cnt=0
for celltype in celltypes_:
    cnt+=1
    if cnt==1:
        title_prefix='B. '
        suptitle='Neighborhood enrichment (T-cells)'
    else:
        title_prefix=None
        suptitle=None
    data=pd.concat([squidpy_nhoodEnrichment_results[squidpy_nhoodEnrichment_results['celltype_1']==celltype][squidpy_nhoodEnrichment_results['celltype_2']=='T-cells'], squidpy_nhoodEnrichment_results[squidpy_nhoodEnrichment_results['celltype_2']==celltype][squidpy_nhoodEnrichment_results['celltype_1']=='T-cells']], axis=0)
    subplot_axis_id=axes[f"C{cnt}"]
    data_=[]
    for condition in conditions:
        data_.append(data[data['condition']==condition])
    data=pd.concat(data_, axis=0)
    _plot_tcells_nhood_enrichment_(p_values_cell_type_squidpy_nhoodEnrichment, squidpy_nhoodEnrichment_results, data, celltype, conditions, radii, subplot_axis_id, title_prefix, suptitle)

# ========== Data for plot D: ==========
celltype='T-cells'
conditions=['AD', 'CTCL']
radii=[1, 2, 3, 4, 5]
title_prefix='D. '
scores = ['degree_centrality', 'closeness_centrality']
data=squidpy_centralityScores_results[squidpy_centralityScores_results['celltypes']==celltype]
data_=[]
for condition in conditions:
    data_.append(data[data['condition']==condition])
data=pd.concat(data_, axis=0)
columns=scores+['condition']
data=data[columns]
subplot_axis_id=axes["D"]
_plot_tcells_centrality_measures_(p_values_cell_type_squidpy_centralityScores, data, celltype, conditions, radii, scores, subplot_axis_id, title_prefix)

# ========== Data for plot E: ==========
celltype='T-cells'
conditions=['AD', 'CTCL']
x_axis_param='local_entropy_3'
y_axis_param='egophily_3'
hue='condition'
title_prefix='E. '
data=cell_results[cell_results['cell_type']==celltype]
data_=[]
for condition in conditions:
    data_.append(data[data['condition']==condition])
data=pd.concat(data_, axis=0)
subplot_axis_id=axes["E"]
legend_loc="upper left"
_generate_scatterplot_(x_axis_param, y_axis_param, data, subplot_axis_id, hue, legend_loc, title_prefix)


# ========== Data for plot F: ==========
celltype='T-cells'
conditions=['AD', 'CTCL']
scores=['local_entropy_3', 'egophily_3']
hue='condition'
legend_locs=['upper left', 'upper right']
data=cell_results[cell_results['cell_type']==celltype]

cnt=0
for score in scores:
    cnt+=1
    if cnt==1:
        title=r"$\bf{"+'F.'+"}$"+f'Distributions of {score} scores\nfor AD and CTCL samples'
    elif cnt==2:
        title_prefix=r"$\bf{"+'G.'+"}$"+f'Distributions of {score} scores\nfor AD and CTCL samples'
    subplot_axis_id=axes[f"F{cnt}"]
    data_=[]
    for condition in conditions:
        data_.append(data[data['condition']==condition])
    data=pd.concat(data_, axis=0)
    legend_index=cnt-1
    legend_loc=legend_locs[legend_index]
    _generate_histplot_(score, data, subplot_axis_id, hue, legend_loc, title)

# ========== Generate, save and show final plot (fig1): ==========
plt.savefig('fig1.pdf', format='pdf', bbox_inches='tight')
fig.tight_layout()
plt.show()

