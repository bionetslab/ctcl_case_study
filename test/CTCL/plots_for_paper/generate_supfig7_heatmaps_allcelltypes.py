import pandas as pd
import os
import seaborn as sns
sns.set_theme(style="whitegrid")
import matplotlib.pyplot as plt
import numpy as np
from _plot_distributions_per_individual_celltypes_ import _plot_distributions_per_individual_celltypes_
from _plot_tcells_nhood_enrichment_ import _plot_tcells_nhood_enrichment_
from _plot_tcells_centrality_measures_ import _plot_tcells_centrality_measures_
from _generate_scatterplot_ import _generate_scatterplot_
from _generate_histplot_ import _generate_histplot_
from utilities import _calculate_mean_column_value_per_category_, _calculate_cell_results_by_sample_id_, _calculate_cell_results_by_patient_id_

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
celltypes=list(set(cell_results.cell_type))

# ========== Create mosaic for plot: ==========
layout = [
    ["local_entropy_3", "local_entropy_3", "local_homophily_3", "local_homophily_3", "egophily_3", "egophily_3"],
    ["0", "degree_centrality", "degree_centrality", "closeness_centrality", "closeness_centrality", "1"]
]
fig, axes = plt.subplot_mosaic(layout, layout='constrained', figsize=(20,10))

# fig, axes = plt.subplots(figsize=(20,10))

# ========== Data for plots supfig7A1, A2, A3, B1, B2, B3: ==========

# --- Calculate mean shout scores per condition: ---
mean_shout_scores=_calculate_mean_column_value_per_category_(conditions, cell_results, 'condition')

# --- Calculate interesting shout scores per celltype per condition: ---
cell_results_by_sample_id=_calculate_cell_results_by_sample_id_(celltypes, cell_results)
cell_results_by_patient_id=_calculate_cell_results_by_patient_id_(celltypes, cell_results)

# --- Generate heatmaps for shout (radius=r) scores: ---
scores = ['mean_local_entropy_3', 'mean_local_homophily_3', 'mean_egophily_3']
for score in scores:
    DF=[]
    for celltype in celltypes:
        df=cell_results_by_patient_id[cell_results_by_patient_id['celltype']==celltype][['patient_ids', score]].T
        df.columns = df.iloc[0]
        df=df.drop(df.index[0])
        df.index=[celltype]
        DF.append(df)
    DF=pd.concat(DF, axis=0)
    DF = DF.astype(np.float32)
    
    # sns.heatmap(DF, mask=DF.isnull())
    ax=sns.heatmap(ax=, DF, mask=DF.isna())
    # plt.savefig('fig1.pdf', format='pdf', bbox_inches='tight')
    # fig.tight_layout()
    # plt.show()
    
    # # OR:
    # DF=cell_results_by_patient_id.pivot(index="celltype", columns="patient_id", values=score)

# # df_by_patient_ids.sort_values(by = 'mean_local_entropy_3')

# # ========== Data for plot A: ==========
# heterogeneity_measure='entropy'
# celltype='T-cells'
# conditions=['AD', 'CTCL']
# radii=[1, 2, 3, 4, 5]
# subplot_axis_id=axes["A"]
# legend_loc='lower right'
# title_prefix='A. '
# # Generate plot-1:
# _plot_distributions_per_individual_celltypes_(p_values_cell_type, cell_results, celltype, heterogeneity_measure, conditions, radii, subplot_axis_id, legend_loc, title_prefix)

# # ========== Data for plot B: ==========
# heterogeneity_measure='egophily'
# celltype='T-cells'
# conditions=['AD', 'CTCL']
# radii=[1, 2, 3, 4, 5]
# subplot_axis_id=axes["B"]
# legend_loc='upper right'
# title_prefix='C. '
# # Generate plot-B:
# _plot_distributions_per_individual_celltypes_(p_values_cell_type, cell_results, celltype, heterogeneity_measure, conditions, radii, subplot_axis_id, legend_loc, title_prefix)


# # ========== Data for plot C: ==========
# celltypes_=['Smooth muscle cells', 'Macrophages', 'Endothelial cells', 'B-cells', 'Fibroblasts', 'Basal keratinocytes']
# conditions=['AD', 'CTCL']
# radii=[1, 2, 3, 4, 5]
# cnt=0
# for celltype in celltypes_:
#     cnt+=1
#     if cnt==1:
#         title_prefix='B. '
#         suptitle='Neighborhood enrichment (T-cells)'
#     else:
#         title_prefix=None
#         suptitle=None
#     data=pd.concat([squidpy_nhoodEnrichment_results[squidpy_nhoodEnrichment_results['celltype_1']==celltype][squidpy_nhoodEnrichment_results['celltype_2']=='T-cells'], squidpy_nhoodEnrichment_results[squidpy_nhoodEnrichment_results['celltype_2']==celltype][squidpy_nhoodEnrichment_results['celltype_1']=='T-cells']], axis=0)
#     subplot_axis_id=axes[f"C{cnt}"]
#     data_=[]
#     for condition in conditions:
#         data_.append(data[data['condition']==condition])
#     data=pd.concat(data_, axis=0)
#     _plot_tcells_nhood_enrichment_(p_values_cell_type_squidpy_nhoodEnrichment, squidpy_nhoodEnrichment_results, data, celltype, conditions, radii, subplot_axis_id, title_prefix, suptitle)

# # ========== Data for plot D: ==========
# celltype='T-cells'
# conditions=['AD', 'CTCL']
# radii=[1, 2, 3, 4, 5]
# title_prefix='D. '
# scores = ['degree_centrality', 'closeness_centrality']
# data=squidpy_centralityScores_results[squidpy_centralityScores_results['celltypes']==celltype]
# data_=[]
# for condition in conditions:
#     data_.append(data[data['condition']==condition])
# data=pd.concat(data_, axis=0)
# columns=scores+['condition']
# data=data[columns]
# subplot_axis_id=axes["D"]
# _plot_tcells_centrality_measures_(p_values_cell_type_squidpy_centralityScores, data, celltype, conditions, radii, scores, subplot_axis_id, title_prefix)

# # ========== Data for plot E: ==========
# celltype='T-cells'
# conditions=['AD', 'CTCL']
# x_axis_param='local_entropy_3'
# y_axis_param='egophily_3'
# hue='condition'
# title_prefix='E. '
# data=cell_results[cell_results['cell_type']==celltype]
# data_=[]
# for condition in conditions:
#     data_.append(data[data['condition']==condition])
# data=pd.concat(data_, axis=0)
# subplot_axis_id=axes["E"]
# legend_loc="upper left"
# _generate_scatterplot_(x_axis_param, y_axis_param, data, subplot_axis_id, hue, legend_loc, title_prefix)


# # ========== Data for plot F: ==========
# celltype='T-cells'
# conditions=['AD', 'CTCL']
# scores=['local_entropy_3', 'egophily_3']
# hue='condition'
# legend_locs=['upper left', 'upper right']
# data=cell_results[cell_results['cell_type']==celltype]

# cnt=0
# for score in scores:
#     cnt+=1
#     if cnt==1:
#         title_prefix='F. '
#     elif cnt==2:
#         title_prefix='G. '
#     subplot_axis_id=axes[f"F{cnt}"]
#     data_=[]
#     for condition in conditions:
#         data_.append(data[data['condition']==condition])
#     data=pd.concat(data_, axis=0)
#     legend_index=cnt-1
#     legend_loc=legend_locs[legend_index]
#     _generate_histplot_(score, data, subplot_axis_id, hue, legend_loc, title_prefix)

# # ========== Generate, save and show final plot (fig1): ==========
# plt.savefig('fig1.pdf', format='pdf', bbox_inches='tight')
# fig.tight_layout()
# plt.show()

