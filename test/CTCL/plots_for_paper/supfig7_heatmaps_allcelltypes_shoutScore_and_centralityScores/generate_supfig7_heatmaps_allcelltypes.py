import pandas as pd
import os
import seaborn as sns
sns.set_theme(style="whitegrid")
import matplotlib.pyplot as plt
import numpy as np
from utilities import _calculate_mean_column_value_per_category_, _calculate_cell_results_by_sample_id_, _calculate_cell_results_by_patient_id_, _calculate_centrality_scores_by_patient_id_, _make_plot_name_from_column_                           

# ========== Load and pre-process data required for generating plots: ==========
cell_results=pd.read_csv(os.path.join('../../results', 'cell_results.csv'))
p_values_cell_type=pd.read_csv(os.path.join('../../results', 'p_values_cell_type.csv'))
squidpy_nhoodEnrichment_results=pd.read_csv(os.path.join('../../results', 'squidpy_nhoodEnrichment_results.csv'))
p_values_cell_type_squidpy_nhoodEnrichment=pd.read_csv(os.path.join('../../results', 'p_values_cell_type_squidpy_nhoodEnrichment.csv'))
squidpy_centralityScores_results=pd.read_csv(os.path.join('../../results', 'squidpy_centralityScores_results.csv'))
p_values_cell_type_squidpy_centralityScores=pd.read_csv(os.path.join('../../results', 'p_values_cell_type_squidpy_centralityScores.csv'))
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
patientID_abbreviatedID_lookup=_make_plot_name_from_column_(cell_results)

# abbreviated_patient_ids_with_condition_appended=_append_condition_
# ========== Create mosaic for plot: ==========
# layout=[
#         ["local_entropy_3", "local_homophily_3"],
#         ["egophily_3", "egophily_3"],
#         ["degree_centrality", "closeness_centrality"]
#         ]
# layout=[
#         ["local_entropy_3"],
#         ["local_homophily_3"],
#         ["egophily_3"],
#         ["degree_centrality"],
#         ["closeness_centrality"]
#         ]
layout = [
    ["local_entropy_3", "local_entropy_3", "local_homophily_3", "local_homophily_3", "egophily_3", "egophily_3"],
    ["degree_centrality", "degree_centrality", "degree_centrality", "closeness_centrality", "closeness_centrality", "closeness_centrality"]
]
fig, axes = plt.subplot_mosaic(layout, layout='constrained', figsize=(30,10))

# ========== Data for plots supfig7A1, A2, A3: ==========

# --- Calculate mean shout scores per condition: ---
mean_shout_scores=_calculate_mean_column_value_per_category_(conditions, cell_results, 'condition')

# --- Calculate interesting shout scores per celltype per condition: ---
cell_results_by_sample_id=_calculate_cell_results_by_sample_id_(celltypes, cell_results)
cell_results_by_patient_id=_calculate_cell_results_by_patient_id_(celltypes, cell_results)

# --- Generate heatmaps for shout (radius=r) scores: ---
scores = ['local_entropy_3', 'local_homophily_3', 'egophily_3']
for score in scores:
    mean_score_string='mean_'+score
    # DF=cell_results_by_patient_id.sort_values(by = mean_score_string).pivot(index="celltype", columns="patient_ids", values=mean_score_string)
    DF=cell_results_by_patient_id.pivot(index="celltype", columns="patient_ids", values=mean_score_string)
    DF=DF.sort_values('T-cells', axis=1)
    DF = pd.concat([DF.loc[['T-cells'],:], DF.drop('T-cells', axis=0)], axis=0)
    DF_with_modified_patientIDs=DF.copy()
    DF_with_modified_patientIDs=DF_with_modified_patientIDs.rename(columns=patientID_abbreviatedID_lookup)
    ax=sns.heatmap(DF_with_modified_patientIDs, mask=DF_with_modified_patientIDs.isna(), ax=axes[score], cbar_kws={"shrink": 0.5, "orientation": "horizontal", "pad":0.02},linewidths=0.8,linecolor="grey")
    ax.set_title(score)
    ax.set_xlabel('Patient id')
    ax.set_ylabel('Cell type')

# ========== Data for plots supfig7B1, B2, B3: ==========

# --- Calculate mean shout scores per condition: ---
mean_squidpy_centralityScores=_calculate_mean_column_value_per_category_(conditions, squidpy_centralityScores_results, 'condition')
# --- Calculate interesting centrality scores per celltype per condition: ---
centrality_scores_by_patient_id=_calculate_centrality_scores_by_patient_id_(celltypes, squidpy_centralityScores_results)

# --- Generate heatmaps for centrality scores: ---
scores = ['degree_centrality', 'closeness_centrality']
for score in scores:
    mean_score_string='mean_'+score
    # DF=centrality_scores_by_patient_id.sort_values(by = mean_score_string).pivot(index="celltype", columns="patient_ids", values=mean_score_string) # See 'Notes' at the end of this document for further explanation on the inbuilt 'pivot()' function.
    DF=centrality_scores_by_patient_id.pivot(index="celltype", columns="patient_ids", values=mean_score_string) # See 'Notes' at the end of this document for further explanation on the inbuilt 'pivot()' function.
    DF=DF.sort_values('T-cells', axis=1)
    DF = pd.concat([DF.loc[['T-cells'],:], DF.drop('T-cells', axis=0)], axis=0)
    DF_with_modified_patientIDs=DF.copy()
    DF_with_modified_patientIDs=DF_with_modified_patientIDs.rename(columns=patientID_abbreviatedID_lookup)
    ax=sns.heatmap(DF_with_modified_patientIDs, mask=DF_with_modified_patientIDs.isna(), ax=axes[score], cbar_kws={"shrink": 0.5, "orientation": "horizontal", "pad":0.02},linewidths=0.8,linecolor="grey")
    ax.set_title(score)
    ax.set_xlabel('Patient id')
    ax.set_ylabel('Cell type')

# # ========== Generate, save and show final plot (supfig7_heatmaps_allcelltypes): ==========
plt.savefig('supfig7_heatmaps_allcelltypes_shoutScore_and_centralityScores.pdf', format='pdf', bbox_inches='tight')
# fig.tight_layout()
plt.show()


# Notes:
# ======

# 1. pivot() function alternative code:
# DF=[]
# for celltype in celltypes:
#     df=cell_results_by_patient_id[cell_results_by_patient_id['celltype']==celltype][['patient_ids', mean_score_string]].T
#     df.columns = df.iloc[0]
#     df=df.drop(df.index[0])
#     df.index=[celltype]
#     DF.append(df)
# DF=pd.concat(DF, axis=0)
# DF = DF.astype(np.float32)

# 2. Sort dataframe rows by a particular column:
# dataframe.sort_values(by = 'column_name')

# Too many 



