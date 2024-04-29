import pandas as pd
import os
import seaborn as sns
sns.set_theme(style="whitegrid")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

# ========== Load and pre-process data required for generating plots: ==========
radii=[5]
conditions_=['Eczema', 'T-Cell Lymphoma', 'Psoriasis']
conditions=['AD', 'CTCL', 'PSO']
conditions_abbreviations_dict=dict(zip(conditions_, conditions))
celltypes_abbreviations_dict={"B-cells":"B-cells",
                              "Basal keratinocytes": "Basal kers",
                              "Endothelial cells": "Endothelial",
                              "Fibroblasts": "Fibroblasts",
                              "Langerhans cells": "Langerhans",
                              "Macrophages": "Macrophages",
                              "Melanocytes": "Melanocytes",
                              "Smooth muscle cells": "SMCs",
                              "Suprabasal keratinocytes": "Suprabasal kers",
                              "T-cells": "T-cells",
                              "Unknown": "Unknown"}

p_values_cell_type=pd.read_csv(os.path.join('../../../results', 'p_values_cell_type.csv'))
p_values_cell_type=p_values_cell_type.replace({"condition_1": conditions_abbreviations_dict})
p_values_cell_type=p_values_cell_type.replace({"condition_2": conditions_abbreviations_dict})
p_values_cell_type_subsampled_patientwise=pd.read_csv(os.path.join('../../../results', 'p_values_cell_type_subsampled_patientwise.csv'))
p_values_cell_type_subsampled_patientwise=p_values_cell_type_subsampled_patientwise.replace({"condition_1": conditions_abbreviations_dict})
p_values_cell_type_subsampled_patientwise=p_values_cell_type_subsampled_patientwise.replace({"condition_2": conditions_abbreviations_dict})

celltypes_of_interest=["T-cells", "T-cells", "T-cells", "T-cells", "T-cells", "T-cells", "Basal keratinocytes", "Basal keratinocytes", "Basal keratinocytes", "Basal keratinocytes", "Basal keratinocytes", "Basal keratinocytes"]
condition_pairs_of_interest=[('CTCL', 'AD'), ('CTCL', 'PSO'), ('AD', 'PSO'), ('CTCL', 'AD'), ('CTCL', 'PSO'), ('AD', 'PSO'), ('CTCL', 'AD'), ('CTCL', 'PSO'), ('AD', 'PSO'), ('CTCL', 'AD'), ('CTCL', 'PSO'), ('AD', 'PSO')]
scores_of_interest=["entropy", "entropy", "entropy", "egophily", "egophily", "egophily", "entropy", "entropy", "entropy", "homophily", "homophily", "homophily"]

layout = [
    ["T-cells_entropy_CTCLvsAD", "T-cells_entropy_CTCLvsAD", "T-cells_entropy_CTCLvsPSO", "T-cells_entropy_CTCLvsPSO"], # "T-cells_egophily_CTCLvsAD", "T-cells_egophily_CTCLvsAD"],
    ["T-cells_entropy_ADvsPSO", "T-cells_entropy_ADvsPSO", "T-cells_egophily_CTCLvsAD", "T-cells_egophily_CTCLvsAD",],
    ["T-cells_egophily_CTCLvsPSO", "T-cells_egophily_CTCLvsPSO", "T-cells_egophily_ADvsPSO", "T-cells_egophily_ADvsPSO",],
    # --------
    ["Basal keratinocytes_entropy_CTCLvsAD", "Basal keratinocytes_entropy_CTCLvsAD", "Basal keratinocytes_entropy_CTCLvsPSO", "Basal keratinocytes_entropy_CTCLvsPSO"],
    ["Basal keratinocytes_entropy_ADvsPSO", "Basal keratinocytes_entropy_ADvsPSO", "Basal keratinocytes_homophily_CTCLvsAD", "Basal keratinocytes_homophily_CTCLvsAD"],
    ["Basal keratinocytes_homophily_CTCLvsPSO", "Basal keratinocytes_homophily_CTCLvsPSO", "Basal keratinocytes_homophily_ADvsPSO", "Basal keratinocytes_homophily_ADvsPSO"]
]
fig, axes = plt.subplot_mosaic(layout, figsize=(25,25))
for count in range(12):
    if scores_of_interest[count]=='egophily':
        local_heterogeneity_measure=[f'{scores_of_interest[count]}_{radius}' for radius in radii]
        local_heterogeneity_score=f'{scores_of_interest[count].capitalize()}'
    else:
        local_heterogeneity_measure=[f'local_{scores_of_interest[count]}_{radius}' for radius in radii]
        local_heterogeneity_score=f'Local {scores_of_interest[count]}'
    if celltypes_of_interest[count]=='T-cells':
        title=f'{local_heterogeneity_score}, {celltypes_of_interest[count]} ({condition_pairs_of_interest[count][0]} vs. {condition_pairs_of_interest[count][1]})'
    elif celltypes_of_interest[count]=='Basal keratinocytes':
        title=f'{local_heterogeneity_score}, {celltypes_of_interest[count].lower()} ({condition_pairs_of_interest[count][0]} vs. {condition_pairs_of_interest[count][1]})'
    xlabel="$-log_{10}$"+"($P_{subsampled}$)"
    data_disease_combinations=[]
    for k in local_heterogeneity_measure:
        if len(p_values_cell_type[p_values_cell_type['cell_type']==celltypes_of_interest[count]][p_values_cell_type['condition_1']==condition_pairs_of_interest[count][0]][p_values_cell_type['condition_2']==condition_pairs_of_interest[count][1]][p_values_cell_type['score']==k]['p_value'])==0:
            p=float(p_values_cell_type[p_values_cell_type['cell_type']==celltypes_of_interest[count]][p_values_cell_type['condition_1']==condition_pairs_of_interest[count][1]][p_values_cell_type['condition_2']==condition_pairs_of_interest[count][0]][p_values_cell_type['score']==k]['p_value'].values[0])
            p=-np.log10(p)
        else:
            p=float(p_values_cell_type[p_values_cell_type['cell_type']==celltypes_of_interest[count]][p_values_cell_type['condition_1']==condition_pairs_of_interest[count][0]][p_values_cell_type['condition_2']==condition_pairs_of_interest[count][1]][p_values_cell_type['score']==k]['p_value'].values[0])
            p=-np.log10(p)
        if len(p_values_cell_type_subsampled_patientwise[p_values_cell_type_subsampled_patientwise['cell_type']==celltypes_of_interest[count]][p_values_cell_type_subsampled_patientwise['condition_1']==condition_pairs_of_interest[count][0]][p_values_cell_type_subsampled_patientwise['condition_2']==condition_pairs_of_interest[count][1]][p_values_cell_type_subsampled_patientwise['score']==k]['p_value'])==0:
            data_=p_values_cell_type_subsampled_patientwise[p_values_cell_type_subsampled_patientwise['cell_type']==celltypes_of_interest[count]][p_values_cell_type_subsampled_patientwise['condition_1']==condition_pairs_of_interest[count][1]][p_values_cell_type_subsampled_patientwise['condition_2']==condition_pairs_of_interest[count][0]][p_values_cell_type_subsampled_patientwise['score']==k]
        else:
            data_=p_values_cell_type_subsampled_patientwise[p_values_cell_type_subsampled_patientwise['cell_type']==celltypes_of_interest[count]][p_values_cell_type_subsampled_patientwise['condition_1']==condition_pairs_of_interest[count][0]][p_values_cell_type_subsampled_patientwise['condition_2']==condition_pairs_of_interest[count][1]][p_values_cell_type_subsampled_patientwise['score']==k]
        data_disease_combinations.append(data_)
    data_disease_combinations=pd.concat(data_disease_combinations, axis=0)
    df=data_disease_combinations.copy()
    df=df[['score', 'p_value']]
    df['minus_log10_p_value'] = -np.log10(df['p_value'])
    df_min=min(df['minus_log10_p_value'])
    df_max=max(df['minus_log10_p_value'])
    axes_str=f'{celltypes_of_interest[count]}_{scores_of_interest[count]}_{condition_pairs_of_interest[count][0]}vs{condition_pairs_of_interest[count][1]}'
    args = dict(ax=axes[axes_str], x='minus_log10_p_value', y='score', data=df)
    p_cutoff=-np.log10(0.05/(3*11*3*5))
    
    if count==3:
        ax=axes[axes_str]
        ax.grid(False)
        ax.text(0.17, 0.35, 'No data shown because:\n$P_{subsampled} = 0 = P_{original}$ for all subsamples', fontsize=23, color='red')
    else:
        ax = sns.boxplot(**args)
        trans = ax.get_xaxis_transform()
        ax.axvline(x=p_cutoff, color='blue')
        ax.axvline(x=p, color='red')
    if count==4:
        ax=axes[axes_str]
        ax.grid(False)
        ax.text(280.17, 0.40, '$P_{subsampled} = 0$', fontsize=17, color='black', rotation=90)
    if count==5:
        ax=axes[axes_str]
        ax.grid(False)
        ax.text(291.17, 0.40, '$P_{subsampled} = 0$', fontsize=17, color='black', rotation=90)

    
    # # This will add title to subplot:
    ax.set_title(title, fontsize=30)
    ax.set_xlabel(xlabel, fontsize=25)
    ax.set_ylabel(None)
    # ax.set_yticks(None)
    ax.tick_params(axis='x', which='major', labelsize=20)
    ax.tick_params(axis='y', which='major', labelsize=20, labelleft=False)
    ax.grid(False)
    # if count==3:
    #     trans = ax.get_xaxis_transform()
    #     ax.text(4.13, 0.89, '$P_{original}=0$', transform=trans, fontsize=23, color='red')
    #     # ax.legend(title='$P_{original}=0$', bbox_to_anchor=(-0.00, -0.00), loc='lower left', fontsize=23, labelcolor='red')

# ========== Generate legend: ===========
legend_elements = [Line2D([0], [0], color='red', lw=4, label='$-log_{10}$'+'($P_{original}$)'),
                   Line2D([0], [0], color='blue', lw=4, label='$-log_{10}$'+'($P_{cutoff}$)')]
plt.legend(handles=legend_elements, bbox_to_anchor=(-0.30, -1.20), loc='lower left', fontsize=23)

    
# ========== Generate, save and show final plot (fig1): ==========
plt.subplots_adjust(wspace=0.4, hspace=0.9)
plt.savefig('fig5_subsampled_patients.pdf', format='pdf', bbox_inches='tight')
plt.show()











































