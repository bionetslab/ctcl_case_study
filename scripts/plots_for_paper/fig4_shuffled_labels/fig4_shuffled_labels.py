import pandas as pd
import os
import seaborn as sns
sns.set_theme(style="whitegrid")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D


# ========== Load and pre-process data required for generating plots: ==========
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
p_values_cell_type_shuffled_labels=pd.read_csv(os.path.join('../../../results', 'p_values_cell_type_shuffled_labels.csv'))
p_values_cell_type_shuffled_labels=p_values_cell_type_shuffled_labels.replace({"condition_1": conditions_abbreviations_dict})
p_values_cell_type_shuffled_labels=p_values_cell_type_shuffled_labels.replace({"condition_2": conditions_abbreviations_dict})

celltypes_of_interest=["T-cells", "T-cells", "T-cells", "T-cells", "T-cells", "T-cells", "Basal keratinocytes", "Basal keratinocytes", "Basal keratinocytes", "Basal keratinocytes", "Basal keratinocytes", "Basal keratinocytes",]
condition_pairs_of_interest=[('CTCL', 'AD'), ('CTCL', 'PSO'), ('AD', 'PSO'), ('CTCL', 'AD'), ('CTCL', 'PSO'), ('AD', 'PSO'), ('CTCL', 'AD'), ('CTCL', 'PSO'), ('AD', 'PSO'), ('CTCL', 'AD'), ('CTCL', 'PSO'), ('AD', 'PSO')]
scores_of_interest=["entropy", "entropy", "entropy", "egophily", "egophily", "egophily", "entropy", "entropy", "entropy", "homophily", "homophily", "homophily"]

# layout = [
#     ["T-cells_entropy_CTCLvsAD", "T-cells_entropy_CTCLvsAD", "T-cells_entropy_CTCLvsPSO", "T-cells_entropy_CTCLvsPSO"], # "T-cells_egophily_CTCLvsAD", "T-cells_egophily_CTCLvsAD"],
#     ["T-cells_entropy_ADvsPSO", "T-cells_entropy_ADvsPSO", "T-cells_egophily_CTCLvsAD", "T-cells_egophily_CTCLvsAD",],
#     ["T-cells_egophily_CTCLvsPSO", "T-cells_egophily_CTCLvsPSO", "T-cells_egophily_ADvsPSO", "T-cells_egophily_ADvsPSO",],
#     # --------
#     ["Basal keratinocytes_entropy_CTCLvsAD", "Basal keratinocytes_entropy_CTCLvsAD", "Basal keratinocytes_entropy_CTCLvsPSO", "Basal keratinocytes_entropy_CTCLvsPSO"],
#     ["Basal keratinocytes_entropy_ADvsPSO", "Basal keratinocytes_entropy_ADvsPSO", "Basal keratinocytes_homophily_CTCLvsAD", "Basal keratinocytes_homophily_CTCLvsAD"],
#     ["Basal keratinocytes_homophily_CTCLvsPSO", "Basal keratinocytes_homophily_CTCLvsPSO", "Basal keratinocytes_homophily_ADvsPSO", "Basal keratinocytes_homophily_ADvsPSO"]
# ]

# layout = [
#     ["T-cells_entropy_CTCLvsAD", "T-cells_egophily_CTCLvsAD", "Basal keratinocytes_homophily_CTCLvsAD"],
#     ["T-cells_entropy_CTCLvsPSO", "T-cells_egophily_CTCLvsAD", "Basal keratinocytes_homophily_CTCLvsAD"],
#     ["T-cells_entropy_ADvsPSO",  "T-cells_egophily_CTCLvsPSO", "Basal keratinocytes_homophily_CTCLvsPSO"],
#     ["Basal keratinocytes_entropy_CTCLvsAD", "T-cells_egophily_CTCLvsPSO", "Basal keratinocytes_homophily_CTCLvsPSO"],
#     ["Basal keratinocytes_entropy_CTCLvsPSO", "T-cells_egophily_ADvsPSO", "Basal keratinocytes_homophily_ADvsPSO"],
#     ["Basal keratinocytes_entropy_ADvsPSO", "T-cells_egophily_ADvsPSO", "Basal keratinocytes_homophily_ADvsPSO"]    
# ]

layout = [
    ["T-cells_entropy_CTCLvsAD", "Basal keratinocytes_entropy_CTCLvsAD", "T-cells_egophily_CTCLvsAD", "Basal keratinocytes_homophily_CTCLvsAD"],
    ["T-cells_entropy_CTCLvsPSO", "Basal keratinocytes_entropy_CTCLvsPSO", "T-cells_egophily_CTCLvsPSO", "Basal keratinocytes_homophily_CTCLvsPSO"],
    ["T-cells_entropy_ADvsPSO", "Basal keratinocytes_entropy_ADvsPSO", "T-cells_egophily_ADvsPSO", "Basal keratinocytes_homophily_ADvsPSO"]    
]

fig, axes = plt.subplot_mosaic(layout, figsize=(30,15))
for count in range(12):
    if scores_of_interest[count]=='egophily':
        # title='Egophily\n'
        local_heterogeneity_measure=f'{scores_of_interest[count]}_5'
    else:
        # title=f'Local {scores_of_interest[count]}\n'
        local_heterogeneity_measure=f'local_{scores_of_interest[count]}_5'
    if count not in [0,3,6,9]:
        title=None
    else:
        if count==0:
            title='Local entropy\n(T cells)\n'
        elif count==3:
            title='Local entropy\n(basal\nkeratinocytes)\n'
        elif count==6:
            title='Egophily\n(T cells)\n'
        else:
            title='Local homophily\n(basal\nkeratinocytes)\n'
    ylabel=f'{condition_pairs_of_interest[count][0]} vs. {condition_pairs_of_interest[count][1]}\n'
    if count not in [0,1,2]:
        ylabel=None
    # if celltypes_of_interest[count]=='T-cells':
    #     ylabel=f'T-cells\n({condition_pairs_of_interest[count][0]} vs. {condition_pairs_of_interest[count][1]})'
    # elif celltypes_of_interest[count]=='Basal keratinocytes':
    #     ylabel=f'Basal\nkeratinocytes\n({condition_pairs_of_interest[count][0]} vs. {condition_pairs_of_interest[count][1]})'
    xlabel="$-log_{10}$"+"($P_{shuffled}$)"
    if count not in [2,5,8,11]:
        xlabel=None
    if len(p_values_cell_type[p_values_cell_type['cell_type']==celltypes_of_interest[count]][p_values_cell_type['condition_1']==condition_pairs_of_interest[count][0]][p_values_cell_type['condition_2']==condition_pairs_of_interest[count][1]][p_values_cell_type['score']==local_heterogeneity_measure]['p_value'])==0:
        p=float(p_values_cell_type[p_values_cell_type['cell_type']==celltypes_of_interest[count]][p_values_cell_type['condition_1']==condition_pairs_of_interest[count][1]][p_values_cell_type['condition_2']==condition_pairs_of_interest[count][0]][p_values_cell_type['score']==local_heterogeneity_measure]['p_value'].values[0])
        p=-np.log10(p)
    else:
        p=float(p_values_cell_type[p_values_cell_type['cell_type']==celltypes_of_interest[count]][p_values_cell_type['condition_1']==condition_pairs_of_interest[count][0]][p_values_cell_type['condition_2']==condition_pairs_of_interest[count][1]][p_values_cell_type['score']==local_heterogeneity_measure]['p_value'].values[0])
        p=-np.log10(p)
    df=p_values_cell_type_shuffled_labels[p_values_cell_type_shuffled_labels['cell_type']==celltypes_of_interest[count]][p_values_cell_type_shuffled_labels['condition_1']==condition_pairs_of_interest[count][1]][p_values_cell_type_shuffled_labels['condition_2']==condition_pairs_of_interest[count][0]][p_values_cell_type_shuffled_labels['score']==local_heterogeneity_measure]
    if len(df)==0: # The (condition_1, condition_2) column pair values are just filled (because the snytzetic data contains pairs in alphabetical order, i.e., condition_1 value always smaller than condition_2 value alphabetically.)
        df=p_values_cell_type_shuffled_labels[p_values_cell_type_shuffled_labels['cell_type']==celltypes_of_interest[count]][p_values_cell_type_shuffled_labels['condition_1']==condition_pairs_of_interest[count][0]][p_values_cell_type_shuffled_labels['condition_2']==condition_pairs_of_interest[count][1]][p_values_cell_type_shuffled_labels['score']==local_heterogeneity_measure]
    df['minus_log10_p_value_adj'] = -np.log10(df['p_value_adj'])
    df['minus_log10_p_value'] = -np.log10(df['p_value'])
    # Sort by adjusted p-values:
    # df=df.sort_values(by=['minus_log10_p_value_adj'])
    # Sort by non-adjusted p-values:
    df=df.sort_values(by=['minus_log10_p_value'])
    df_=df.copy()
    # =================== Family-wise error correction (FWER) section ===================
    percentile=0.95
    scaling_ratio=1-percentile
    no_of_elements_to_be_dropped=int(round(scaling_ratio*len(df),1))
    # df = df.iloc[no_of_elements_to_be_dropped:]
    df[df.minus_log10_p_value > np.percentile(df.minus_log10_p_value,5)]
    # ===================================================================================
    axes_str=f'{celltypes_of_interest[count]}_{scores_of_interest[count]}_{condition_pairs_of_interest[count][0]}vs{condition_pairs_of_interest[count][1]}'
    ax=sns.histplot(ax=axes[axes_str], data=df, x="minus_log10_p_value", kde=True)
    ax.axvline(x=p, color='red')
    ax.set_title(title, fontsize=30)
    ax.set_xlabel(xlabel, fontsize=30) 
    ax.set_ylabel(ylabel, fontsize=30)
    ax.tick_params(axis='x', which='major', labelsize=20)
    ax.tick_params(axis='y', which='major', labelsize=20)
    ax.grid(False)
    if count==3:
        trans = ax.get_xaxis_transform()
        ax.text(1.85, 0.89, '$P_{original}=0$', transform=trans, fontsize=23, color='red')
        # ax.legend(title='$P_{original}=0$', bbox_to_anchor=(-0.00, -0.00), loc='lower left', fontsize=23, labelcolor='red')
    if count==9:
        ax=axes[axes_str]
        ax.grid(False)
        ax.text(5.50, 18.50, 'Count', fontsize=25, color='black', rotation=90)
    elif count==10:
        ax=axes[axes_str]
        ax.grid(False)
        ax.text(30.00, 18.50, 'Count', fontsize=25, color='black', rotation=90)
    elif count==11:
        ax=axes[axes_str]
        ax.grid(False)
        ax.text(9.25, 18.50, 'Count', fontsize=25, color='black', rotation=90)
    
    
# ========== Generate legend: ===========
legend_elements = [Line2D([0], [0], color='red', lw=4, label='$-log_{10}$'+'($P_{original}$)')]
plt.legend(handles=legend_elements, bbox_to_anchor=(-1.70, -0.75), loc='lower left', fontsize=27)
   
# ========== Generate, save and show final plot (fig1): ==========
# plt.subplots_adjust(wspace=0.3, hspace=0.6)
plt.savefig('fig4_shuffled_labels.pdf', format='pdf', bbox_inches='tight')
plt.show()





