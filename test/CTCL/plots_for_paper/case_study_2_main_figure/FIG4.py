import pandas as pd
import os
import seaborn as sns
sns.set_theme(style="whitegrid")
import matplotlib.pyplot as plt
import numpy as np
import itertools
from decimal import Decimal

# ========== Load and pre-process data required for generating plots: ==========
celltypes_of_interest=["T-cells", "Basal keratinocytes"]
layout = [
    ["T-cells", "T-cells", "Basal keratinocytes", "Basal keratinocytes"],
    ["T-cells", "T-cells", "Basal keratinocytes", "Basal keratinocytes"],
]

# celltypes_of_interest=["B-cells", "Basal keratinocytes",
#                        "Endothelial cells", "Fibroblasts",
#                        "Langerhans cells", "Macrophages",
#                        "Melanocytes", "Smooth muscle cells",
#                        "Suprabasal keratinocytes", "T-cells",
#                        "Unknown"]
# layout = [
#     ["B-cells", "B-cells", "Basal keratinocytes", "Basal keratinocytes", "Endothelial cells", "Endothelial cells", "Fibroblasts", "Fibroblasts", "Langerhans cells", "Langerhans cells", "Macrophages", "Macrophages"],
#     ["B-cells", "B-cells", "Basal keratinocytes", "Basal keratinocytes", "Endothelial cells", "Endothelial cells", "Fibroblasts", "Fibroblasts", "Langerhans cells", "Langerhans cells", "Macrophages", "Macrophages"],
#     [".", "Melanocytes", "Melanocytes", "Smooth muscle cells", "Smooth muscle cells", "Suprabasal keratinocytes", "Suprabasal keratinocytes", "T-cells", "T-cells", "Unknown", "Unknown", "."],
#     [".", "Melanocytes", "Melanocytes", "Smooth muscle cells", "Smooth muscle cells", "Suprabasal keratinocytes", "Suprabasal keratinocytes", "T-cells", "T-cells", "Unknown", "Unknown", "."]
# ]

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

p_values_cell_type=pd.read_csv(os.path.join('../../results', 'p_values_cell_type.csv'))
p_values_cell_type=p_values_cell_type.replace({"condition_1": conditions_abbreviations_dict})
p_values_cell_type=p_values_cell_type.replace({"condition_2": conditions_abbreviations_dict})
p_values_cell_type_shuffled_labels=pd.read_csv(os.path.join('../../results', 'p_values_cell_type_shuffled_labels.csv'))
p_values_cell_type_shuffled_labels=p_values_cell_type_shuffled_labels.replace({"condition_1": conditions_abbreviations_dict})
p_values_cell_type_shuffled_labels=p_values_cell_type_shuffled_labels.replace({"condition_2": conditions_abbreviations_dict})
celltypes_=list(np.unique(p_values_cell_type['cell_type']))
heterogeneity_measures=['entropy', 'homophily', 'egophily']
conditions_=list(set(np.unique(p_values_cell_type.condition_1)).union(set(np.unique(p_values_cell_type.condition_2))))
conditions_abbreviations_dict=dict(zip(conditions_, ['AD', 'PSO', 'CTCL']))
disease_combinations=list(itertools.combinations(conditions_, 2))
radii=[1, 2, 3, 4, 5]
no_of_rows=len(radii)
no_of_cols=len(disease_combinations)

for heterogeneity_measure in heterogeneity_measures:
    fig = plt.figure(figsize=(25, 25))
    if heterogeneity_measure=='egophily':
        fig.suptitle('Egophily ($r=5$)', fontsize=15)
        local_heterogeneity_measure=f'{heterogeneity_measure}_5'
    else:
        fig.suptitle(f'Local {heterogeneity_measure} ($r=5$)', fontsize=15)
        local_heterogeneity_measure=f'local_{heterogeneity_measure}_5'
    subfigs = fig.subfigures(3, 1, hspace=0.00)
    count=-1
    for dis_comb in disease_combinations:
        count+=1
        if count==0:
            title_str=r"$\bf{" + 'Local\:' + heterogeneity_measure + "}$"+ r' $\bf{(}$' + r"$\bf{" + 'r=5' + ')' + "}$" + f'\n({dis_comb[0]}, {dis_comb[1]})'
            subfigs[count].suptitle(title_str, fontsize=30, y=1.05)
        else:
            subfigs[count].suptitle(f'({dis_comb[0]}, {dis_comb[1]})', fontsize=30, y=0.99)
        textstr_ylabel = 'Count'
        props = dict(boxstyle='round', facecolor='white', alpha=0.00)
        subfigs[count].text(0.08, 0.40, textstr_ylabel, fontsize=25,
                verticalalignment='top', bbox=props,
                rotation = 90, 
                rotation_mode = 'anchor')
        textstr_xlabel = "$-log_{10}$"+"($p_{shuffled}$)"
        props = dict(boxstyle='round', facecolor='white', alpha=0.00)
        subfigs[count].text(0.45, 0.08, textstr_xlabel, fontsize=25,
                verticalalignment='top', bbox=props,
                rotation = 0, 
                rotation_mode = 'anchor')
        
        
        # axes = subfigs[count].subplot_mosaic(layout, sharex=True, sharey=True)
        axes = subfigs[count].subplot_mosaic(layout)
        if count==0:
            subfigs[count].set_facecolor('0.85')
        elif count==1:
            subfigs[count].set_facecolor('0.75')
        elif count==2:
            subfigs[count].set_facecolor('0.85')
        
        for i in celltypes_:
            if i not in celltypes_of_interest:
                continue
            if len(p_values_cell_type[p_values_cell_type['cell_type']==i][p_values_cell_type['condition_1']==dis_comb[0]][p_values_cell_type['condition_2']==dis_comb[1]][p_values_cell_type['score']==local_heterogeneity_measure]['p_value'])==0:
                p=float(p_values_cell_type[p_values_cell_type['cell_type']==i][p_values_cell_type['condition_1']==dis_comb[1]][p_values_cell_type['condition_2']==dis_comb[0]][p_values_cell_type['score']==local_heterogeneity_measure]['p_value'].values[0])
                p=-np.log10(p)
            else:
                p=float(p_values_cell_type[p_values_cell_type['cell_type']==i][p_values_cell_type['condition_1']==dis_comb[0]][p_values_cell_type['condition_2']==dis_comb[1]][p_values_cell_type['score']==local_heterogeneity_measure]['p_value'].values[0])
                p=-np.log10(p)
            df=p_values_cell_type_shuffled_labels[p_values_cell_type_shuffled_labels['cell_type']==i][p_values_cell_type_shuffled_labels['condition_1']==dis_comb[1]][p_values_cell_type_shuffled_labels['condition_2']==dis_comb[0]][p_values_cell_type_shuffled_labels['score']==local_heterogeneity_measure]
            if len(df)==0: # The (condition_1, condition_2) column pair values are just filled (because the snytzetic data contains pairs in alphabetical order, i.e., condition_1 value always smaller than condition_2 value alphabetically.)
                df=p_values_cell_type_shuffled_labels[p_values_cell_type_shuffled_labels['cell_type']==i][p_values_cell_type_shuffled_labels['condition_1']==dis_comb[0]][p_values_cell_type_shuffled_labels['condition_2']==dis_comb[1]][p_values_cell_type_shuffled_labels['score']==local_heterogeneity_measure]
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
            ax=sns.histplot(ax=axes[i], data=df, x="minus_log10_p_value", kde=True)
            # ax.axvline(x=p, color='red')
            text_str=str('%.2E' % Decimal(p))                  
            # Show texts
            # ax.text(0.1, 0.5, 'p-value of actual data={text_str}', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
            ax.text(0.20, 0.89, '$p_{actual}$='+text_str, transform=ax.transAxes, fontsize=20)
            
            # This will add title to subplot:
            ax.set_title(celltypes_abbreviations_dict[i], fontsize=30)
            # This will add label to x-axis:
            # ax.set_xlabel("-log10($p_{shuffled}$)\n", fontsize=25) 
            ax.set_xlabel(None, fontsize=25) 
            # This will add label to y-axis:
            # ax.set_ylabel("Count", fontsize=25)
            ax.set_ylabel(None, fontsize=25)
            ax.tick_params(axis='x', which='major', labelsize=20)
            ax.tick_params(axis='y', which='major', labelsize=20)
    # ========== Generate, save and show final plot (fig1): ==========
    plt.subplots_adjust(wspace=0.49, hspace=1.8)
    plt.savefig(f'{heterogeneity_measure}_shuffled_labels.pdf', format='pdf', bbox_inches='tight')
    plt.show()





