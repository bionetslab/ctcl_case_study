import pandas as pd
import os
import seaborn as sns
sns.set_theme(style="whitegrid")
import matplotlib.pyplot as plt
import numpy as np
import itertools
from decimal import Decimal

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

p_values_cell_type=pd.read_csv(os.path.join('../../results', 'p_values_cell_type.csv'))
p_values_cell_type=p_values_cell_type.replace({"condition_1": conditions_abbreviations_dict})
p_values_cell_type=p_values_cell_type.replace({"condition_2": conditions_abbreviations_dict})
p_values_cell_type_subsampled_patientwise=pd.read_csv(os.path.join('../../results', 'p_values_cell_type_subsampled_patientwise.csv'))
p_values_cell_type_subsampled_patientwise=p_values_cell_type_subsampled_patientwise.replace({"condition_1": conditions_abbreviations_dict})
p_values_cell_type_subsampled_patientwise=p_values_cell_type_subsampled_patientwise.replace({"condition_2": conditions_abbreviations_dict})
celltypes_=list(np.unique(p_values_cell_type['cell_type']))
heterogeneity_measures=['entropy', 'homophily', 'egophily']
conditions_=list(set(np.unique(p_values_cell_type.condition_1)).union(set(np.unique(p_values_cell_type.condition_2))))
conditions_abbreviations_dict=dict(zip(conditions_, ['AD', 'PSO', 'CTCL']))
disease_combinations=list(itertools.combinations(conditions_, 2))
radii=[1, 5]
no_of_rows=len(radii)
no_of_cols=len(disease_combinations)

layout = [
    ["B-cells", "B-cells", "Basal keratinocytes", "Basal keratinocytes", "Endothelial cells", "Endothelial cells", "Fibroblasts", "Fibroblasts", "Langerhans cells", "Langerhans cells", "Macrophages", "Macrophages"],
    ["B-cells", "B-cells", "Basal keratinocytes", "Basal keratinocytes", "Endothelial cells", "Endothelial cells", "Fibroblasts", "Fibroblasts", "Langerhans cells", "Langerhans cells", "Macrophages", "Macrophages"],
    [".", "Melanocytes", "Melanocytes", "Smooth muscle cells", "Smooth muscle cells", "Suprabasal keratinocytes", "Suprabasal keratinocytes", "T-cells", "T-cells", "Unknown", "Unknown", "."],
    [".", "Melanocytes", "Melanocytes", "Smooth muscle cells", "Smooth muscle cells", "Suprabasal keratinocytes", "Suprabasal keratinocytes", "T-cells", "T-cells", "Unknown", "Unknown", "."]
]

for heterogeneity_measure in heterogeneity_measures:
    fig = plt.figure(figsize=(25, 25))
    if heterogeneity_measure=='egophily':
        fig.suptitle('Egophily ($r=5$)', fontsize=15)
        local_heterogeneity_measure=[f'{heterogeneity_measure}_{radius}' for radius in radii]
    else:
        fig.suptitle(f'Local {heterogeneity_measure} ($r=5$)', fontsize=15)
        local_heterogeneity_measure=[f'local_{heterogeneity_measure}_{radius}' for radius in radii]
    subfigs = fig.subfigures(3, 1, hspace=0.00)
    count=-1
    data_heterogeneity_measure=[]
    for dis_comb in disease_combinations:
        count+=1
        if count==0:
            title_str=r"$\bf{" + 'Local\:' + heterogeneity_measure + "}$"+ r' $\bf{(}$' + r"$\bf{" + 'r=5' + ')' + "}$" + f'\n{dis_comb}'
            subfigs[count].suptitle(title_str, fontsize=30, y=1.05)
        else:
            subfigs[count].suptitle(str(dis_comb), fontsize=30)
        axes = subfigs[count].subplot_mosaic(layout)
        
        if count==0:
            subfigs[count].set_facecolor('0.85')
        elif count==1:
            subfigs[count].set_facecolor('0.75')
        elif count==2:
            subfigs[count].set_facecolor('0.85')
        
        for i in celltypes_:
            data_disease_combinations=[]
            for k in local_heterogeneity_measure:
                if len(p_values_cell_type_subsampled_patientwise[p_values_cell_type_subsampled_patientwise['cell_type']==i][p_values_cell_type_subsampled_patientwise['condition_1']==dis_comb[0]][p_values_cell_type_subsampled_patientwise['condition_2']==dis_comb[1]][p_values_cell_type_subsampled_patientwise['score']==k]['p_value'])==0:
                    data_=p_values_cell_type_subsampled_patientwise[p_values_cell_type_subsampled_patientwise['cell_type']==i][p_values_cell_type_subsampled_patientwise['condition_1']==dis_comb[1]][p_values_cell_type_subsampled_patientwise['condition_2']==dis_comb[0]][p_values_cell_type_subsampled_patientwise['score']==k]
                else:
                    data_=p_values_cell_type_subsampled_patientwise[p_values_cell_type_subsampled_patientwise['cell_type']==i][p_values_cell_type_subsampled_patientwise['condition_1']==dis_comb[0]][p_values_cell_type_subsampled_patientwise['condition_2']==dis_comb[1]][p_values_cell_type_subsampled_patientwise['score']==k]
                data_disease_combinations.append(data_)
            data_disease_combinations=pd.concat(data_disease_combinations, axis=0)
            data_heterogeneity_measure.append(data_disease_combinations)
            df=data_disease_combinations.copy()
            df=df[['score', 'p_value']]
            args = dict(ax=axes[i], x='score', y='p_value', data=df)
            ax = sns.boxplot(**args)
            medians_ = df.groupby(['score'])['p_value'].median()
            medians=[]
            medians_str=[]
            for median_ in medians_:
                medians.append(float('%.2E' % Decimal(median_)))
                medians_str.append(str('%.2E' % Decimal(median_)))
            vertical_offset = df['p_value'].median() * 0.05 # offset from median for display
            xticklabels=[f'$r={radius}$' for radius in radii]
            ax.set_xticklabels(xticklabels, fontsize=20)
            
            for xtick in ax.get_xticks():
                ax.text(xtick,medians[xtick] + 1.5*vertical_offset,medians_str[xtick], 
                        horizontalalignment='center',size=19,color='black',weight='semibold')
            # This will add title to subplot:
            ax.set_title(celltypes_abbreviations_dict[i], fontsize=30)
            # This will add label to x-axis:
            ax.set_xlabel(None, fontsize=25) 
            # This will add label to y-axis:
            ax.set_ylabel(None, fontsize=25)
            ax.tick_params(axis='x', which='major', labelsize=20)
            ax.tick_params(axis='y', which='major', labelsize=20)
    
    # ========== Generate, save and show final plot (fig1): ==========
    # plt.subplots_adjust(wspace=0.02, hspace=0.5)
    plt.subplots_adjust(wspace=1.50, hspace=1.8)
    plt.savefig(f'{heterogeneity_measure}_subsampled_patientwise.pdf', format='pdf', bbox_inches='tight')
    # fig.tight_layout(pad=100.0)
    plt.show()



