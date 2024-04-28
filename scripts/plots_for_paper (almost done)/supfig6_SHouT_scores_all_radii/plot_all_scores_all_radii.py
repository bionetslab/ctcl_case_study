import pandas as pd
import os
import seaborn as sns
sns.set_theme(style="whitegrid")
import matplotlib.pyplot as plt
# import numpy as np
# from matplotlib.lines import Line2D
from statannotations.Annotator import Annotator
import itertools
from decimal import Decimal


# read in cell results and p-values per score per condition pair per cell type:
cell_results=pd.read_csv(os.path.join('../../../results', 'cell_results.csv'))
p_values_cell_type=pd.read_csv(os.path.join('../../../results', 'p_values_cell_type.csv'))

palette={"AD":(0.12156862745098039, 0.4666666666666667, 0.7058823529411765),
         "PSO": (0.17254901960784313, 0.6274509803921569, 0.17254901960784313),
         "CTCL": (1.0, 0.4980392156862745, 0.054901960784313725)}

# ========== Load and pre-process data required for generating plots: ==========
radii=[1, 2, 3, 4, 5]
conditions_=['Eczema', 'T-Cell Lymphoma', 'Psoriasis']
conditions=['AD', 'CTCL', 'PSO']
conditions_abbreviations_dict=dict(zip(conditions_, conditions))
cell_results=cell_results.replace({"condition": conditions_abbreviations_dict})
p_values_cell_type=p_values_cell_type.replace({"condition_1": conditions_abbreviations_dict})
p_values_cell_type=p_values_cell_type.replace({"condition_2": conditions_abbreviations_dict})
celltypes=sorted(list(set(cell_results.cell_type)))
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

celltypes_of_interest=["T-cells", "T-cells", "Basal keratinocytes", "Basal keratinocytes"]
# condition_pairs_of_interest=[('CTCL', 'AD'), ('CTCL', 'PSO'), ('AD', 'PSO'), ('CTCL', 'AD'), ('CTCL', 'PSO'), ('AD', 'PSO'), ('CTCL', 'AD'), ('CTCL', 'PSO'), ('AD', 'PSO'), ('CTCL', 'AD'), ('CTCL', 'PSO'), ('AD', 'PSO')]
scores_of_interest=["entropy", "egophily", "entropy", "homophily"]

layout = [
    ["T-cells_entropy"], ["T-cells_egophily"],
    # --------
    ["Basal keratinocytes_entropy"], ["Basal keratinocytes_homophily"]
]
disease_combinations=list(itertools.combinations(conditions, 2))
fig, axes = plt.subplot_mosaic(layout, figsize=(25,25))
title_loc='left'
xlabel=None
legend_loc=None
for count in range(4):
    if count==0:
        legend=True
    else:
        legend=None
    axes_str=f'{celltypes_of_interest[count]}_{scores_of_interest[count]}'
    celltype=celltypes_of_interest[count]
    if scores_of_interest[count]=='egophily':
        local_heterogeneity_measure=[f'{scores_of_interest[count]}_{radius}' for radius in radii]
        local_heterogeneity_score=f'{scores_of_interest[count].capitalize()}'
        var_name=f'{scores_of_interest[count]}_measure'
        value_name=f'{scores_of_interest[count]}_score'
        ylabel='Egophily'
    else:
        local_heterogeneity_measure=[f'local_{scores_of_interest[count]}_{radius}' for radius in radii]
        local_heterogeneity_score=f'Local {scores_of_interest[count]}'
        var_name=f'local_{scores_of_interest[count]}_measure'
        value_name=f'local_{scores_of_interest[count]}_score'
        ylabel=f'Local {scores_of_interest[count]}'
    
    if celltypes_of_interest[count]=='T-cells':
            plot_title=f'{local_heterogeneity_score}, {celltypes_of_interest[count]}'
    elif celltypes_of_interest[count]=='Basal keratinocytes':
            plot_title=f'{local_heterogeneity_score}, {celltypes_of_interest[count].lower()}'
    
    columns=local_heterogeneity_measure+['condition']
    data_=cell_results[cell_results['cell_type']==celltype][columns]
    data=data_.melt('condition', var_name=var_name, value_name=value_name)
    if palette:
        args = dict(x=var_name, y=value_name, data=data, hue="condition", hue_order=list(conditions), order=local_heterogeneity_measure, palette=palette)
    else:
        args = dict(x=var_name, y=value_name, data=data, hue="condition", hue_order=list(conditions), order=local_heterogeneity_measure)
    pairs=[]
    pvals=[]
    for j in local_heterogeneity_measure:
        for k in disease_combinations:
            pairs.append(((j, k[0]), (j, k[1])))
            if len(p_values_cell_type[p_values_cell_type['cell_type']==celltype][p_values_cell_type['condition_1']==k[0]][p_values_cell_type['condition_2']==k[1]][p_values_cell_type['score']==j]['p_value_adj'])==0:
                p_=str(p_values_cell_type[p_values_cell_type['cell_type']==celltype][p_values_cell_type['condition_1']==k[1]][p_values_cell_type['condition_2']==k[0]][p_values_cell_type['score']==j]['p_value_adj'].values[0])
            else:
                p_=str(p_values_cell_type[p_values_cell_type['cell_type']==celltype][p_values_cell_type['condition_1']==k[0]][p_values_cell_type['condition_2']==k[1]][p_values_cell_type['score']==j]['p_value_adj'].values[0])
            if float(p_)<0.05:
                pvals.append('%.2E' % Decimal(p_))
            else:
                pvals.append('$ns$')
    # pairs_pvals_dict=dict(zip(pairs, pvals))
    # fig, axes_ = plt.subplots(figsize=(20,10))
    # sns.set(font_scale = 1.2)
    # sns.set_style("white")
    # ax = sns.violinplot(**args, cut=0)
    # annot = Annotator(ax, pairs, **args)
    # annot.set_custom_annotations(pvals)
    # annot.annotate()
    # plt.close()
    # pairs_corrected=[]
    # pvals_corrected=[]
    # for j in range(len(annot.__dict__['annotations'])):
    #     pair_1=annot.__dict__['annotations'][j].__dict__['structs'][0]['group']
    #     pair_2=annot.__dict__['annotations'][j].__dict__['structs'][1]['group']   
    #     pairs_corrected.append((pair_1, pair_2))
    #     pvals_corrected.append(pairs_pvals_dict[(pair_1, pair_2)])
    sns.set(font_scale = 1.2)
    sns.set_style("white")
    ax = sns.violinplot(ax=axes[axes_str], **args, cut=0)
    title_prefix=None
    if title_prefix:
        if plot_title:
            plot_title=r"$\bf{" + title_prefix + "}$" + plot_title
        else:
            plot_title=r"$\bf{" + title_prefix + "}$"
    ax.set_title(plot_title, loc=title_loc, fontsize=30) # fontsize=25, pad=20
    ax.set_xlabel(xlabel, fontsize=27) # fontsize=25, labelpad=20
    ax.set_ylabel(ylabel, fontsize=27, labelpad=5) # fontsize=25, labelpad=20
    # # ax.set(xticklabels=[]) # size=20
    # annot = Annotator(ax, pairs, **args)
    # annot.configure(text_format='simple', loc='inside', verbose=2, fontsize=20) # fontsize=25
    # annot.set_custom_annotations(pvals_corrected)
    # annot.annotate()
    # ax.legend(title='Condition', bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
    # ax.legend(title=None, loc=legend_loc, borderaxespad=0, fontsize=20)
    ax.get_legend().set_visible(False)
    # ax.legend.remove()
    yticks=[]
    for j in ax.get_yticks():
        yticks.append(round(j,1))
    ax.set_yticklabels(yticks, size = 23) # size = 20
    xticks=[]
    for r in radii:
        xticks.append(f'$r={r}$')
    ax.set_xticklabels(xticks, size = 23) # size = 20
    # ax.tick_params(axis='x', which='major', labelsize=20)
    ax.grid(False)
    if legend:
        sns.move_legend(ax, loc=legend_loc,
                        title=None, fontsize=20,
                        bbox_to_anchor=(1.0, 1.05))
    else:
        ax.legend([],[], frameon=False)
    
    
# ========== Generate, save and show final plot (fig1): ==========
# plt.subplots_adjust(wspace=0.02, hspace=0.5)
plt.subplots_adjust(hspace=0.50)
plt.savefig('plots_all_scores_all_radii.pdf', format='pdf', bbox_inches='tight')
# fig.tight_layout(pad=100.0)
plt.show()















