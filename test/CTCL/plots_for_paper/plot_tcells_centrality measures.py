import pandas as pd
import scipy.stats as stats
import os
import itertools as itt
from IPython.display import set_matplotlib_formats
import seaborn as sns
sns.set_theme(style="whitegrid")
import matplotlib.pyplot as plt
import numpy as np
import itertools
from statannotations.Annotator import Annotator
from decimal import Decimal

if __name__ == '__main__':
    var_name='centrality_measure'
    value_name='centrality_score'
    p_values_cell_type=pd.read_csv(os.path.join('../results', 'p_values_cell_type_squidpy_centralityScores.csv'))
    squidpy_centralityScores_results=pd.read_csv(os.path.join('../results', 'squidpy_centralityScores_results.csv'))
    celltypes_=list(set(p_values_cell_type['celltypes']))
    scores = ['degree_centrality', 'average_clustering', 'closeness_centrality']
    conditions_=list(np.unique(squidpy_centralityScores_results.condition))
    conditions_abbreviations_dict=dict(zip(conditions_, ['AD', 'PSO', 'CTCL']))
    disease_combinations=list(itertools.combinations(conditions_, 2))
    for i in celltypes_:
        data=squidpy_centralityScores_results[squidpy_centralityScores_results['celltypes']==i]
        columns=scores+['condition']
        data=data[columns]
        data=data.melt('condition', var_name=var_name, value_name=value_name)
        args = dict(x=var_name, y=value_name, data=data, hue="condition", hue_order=list(conditions_), order=scores)
        pairs=[]
        pvals=[]
        for k in scores:
            for j in disease_combinations:
                pairs.append(((k, j[0]), (k, j[1])))
                if len(p_values_cell_type[p_values_cell_type['condition_1']==j[0]][p_values_cell_type['condition_2']==j[1]][p_values_cell_type['score']==k])==0:
                    p_=str(p_values_cell_type[p_values_cell_type['condition_1']==j[1]][p_values_cell_type['condition_2']==j[0]][p_values_cell_type['score']==k]['p_value'][p_values_cell_type['celltypes']==i].values[0])
                else:
                    p_=str(p_values_cell_type[p_values_cell_type['condition_1']==j[0]][p_values_cell_type['condition_2']==j[1]][p_values_cell_type['score']==k]['p_value'][p_values_cell_type['celltypes']==i].values[0])
                if float(p_)<0.05:
                    pvals.append('%.2E' % Decimal(p_))
                else:
                    pvals.append('$ns$')
        plot_title=f'Centrality scores ({i})'
        pairs_pvals_dict=dict(zip(pairs, pvals))
        fig, axes = plt.subplots(figsize=(20,10))
        sns.set(font_scale = 1.2)
        sns.set_style("white")
        ax = sns.violinplot(**args, cut=0)
        annot = Annotator(ax, pairs, **args)
        annot.set_custom_annotations(pvals)
        annot.annotate()
        pairs_corrected=[]
        pvals_corrected=[]
        for j in range(len(annot.__dict__['annotations'])):
            pair_1=annot.__dict__['annotations'][j].__dict__['structs'][0]['group']
            pair_2=annot.__dict__['annotations'][j].__dict__['structs'][1]['group']   
            pairs_corrected.append((pair_1, pair_2))
            pvals_corrected.append(pairs_pvals_dict[(pair_1, pair_2)])
        fig, axes = plt.subplots(figsize=(20,10))
        sns.set(font_scale = 1.2)
        sns.set_style("white")
        ax = sns.violinplot(**args, cut=0)
        annot = Annotator(ax, pairs, **args)
        annot.configure(text_format='simple', loc='inside', verbose=2, fontsize=25)
        annot.set_custom_annotations(pvals_corrected)
        annot.annotate()
        plt.xlabel('', fontsize=25, labelpad=20)
        plt.ylabel(f'Centrality scores', fontsize=25, labelpad=20)
        plt.title(plot_title, fontsize=25, pad=20)
        xticks=[]
        for j in ax.get_yticks():
            xticks.append(round(j,1))
        ax.set_yticklabels(xticks, size = 20)
        plt.savefig(f'{plot_title}.pdf', format='pdf', bbox_inches='tight')
        fig.tight_layout()
        plt.show()
        
        
        