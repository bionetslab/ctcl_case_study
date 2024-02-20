import pandas as pd
import scipy.stats as stats
import os
import itertools as itt
from IPython.display import set_matplotlib_formats
import seaborn as sns
sns.set_theme(style="whitegrid")
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import itertools
from statannotations.Annotator import Annotator
from decimal import Decimal
import numpy as np

if __name__ == '__main__':
    p_values_cell_type=pd.read_csv(os.path.join('results', 'p_values_cell_type.csv'))
    p_values_cell_type_shuffled_labels=pd.read_csv(os.path.join('results', 'p_values_cell_type_shuffled_labels.csv'))
    celltypes_=list(np.unique(p_values_cell_type['cell_type']))
    heterogeneity_measures=['entropy', 'homophily', 'egophily']
    conditions_=list(set(np.unique(p_values_cell_type.condition_1)).union(set(np.unique(p_values_cell_type.condition_2))))
    conditions_abbreviations_dict=dict(zip(conditions_, ['AD', 'PSO', 'CTCL']))
    disease_combinations=list(itertools.combinations(conditions_, 2))
    radii=[1, 2, 3, 4, 5]
    no_of_rows=len(radii)
    no_of_cols=len(disease_combinations)
    for i in celltypes_:
        for heterogeneity_measure in heterogeneity_measures:
            matplotlib.rc('xtick', labelsize=40) 
            matplotlib.rc('ytick', labelsize=40) 
            fig, axes = plt.subplots(no_of_rows, no_of_cols, figsize=(40, 50), sharey=True, constrained_layout=True)
            plt.subplots_adjust(hspace=10, wspace=10)
            plt_count=-1
            row=-1
            for radius in radii:
                row+=1
                if heterogeneity_measure=='egophily':
                    local_heterogeneity_measure=f'{heterogeneity_measure}_{radius}'
                    suptitle=f'\n{heterogeneity_measure} ({i})'
                else:
                    local_heterogeneity_measure=f'local_{heterogeneity_measure}_{radius}'
                    suptitle=f'\nlocal {heterogeneity_measure} ({i})'
                col=-1
                for dis_comb in disease_combinations:
                    plt_count+=1
                    col+=1
                    if len(p_values_cell_type[p_values_cell_type['cell_type']==i][p_values_cell_type['condition_1']==dis_comb[0]][p_values_cell_type['condition_2']==dis_comb[1]][p_values_cell_type['score']==local_heterogeneity_measure]['p_value_adj'])==0:
                        df=p_values_cell_type_shuffled_labels[p_values_cell_type_shuffled_labels['cell_type']==i][p_values_cell_type_shuffled_labels['condition_1']==dis_comb[1]][p_values_cell_type_shuffled_labels['condition_2']==dis_comb[0]][p_values_cell_type_shuffled_labels['score']==local_heterogeneity_measure]
                        p=float(p_values_cell_type[p_values_cell_type['cell_type']==i][p_values_cell_type['condition_1']==dis_comb[1]][p_values_cell_type['condition_2']==dis_comb[0]][p_values_cell_type['score']==local_heterogeneity_measure]['p_value_adj'].values[0])
                    else:
                        df=p_values_cell_type_shuffled_labels[p_values_cell_type_shuffled_labels['cell_type']==i][p_values_cell_type_shuffled_labels['condition_1']==dis_comb[0]][p_values_cell_type_shuffled_labels['condition_2']==dis_comb[1]][p_values_cell_type_shuffled_labels['score']==local_heterogeneity_measure]
                        p=float(p_values_cell_type[p_values_cell_type['cell_type']==i][p_values_cell_type['condition_1']==dis_comb[0]][p_values_cell_type['condition_2']==dis_comb[1]][p_values_cell_type['score']==local_heterogeneity_measure]['p_value_adj'].values[0])
                    df['minus_log10_p_value_adj'] = -np.log10(df['p_value_adj'])
                    ax=sns.histplot(ax=axes[row, col], data=df, x="minus_log10_p_value_adj", kde=True)
                    ax.axvline(x=p, color='red')
                    # This will add title to subplot:
                    ax.set_title(str(dis_comb), fontsize=60)
                    # This will add label to x-axis:
                    ax.set_xlabel("-log10(p_value_adjusted)\n", fontsize=50) 
                    # This will add label to y-axis:
                    ax.set_ylabel(f"$radius={radius}$", fontsize=50)
            fig.suptitle(suptitle, weight='bold', fontsize=80, y=1.02)
            plt.xlabel('radius') # fontsize=25, labelpad=20
            plt.ylabel(f'{heterogeneity_measure} score') # fontsize=25, labelpad=20
            plt.tight_layout(pad=1.08, h_pad=1.08, w_pad=1.08)
            plt.savefig(f'{i}_{heterogeneity_measure}_shuffled_labels.jpg', format='jpg', bbox_inches='tight')
            plt.show()
    
    