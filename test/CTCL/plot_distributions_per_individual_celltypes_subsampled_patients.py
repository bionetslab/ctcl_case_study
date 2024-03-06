import pandas as pd
import scipy.stats as stats
import os
import itertools as itt
from IPython.display import set_matplotlib_formats
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import itertools
from statannotations.Annotator import Annotator
from decimal import Decimal

if __name__ == '__main__':
    p_values_cell_type_subsampled_patientwise=pd.read_csv(os.path.join('results', 'p_values_cell_type_subsampled_patientwise.csv'))    
    celltypes_=list(np.unique(p_values_cell_type_subsampled_patientwise['cell_type']))
    heterogeneity_measures=['entropy', 'homophily', 'egophily']
    conditions_=list(set(np.unique(p_values_cell_type_subsampled_patientwise.condition_1)).union(set(np.unique(p_values_cell_type_subsampled_patientwise.condition_2))))
    conditions_abbreviations_dict=dict(zip(conditions_, ['AD', 'PSO', 'CTCL']))
    disease_combinations=list(itertools.combinations(conditions_, 2))
    radii=[1, 2, 3, 4, 5]
    for i in celltypes_:
        for heterogeneity_measure in heterogeneity_measures:
            if heterogeneity_measure=='egophily':
                if not(i in ['T-cells', 'Fibroblasts', 'Suprabasal keratinocytes']):
                    continue
                local_heterogeneity_measure=[f'{heterogeneity_measure}_{radius}' for radius in radii]
                var_name=f'{heterogeneity_measure}_measure'
                value_name=f'{heterogeneity_measure}_score'
                suptitle=f'\n{heterogeneity_measure} ({i})'
            else:
                if heterogeneity_measure=='entropy':
                    if not(i in ['Macrophages', 'Basal keratinocytes', 'Endothelial cells', 'T-cells']):
                        continue
                elif heterogeneity_measure=='homophily':
                    if not(i in ['T-cells', 'Basal keratinocytes']):
                        continue
                local_heterogeneity_measure=[f'local_{heterogeneity_measure}_{radius}' for radius in radii]
                var_name=f'local_{heterogeneity_measure}_measure'
                value_name=f'local_{heterogeneity_measure}_score'
                suptitle=f'\nlocal {heterogeneity_measure} ({i})'
            sns.set(font_scale = 1.2)
            sns.set_style("white")
            fig, axes = plt.subplots(len(disease_combinations), 1, figsize=(10, 10), constrained_layout=True)
            plt.subplots_adjust(hspace=10, wspace=10)
            row=-1
            data_heterogeneity_measure=[]
            for j in disease_combinations:
                row+=1
                data_disease_combinations=[]
                for k in local_heterogeneity_measure:
                    if len(p_values_cell_type_subsampled_patientwise[p_values_cell_type_subsampled_patientwise['cell_type']==i][p_values_cell_type_subsampled_patientwise['condition_1']==j[0]][p_values_cell_type_subsampled_patientwise['condition_2']==j[1]][p_values_cell_type_subsampled_patientwise['score']==k]['p_value_adj'])==0:
                        data_=p_values_cell_type_subsampled_patientwise[p_values_cell_type_subsampled_patientwise['cell_type']==i][p_values_cell_type_subsampled_patientwise['condition_1']==j[1]][p_values_cell_type_subsampled_patientwise['condition_2']==j[0]][p_values_cell_type_subsampled_patientwise['score']==k]
                    else:
                        data_=p_values_cell_type_subsampled_patientwise[p_values_cell_type_subsampled_patientwise['cell_type']==i][p_values_cell_type_subsampled_patientwise['condition_1']==j[0]][p_values_cell_type_subsampled_patientwise['condition_2']==j[1]][p_values_cell_type_subsampled_patientwise['score']==k]
                    data_disease_combinations.append(data_)
                data_disease_combinations=pd.concat(data_disease_combinations, axis=0)
                data_heterogeneity_measure.append(data_disease_combinations)
                df=data_disease_combinations.copy()
                df=df[['score', 'p_value_adj']]
                # df['minus_log10_p_value_adj'] = -np.log10(df['p_value_adj'])
                # args = dict(ax=axes[row], x='score', y='minus_log10_p_value_adj', data=df)
                args = dict(ax=axes[row], x='score', y='p_value_adj', data=df)
                ax = sns.boxplot(**args)
                title=str(j)
                ax.set_title(f'{title}')
                ax.set_xlabel("radius\n")
                # ax.set_ylabel("-log10(p_value_adjusted)\n")
                ax.set_ylabel("p_value_adjusted\n")
                ax.set_xticklabels(radii) # , size=20
            fig.suptitle(suptitle, weight='bold', fontsize=15, y=1.001)
            plt.xlabel('radius') # fontsize=25, labelpad=20
            plt.ylabel(f'{heterogeneity_measure} score') # fontsize=25, labelpad=20
            plt.tight_layout(pad=1.08, h_pad=1.08, w_pad=1.08)
            plt.savefig(f'{i}_{heterogeneity_measure}_subsampled_patients.pdf', format='pdf', bbox_inches='tight')
            plt.show()
            data_heterogeneity_measure=pd.concat(data_heterogeneity_measure, axis=0)
            
