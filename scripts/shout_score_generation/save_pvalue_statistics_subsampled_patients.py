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
    p_values_cell_type_subsampled_patientwise=pd.read_csv(os.path.join('../../results', 'p_values_cell_type_subsampled_patientwise.csv'))    
    sample_results=pd.read_csv(os.path.join('../../results', 'sample_results.csv'))
    celltypes_=list(np.unique(p_values_cell_type_subsampled_patientwise['cell_type']))
    heterogeneity_measures=['entropy', 'homophily', 'egophily']
    conditions_=list(set(np.unique(p_values_cell_type_subsampled_patientwise.condition_1)).union(set(np.unique(p_values_cell_type_subsampled_patientwise.condition_2))))
    conditions_abbreviations_dict=dict(zip(conditions_, ['AD', 'PSO', 'CTCL']))
    disease_combinations=list(itertools.combinations(conditions_, 2))
    radii=[1, 2, 3, 4, 5]
    pvalue_statistics=[]
    for i in celltypes_:
        data_all_heterogeneity_measures=[]
        all_heterogeneity_measures_statistics=[]
        for heterogeneity_measure in heterogeneity_measures:
            if heterogeneity_measure=='egophily':
                local_heterogeneity_measure=[f'{heterogeneity_measure}_{radius}' for radius in radii]
            else:
                local_heterogeneity_measure=[f'local_{heterogeneity_measure}_{radius}' for radius in radii]
            data_heterogeneity_measure=[]
            data_heterogeneity_measure_stats=[]
            for j in disease_combinations:
                data_disease_combinations=[]
                data_disease_combinations_stats=[]
                for k in local_heterogeneity_measure:
                    if len(p_values_cell_type_subsampled_patientwise[p_values_cell_type_subsampled_patientwise['cell_type']==i][p_values_cell_type_subsampled_patientwise['condition_1']==j[0]][p_values_cell_type_subsampled_patientwise['condition_2']==j[1]][p_values_cell_type_subsampled_patientwise['score']==k]['p_value_adj'])==0:
                        data_=p_values_cell_type_subsampled_patientwise[p_values_cell_type_subsampled_patientwise['cell_type']==i][p_values_cell_type_subsampled_patientwise['condition_1']==j[1]][p_values_cell_type_subsampled_patientwise['condition_2']==j[0]][p_values_cell_type_subsampled_patientwise['score']==k]
                    else:
                        data_=p_values_cell_type_subsampled_patientwise[p_values_cell_type_subsampled_patientwise['cell_type']==i][p_values_cell_type_subsampled_patientwise['condition_1']==j[0]][p_values_cell_type_subsampled_patientwise['condition_2']==j[1]][p_values_cell_type_subsampled_patientwise['score']==k]
                    data_stats=pd.DataFrame()
                    data_stats['condition_1']=[list(data_['condition_1'])[0]]
                    data_stats['condition_2']=[list(data_['condition_2'])[0]]
                    data_stats['score']=[list(data_['score'])[0]] 
                    data_stats['cell_type']=[list(data_['cell_type'])[0]]
                    # Compute dataframe of mean p_value per radius per diseaseCombination:
                    data_stats['p_value_mean']=[(pd.DataFrame(data_.mean()).T)['p_value'].values[0]]
                    data_stats['p_value_adj_mean']=[(pd.DataFrame(data_.mean()).T)['p_value_adj'].values[0]]
                    data_stats['statistic_mean']=[(pd.DataFrame(data_.mean()).T)['statistic'].values[0]]
                    # Compute dataframe of median p_value per radius per diseaseCombination:
                    data_stats['p_value_median']=[(pd.DataFrame(data_.median()).T)['p_value'].values[0]]
                    data_stats['p_value_adj_median']=[(pd.DataFrame(data_.median()).T)['p_value_adj'].values[0]]
                    data_stats['statistic_median']=[(pd.DataFrame(data_.median()).T)['statistic'].values[0]]
                    # Compute dataframe of spread of p_value per radius per diseaseCombination:
                    data_stats['p_value_spread']=[(pd.DataFrame(data_.max()).T)['p_value'].values[0]-(pd.DataFrame(data_.min()).T)['p_value'].values[0]]
                    data_stats['p_value_adj_spread']=[(pd.DataFrame(data_.max()).T)['p_value'].values[0]-(pd.DataFrame(data_.min()).T)['p_value_adj'].values[0]]
                    data_stats['statistic_spread']=[(pd.DataFrame(data_.max()).T)['p_value'].values[0]-(pd.DataFrame(data_.min()).T)['statistic'].values[0]]
                    data_disease_combinations.append(data_)
                    data_disease_combinations_stats.append(data_stats)
                data_disease_combinations=pd.concat(data_disease_combinations, axis=0)
                data_disease_combinations_stats=pd.concat(data_disease_combinations_stats, axis=0)
                data_heterogeneity_measure.append(data_disease_combinations)
                data_heterogeneity_measure_stats.append(data_disease_combinations_stats)
                df=data_disease_combinations.copy()
            data_heterogeneity_measure=pd.concat(data_heterogeneity_measure, axis=0)
            data_heterogeneity_measure_stats=pd.concat(data_heterogeneity_measure_stats, axis=0)
            data_all_heterogeneity_measures.append(data_heterogeneity_measure)
            all_heterogeneity_measures_statistics.append(data_heterogeneity_measure_stats)
        data_all_heterogeneity_measures=pd.concat(data_all_heterogeneity_measures, axis=0)
        all_heterogeneity_measures_statistics=pd.concat(all_heterogeneity_measures_statistics, axis=0)
        pvalue_statistics.append(all_heterogeneity_measures_statistics)
    pvalue_statistics=pd.concat(pvalue_statistics, axis=0)
    pvalue_statistics.to_csv('../../results/pvalue_statistics_subsampled_patients.csv')
    
