import pandas as pd
import scipy.stats as stats
import itertools as itt

def _run_differential_analyses_on_dict_(dictionary, categories, key_category_lookup):
    p_values_allCombinations={}
    data_allCombinations={}
    category_combinations=list(itt.combinations(categories, 2))
    for combination in category_combinations:
        p_values_all_keys=[]
        data_all_keys=[]
        for cell, list_of_ratios in dictionary.items():
            data=pd.DataFrame()
            data['celltype']=[cell]*len(list_of_ratios)
            data['value']=list_of_ratios
            data['category']=key_category_lookup[cell]
            data=pd.concat([data[data['category']==combination[0]], data[data['category']==combination[1]]], axis=0)
            values_0=list(data[data['category']==combination[0]]['value'])
            values_1=list(data[data['category']==combination[1]]['value'])
            U, p = stats.mannwhitneyu(x=values_0, y=values_1)
            p_value=pd.DataFrame()
            p_value['celltype']=[cell]
            p_value['condition_1']=[combination[0]]
            p_value['condition_2']=[combination[1]]
            p_value['U']=[U]
            p_value['p']=[p]
            p_values_all_keys.append(p_value)
            data_all_keys.append(data)
        p_values_all_keys=pd.concat(p_values_all_keys, axis=0)
        p_values_all_keys=p_values_all_keys.sort_values(by=['p'])
        p_values_all_keys=p_values_all_keys.reset_index(drop=True)
        data_all_keys=pd.concat(data_all_keys, axis=0)
        data_all_keys=data_all_keys.reset_index(drop=True)
        p_values_allCombinations[combination]=p_values_all_keys
        data_allCombinations[combination]=data_all_keys
    return data_allCombinations, p_values_allCombinations


def _run_differential_analyses_on_dataframe_(dataframe, categories, category_name, value_name):
    p_values_allCombinations=[]
    data_allCombinations={}
    category_combinations=list(itt.combinations(categories, 2))
    for combination in category_combinations:
        data_1=dataframe[dataframe[category_name]==combination[0]]
        data_2=dataframe[dataframe[category_name]==combination[1]]
        data=pd.concat([data_1, data_2], axis=0)
        values_0=list(data_1[value_name])
        values_1=list(data_2[value_name])
        U, p = stats.mannwhitneyu(x=values_0, y=values_1)
        p_value=pd.DataFrame()
        p_value['combination']=[str(combination)]
        p_value['U']=[U]
        p_value['p']=[p]
        print(p_value)
        p_values_allCombinations.append(p_value)
        data=data.reset_index(drop=True)
        data_allCombinations[combination]=data
    p_values_allCombinations=pd.concat(p_values_allCombinations, axis=0)
    p_values_allCombinations=p_values_allCombinations.reset_index(drop=True)
    return data_allCombinations, p_values_allCombinations
    


