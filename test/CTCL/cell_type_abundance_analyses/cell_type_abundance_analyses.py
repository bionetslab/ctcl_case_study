import os
import scanpy as sc
import sys
sys.path.insert(0, '../../')
from _utilities_ import _transfer_values_from_local_to_global_dict_, _save_dict_of_dataframes_to_csv_files_
from _statistical_tests_ import _run_differential_analyses_on_dict_, _run_differential_analyses_on_dataframe_
import pandas as pd

if __name__ == '__main__':
    # Prepare spatial data for running the algorithm:
    cell_type_counts_allSamples={}
    cell_type_ratios_allSamples={}
    celltypes_conditions_dict={} # We need this because there is a possibility that not all samples contain all types of cells. 
    cell_counts_allSamples=[]
    list_of_conditions_across_samples=[]
    for filename in os.listdir('../data'):
        filename = os.fsdecode(filename)
        if filename.endswith('.h5ad'):
            cell_types_allSamples=list(cell_type_counts_allSamples.keys())
            adata = sc.read_h5ad('../data/'+filename)
            sample_condition=adata.uns['Group'][0]
            list_of_conditions_across_samples.append(sample_condition)
            cell_types_in_sample=list(set(adata.obs['celltype']))
            cells_in_sample=list(adata.obs['celltype'])
            cell_type_counts_in_sample = {}
            cell_type_ratios_in_sample = {}
            condition_in_sample={}
            cell_counts_in_sample = len(adata.to_df())
            cell_counts_allSamples.append(cell_counts_in_sample)
            for cell in cells_in_sample:
                condition_in_sample[cell] = sample_condition
                cell_type_counts_in_sample[cell] = cells_in_sample.count(cell)
                cell_type_ratios_in_sample[cell] = cells_in_sample.count(cell)/cell_counts_in_sample
            celltype_counts_allSamples=_transfer_values_from_local_to_global_dict_(cell_type_counts_in_sample, cell_type_counts_allSamples)
            celltype_ratios_allSamples=_transfer_values_from_local_to_global_dict_(cell_type_ratios_in_sample, cell_type_ratios_allSamples)
            celltypes_conditions_dict=_transfer_values_from_local_to_global_dict_(condition_in_sample, celltypes_conditions_dict)
    cell_counts_allSamples_df=pd.DataFrame()
    cell_counts_allSamples_df['cell_count']=cell_counts_allSamples
    cell_counts_allSamples_df['condition']=list_of_conditions_across_samples
    conditions=list(set(list_of_conditions_across_samples))
    condition_abbreviations_dict={'Eczema': 'AD', 'Psoriasis': 'PSO', 'T-Cell Lymphoma': 'CTCL'}
    list_of_conditions_across_samples=[condition_abbreviations_dict.get(item,item)  for item in list_of_conditions_across_samples]
    cell_type_ratios_allCombinations, p_values_cell_type_ratios_allCombinations=_run_differential_analyses_on_dict_(cell_type_ratios_allSamples, conditions, celltypes_conditions_dict)
    cell_type_counts_allCombinations, p_values_cell_type_counts_allCombinations=_run_differential_analyses_on_dict_(cell_type_counts_allSamples, conditions, celltypes_conditions_dict)
    cell_counts_allCombinations, p_values_cell_counts_allCombinations=_run_differential_analyses_on_dataframe_(cell_counts_allSamples_df, conditions, 'condition', 'cell_count')
    dicts_to_be_saved_as_csv=[cell_type_ratios_allCombinations, p_values_cell_type_ratios_allCombinations, cell_type_counts_allCombinations, p_values_cell_type_counts_allCombinations, cell_counts_allCombinations]
    filename_prefixes=['cell_type_ratios', 'p_values_cell_type_ratios', 'cell_type_counts', 'p_values_cell_type_counts', 'cell_counts'] 
    filenamePrefixes_dicts=dict(zip(filename_prefixes, dicts_to_be_saved_as_csv))
    for filename_prefix, dictionary in filenamePrefixes_dicts.items():
        _save_dict_of_dataframes_to_csv_files_(dictionary, filename_prefix)
    p_values_cell_counts_allCombinations.to_csv('p_values_cell_counts_allCombinations.csv')


