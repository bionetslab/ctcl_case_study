import os
import sys
sys.path.insert(0, '../../')
from _utilities_ import _substitute_value_in_dataframe_from_dict_
import pandas as pd
import seaborn as sns
sns.set_theme(style="whitegrid")
import matplotlib.pyplot as plt
from _plots_ import _generate_violinplot_, _generate_histplot_
import itertools as itt
from decimal import Decimal
from statannotations.Annotator import Annotator

if __name__ == '__main__':
    cell_counts={}
    p_value_cell_counts=[]
    for filename in os.listdir('cell_count'):
        filename = os.fsdecode(filename)
        if filename.startswith('cell_counts_'):
            cell_count = pd.read_csv('cell_count/'+filename)
            conditions_in_cell_count=list(set(cell_count.condition))
            for condition in conditions_in_cell_count:
                if not(condition in cell_counts.keys()):
                    cell_counts[condition]=cell_count[cell_count['condition']==condition]
        elif filename.startswith('p_values_cell_counts_'):
            p_value_cell_count = pd.read_csv('cell_count/'+filename)
            p_value_cell_counts.append(p_value_cell_count)
    p_value_cell_counts_df=pd.concat(p_value_cell_counts, axis=0)
    conditions_abbreviations_lookup={'Eczema': 'AD', 'Psoriasis': 'PSO', 'T-Cell Lymphoma': 'CTCL'}
    # p_value_cell_counts_df=_substitute_value_in_dataframe_from_dict_(p_value_cell_counts_df, 'condition', conditions_abbreviations_lookup)
    # p_value_cell_counts_df=_substitute_value_in_dataframe_from_dict_(p_value_cell_counts_df, 'condition_2', conditions_abbreviations_lookup)
    pvalue_combinations_list=list(p_value_cell_counts_df['combination'])
    new_list=[]
    for i in pvalue_combinations_list:
        exec('element=('+i+'[0], '+i+'[1])')
        new_list.append(element)
    df=pd.DataFrame(new_list, columns =['condition_1', 'condition_2'])
    p_value_cell_counts_df=p_value_cell_counts_df.drop(columns=['combination'])
    p_value_cell_counts_df=pd.concat([p_value_cell_counts_df, df], axis=1)
    p_value_cell_counts_df=_substitute_value_in_dataframe_from_dict_(p_value_cell_counts_df, 'condition_1', conditions_abbreviations_lookup)
    p_value_cell_counts_df=_substitute_value_in_dataframe_from_dict_(p_value_cell_counts_df, 'condition_2', conditions_abbreviations_lookup)
    cell_counts_df=[]
    for condition, ratios_dataframe in cell_counts.items():
        cell_counts_df.append(ratios_dataframe)
    cell_counts_df=pd.concat(cell_counts_df, axis=0)
    cell_counts_df=_substitute_value_in_dataframe_from_dict_(cell_counts_df, 'condition', conditions_abbreviations_lookup)
    conditions=list(sorted(conditions_abbreviations_lookup.values()))
    
    # ===== Generate cell count plots for all celltypes: =====
    # Generate mosaic:
    layout = [
        ["X1", "X2"]
    ]
    fig, axes = plt.subplot_mosaic(layout, layout='constrained', figsize=(25,12.5))
    title_prefix_='A'
    # Generate violinplot:
    subplot_axis_id=axes['X1']        
    legend_loc='lower right'
    title_prefix=f'{title_prefix_}1. '
    # plot_title=''
    data=cell_counts_df
    data=data.reset_index(drop=True)
    p_value=p_value_cell_counts_df
    p_value=p_value.reset_index(drop=True)
    _generate_violinplot_(data, conditions, 'condition', 'cell_count', p_value, subplot_axis_id, legend_loc, title_prefix=title_prefix, ylabel='Cell count')
    
    # Generate histplot:
    subplot_axis_id=axes['X2']
    legend_loc='upper right'
    title_prefix=f'{title_prefix_}2. '
    # plot_title=''
    data=cell_counts_df
    data=data.reset_index(drop=True)
    _generate_histplot_('cell_count', data, subplot_axis_id, 'condition', legend_loc, title_prefix=title_prefix, xlabel='Cell count')
    title_prefix_ = chr(ord(title_prefix_) + 1)
    
    # ========== Generate, save and show final plot: ==========
    plt.savefig('cell_count_abundance_allConditions.pdf', format='pdf', bbox_inches='tight')
    fig.tight_layout()
    plt.show()
        
        
        
        
        
        
        
        
        
        
    
    
    