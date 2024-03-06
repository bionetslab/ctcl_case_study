import os
import sys
sys.path.insert(0, '../../')
from _utilities_ import _substitute_value_in_dataframe_from_dict_
import pandas as pd
import seaborn as sns
sns.set_theme(style="whitegrid")
import matplotlib.pyplot as plt
from _plots_ import _generate_violinplot_, _generate_histplot_

if __name__ == '__main__':
    cell_type_counts={}
    p_value_cell_type_counts=[]
    for filename in os.listdir('celltype_count'):
        filename = os.fsdecode(filename)
        if filename.startswith('cell_type_counts_'):
            cell_type_count = pd.read_csv('celltype_count/'+filename)
            conditions_in_cell_type_count=list(set(cell_type_count.category))
            for condition in conditions_in_cell_type_count:
                if not(condition in cell_type_counts.keys()):
                    cell_type_counts[condition]=cell_type_count[cell_type_count['category']==condition]
        elif filename.startswith('p_values_cell_type_counts_'):
            p_value_cell_type_count = pd.read_csv('celltype_count/'+filename)
            p_value_cell_type_counts.append(p_value_cell_type_count)
    p_value_cell_type_counts_df=pd.concat(p_value_cell_type_counts, axis=0)
    conditions_abbreviations_lookup={'Eczema': 'AD', 'Psoriasis': 'PSO', 'T-Cell Lymphoma': 'CTCL'}
    p_value_cell_type_counts_df=_substitute_value_in_dataframe_from_dict_(p_value_cell_type_counts_df, 'condition_1', conditions_abbreviations_lookup)
    p_value_cell_type_counts_df=_substitute_value_in_dataframe_from_dict_(p_value_cell_type_counts_df, 'condition_2', conditions_abbreviations_lookup)
    cell_type_counts_df=[]
    for condition, counts_dataframe in cell_type_counts.items():
        cell_type_counts_df.append(counts_dataframe)
    cell_type_counts_df=pd.concat(cell_type_counts_df, axis=0)
    cell_type_counts_df=_substitute_value_in_dataframe_from_dict_(cell_type_counts_df, 'category', conditions_abbreviations_lookup)
    conditions=list(sorted(conditions_abbreviations_lookup.values()))
    
    # ===== Generate celltype ratio plots for the most important celltypes: =====
    celltypes=[
        'Macrophages', 'Suprabasal keratinocytes',
        'Smooth muscle cells', 'Fibroblasts',
        'T-cells', 'Unknown', 'B-cells',
        'Melanocytes', 'Basal keratinocytes'
    ]
    # Generate mosaic:
    layout = [
        ["Macrophages 1", "Macrophages 2"],
        ["Suprabasal keratinocytes 1", "Suprabasal keratinocytes 2"],
        ["Smooth muscle cells 1", "Smooth muscle cells 2"],
        ["Fibroblasts 1", "Fibroblasts 2"],
        ["T-cells 1", "T-cells 2"],
        ["Unknown 1", "Unknown 2"],
        ["B-cells 1", "B-cells 2"],
        ["Melanocytes 1", "Melanocytes 2"],
        ["Basal keratinocytes 1", "Basal keratinocytes 2"]
    ]
    fig, axes = plt.subplot_mosaic(layout, layout='constrained', figsize=(10,25)) # celltype ratio figsize = (20, 15)
    title_prefix_='A'
    for celltype in celltypes:
        # Generate violinplot:
        subplot_axis_id=axes[f'{celltype} 1']
        legend_loc='lower right'
        title_prefix=f'{title_prefix_}1. '
        plot_title=celltype
        data=cell_type_counts_df[cell_type_counts_df['celltype']==celltype]
        data=data.reset_index(drop=True)
        p_value=p_value_cell_type_counts_df[p_value_cell_type_counts_df['celltype']==celltype]
        p_value=p_value.reset_index(drop=True)
        _generate_violinplot_(data, conditions, 'category', 'value', p_value, subplot_axis_id, legend_loc, title_prefix=title_prefix, plot_title=plot_title, ylabel='Celltype count')
        
        # Generate histplot:
        subplot_axis_id=axes[f'{celltype} 2']
        legend_loc='upper right'
        title_prefix=f'{title_prefix_}2. '
        plot_title=celltype
        data=cell_type_counts_df[cell_type_counts_df['celltype']==celltype]
        data=data.reset_index(drop=True)
        _generate_histplot_('value', data, subplot_axis_id, 'category', legend_loc, title_prefix=title_prefix, plot_title=plot_title, xlabel='Celltype ratio')
        title_prefix_ = chr(ord(title_prefix_) + 1)
    
    # ========== Generate, save and show final plot: ==========
    plt.savefig('celltype_count_abundance_allConditions.pdf', format='pdf', bbox_inches='tight')
    fig.tight_layout()
    plt.show()
        
        
        
        
        
        
        
        
        
        
    
    
    