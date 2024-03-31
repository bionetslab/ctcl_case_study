import pandas as pd
import os
import seaborn as sns
sns.set_theme(style="whitegrid")
import matplotlib.pyplot as plt
from _utilities_ import _substitute_value_in_dataframe_from_dict_
from _plots_ import _generate_violinplot_

# ========== Load and pre-process data required for generating plots: ==========
fig = plt.figure(figsize=(35, 30))
# subfigs = fig.subfigures(3, 1, hspace=0.07)
subfigs = fig.subfigures(1, 1)
subfigs_nested1=subfigs.subfigures(2, 1, wspace=0.00, height_ratios=[0.4, 0.6])
subfigs_nested2 = subfigs_nested1[0].subfigures(1, 2, wspace=0.00, width_ratios=[0.3, 0.4])

# ========== Create mosaic for: IA. Cell counts ==========
layout = [
    ["cell count"]
]
axes = subfigs_nested2[0].subplot_mosaic(layout)
subfigs_nested2[0].set_facecolor('0.85')
subfigs_nested2[0].suptitle('Cell counts', fontsize=40)

cell_counts={}
p_value_cell_counts=[]
for filename in os.listdir('../../cell_type_abundance_analyses/cell_count'):
    filename = os.fsdecode(filename)
    if filename.startswith('cell_counts_'):
        cell_count = pd.read_csv('../../cell_type_abundance_analyses/cell_count/'+filename)
        conditions_in_cell_count=list(set(cell_count.condition))
        for condition in conditions_in_cell_count:
            if not(condition in cell_counts.keys()):
                cell_counts[condition]=cell_count[cell_count['condition']==condition]
    elif filename.startswith('p_values_cell_counts_'):
        p_value_cell_count = pd.read_csv('../../cell_type_abundance_analyses/cell_count/'+filename)
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
title_prefix='A'
plot_title="All cell types"
# Generate violinplot:
subplot_axis_id=axes["cell count"]        
legend_loc='lower right'
# plot_title=''
data=cell_counts_df
data=data.reset_index(drop=True)
p_value=p_value_cell_counts_df
p_value=p_value.reset_index(drop=True)
_generate_violinplot_(data, conditions, 'condition', 'cell_count', p_value, subplot_axis_id, legend_loc, title_prefix=title_prefix, ylabel='Cell count', plot_title=plot_title)

# ========== Create mosaic for: IB. Celltype ratios ==========
layout = [
    ["celltype ratio: Endothelial cells", "celltype ratio: Endothelial cells", "celltype ratio: Macrophages", "celltype ratio: Macrophages", "celltype ratio: T-cells", "celltype ratio: T-cells", "celltype ratio: Fibroblasts", "celltype ratio: Fibroblasts"],
    ["celltype ratio: Endothelial cells", "celltype ratio: Endothelial cells", "celltype ratio: Macrophages", "celltype ratio: Macrophages", "celltype ratio: T-cells", "celltype ratio: T-cells", "celltype ratio: Fibroblasts", "celltype ratio: Fibroblasts"]
]
axes = subfigs_nested2[1].subplot_mosaic(layout)
subfigs_nested2[1].set_facecolor('0.95')
subfigs_nested2[1].suptitle('Cell type ratios', fontsize=40)

cell_type_ratios={}
p_value_cell_type_ratios=[]
for filename in os.listdir('../../cell_type_abundance_analyses/celltype_ratio'):
    filename = os.fsdecode(filename)
    if filename.startswith('cell_type_ratios_'):
        print(filename)
        cell_type_ratio = pd.read_csv('../../cell_type_abundance_analyses/celltype_ratio/'+filename)
        conditions_in_cell_type_ratio=list(set(cell_type_ratio.category))
        for condition in conditions_in_cell_type_ratio:
            if not(condition in cell_type_ratios.keys()):
                cell_type_ratios[condition]=cell_type_ratio[cell_type_ratio['category']==condition]
    elif filename.startswith('p_values_cell_type_ratios_'):
        p_value_cell_type_ratio = pd.read_csv('../../cell_type_abundance_analyses/celltype_ratio/'+filename)
        p_value_cell_type_ratios.append(p_value_cell_type_ratio)
p_value_cell_type_ratios_df=pd.concat(p_value_cell_type_ratios, axis=0)
conditions_abbreviations_lookup={'Eczema': 'AD', 'Psoriasis': 'PSO', 'T-Cell Lymphoma': 'CTCL'}
p_value_cell_type_ratios_df=_substitute_value_in_dataframe_from_dict_(p_value_cell_type_ratios_df, 'condition_1', conditions_abbreviations_lookup)
p_value_cell_type_ratios_df=_substitute_value_in_dataframe_from_dict_(p_value_cell_type_ratios_df, 'condition_2', conditions_abbreviations_lookup)
cell_type_ratios_df=[]
for condition, ratios_dataframe in cell_type_ratios.items():
    cell_type_ratios_df.append(ratios_dataframe)
cell_type_ratios_df=pd.concat(cell_type_ratios_df, axis=0)
cell_type_ratios_df=_substitute_value_in_dataframe_from_dict_(cell_type_ratios_df, 'category', conditions_abbreviations_lookup)
conditions=list(sorted(conditions_abbreviations_lookup.values()))

# ===== Generate celltype ratio plots for the most important celltypes: =====
celltypes=['Endothelial cells', 'Macrophages', 'T-cells', 'Fibroblasts']
title_prefix_='A'
for celltype in celltypes:
    # Generate violinplot:
    title_prefix=f'{title_prefix_}'
    subplot_axis_id=axes[f'celltype ratio: {celltype}']
    legend_loc='lower right'
    plot_title=celltype
    data=cell_type_ratios_df[cell_type_ratios_df['celltype']==celltype]
    data=data.reset_index(drop=True)
    p_value=p_value_cell_type_ratios_df[p_value_cell_type_ratios_df['celltype']==celltype]
    p_value=p_value.reset_index(drop=True)
    _generate_violinplot_(data, conditions, 'category', 'value', p_value, subplot_axis_id, legend_loc, title_prefix=title_prefix, plot_title=plot_title, ylabel='Celltype ratio')
    title_prefix_ = chr(ord(title_prefix_) + 1)



# ========== Create mosaic for: IA. Cell counts ==========
layout = [
    ["celltype count: Macrophages", "celltype count: Macrophages", "celltype count: Suprabasal keratinocytes", "celltype count: Suprabasal keratinocytes", "celltype count: Smooth muscle cells", "celltype count: Smooth muscle cells", "celltype count: Fibroblasts", "celltype count: Fibroblasts", "celltype count: T-cells", "celltype count: T-cells"],
    ["celltype count: Macrophages", "celltype count: Macrophages", "celltype count: Suprabasal keratinocytes", "celltype count: Suprabasal keratinocytes", "celltype count: Smooth muscle cells", "celltype count: Smooth muscle cells", "celltype count: Fibroblasts", "celltype count: Fibroblasts", "celltype count: T-cells", "celltype count: T-cells"],
    [".", "celltype count: Unknown", "celltype count: Unknown", "celltype count: B-cells", "celltype count: B-cells", "celltype count: Melanocytes", "celltype count: Melanocytes", "celltype count: Basal keratinocytes", "celltype count: Basal keratinocytes", "."],
    [".", "celltype count: Unknown", "celltype count: Unknown", "celltype count: B-cells", "celltype count: B-cells", "celltype count: Melanocytes", "celltype count: Melanocytes", "celltype count: Basal keratinocytes", "celltype count: Basal keratinocytes", "."]
]
axes = subfigs_nested1[1].subplot_mosaic(layout)
subfigs_nested1[1].set_facecolor('0.75')
subfigs_nested1[1].suptitle('Cell type counts', fontsize=40)

cell_type_counts={}
p_value_cell_type_counts=[]
for filename in os.listdir('../../cell_type_abundance_analyses/celltype_count'):
    filename = os.fsdecode(filename)
    if filename.startswith('cell_type_counts_'):
        cell_type_count = pd.read_csv('../../cell_type_abundance_analyses/celltype_count/'+filename)
        conditions_in_cell_type_count=list(set(cell_type_count.category))
        for condition in conditions_in_cell_type_count:
            if not(condition in cell_type_counts.keys()):
                cell_type_counts[condition]=cell_type_count[cell_type_count['category']==condition]
    elif filename.startswith('p_values_cell_type_counts_'):
        p_value_cell_type_count = pd.read_csv('../../cell_type_abundance_analyses/celltype_count/'+filename)
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

title_prefix_='A'
for celltype in celltypes:
    # Generate violinplot:
    subplot_axis_id=axes[f'celltype count: {celltype}']
    legend_loc='lower right'
    title_prefix=f'{title_prefix_}'
    plot_title=celltype
    data=cell_type_counts_df[cell_type_counts_df['celltype']==celltype]
    data=data.reset_index(drop=True)
    p_value=p_value_cell_type_counts_df[p_value_cell_type_counts_df['celltype']==celltype]
    p_value=p_value.reset_index(drop=True)
    _generate_violinplot_(data, conditions, 'category', 'value', p_value, subplot_axis_id, legend_loc, title_prefix=title_prefix, plot_title=plot_title, ylabel='Celltype count')
    
    title_prefix_ = chr(ord(title_prefix_) + 1)



# ========== Generate, save and show final plot (fig1): ==========
# plt.subplots_adjust(wspace=0.02, hspace=0.5)
plt.subplots_adjust(wspace=1.00, hspace=1.8)
plt.savefig('cell_abundance_analyses.pdf', format='pdf', bbox_inches='tight')
# fig.tight_layout(pad=100.0)
plt.show()

