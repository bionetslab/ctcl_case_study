import pandas as pd
import os
import seaborn as sns
sns.set_theme(style="whitegrid")
import matplotlib.pyplot as plt
from _utilities_ import _substitute_value_in_dataframe_from_dict_
from _generate_plots_for_celltype_ratio_analysis_ import _generate_violinplot_

# ========== Create mosaic for plot: ==========
palette={"AD":(0.12156862745098039, 0.4666666666666667, 0.7058823529411765),
         "PSO": (0.17254901960784313, 0.6274509803921569, 0.17254901960784313),
         "CTCL": (1.0, 0.4980392156862745, 0.054901960784313725)}
layout = [
    ["B-cells", "B-cells", "Basal keratinocytes", "Basal keratinocytes", "Endothelial cells", "Endothelial cells"],
    ["Fibroblasts", "Fibroblasts", "Langerhans cells", "Langerhans cells", "Macrophages", "Macrophages"],
    ["Melanocytes", "Melanocytes", "Smooth muscle cells", "Smooth muscle cells", "Suprabasal keratinocytes", "Suprabasal keratinocytes"],
    [".", "T-cells", "T-cells", "Unknown", "Unknown", "."]
]

fig, axes = plt.subplot_mosaic(layout, figsize=(25,25))

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
celltypes=sorted(list(set(cell_type_ratios_df.celltype)))

# ===== Generate celltype ratio plots for the most important celltypes: =====
conditions=['AD', 'PSO', 'CTCL']
for celltype in celltypes:
    # Generate violinplot:
    subplot_axis_id=axes[celltype]
    legend_loc='lower right'
    plot_title=celltype
    data=cell_type_ratios_df[cell_type_ratios_df['celltype']==celltype]
    data=data.reset_index(drop=True)
    p_value=p_value_cell_type_ratios_df[p_value_cell_type_ratios_df['celltype']==celltype]
    p_value=p_value.reset_index(drop=True)
    _generate_violinplot_(data, conditions, 'category', 'value', p_value, subplot_axis_id, legend_loc, plot_title=plot_title, ylabel='Celltype ratio', palette=palette)

# ========== Generate, save and show final plot (fig1): ==========
# plt.subplots_adjust(wspace=0.02, hspace=0.5)
plt.subplots_adjust(wspace=1.00, hspace=0.3)
plt.savefig('celltype_ratio_abundance_analyses_all_celltypes.pdf', format='pdf', bbox_inches='tight')
# fig.tight_layout(pad=100.0)
plt.show()

