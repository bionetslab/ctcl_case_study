import pandas as pd
import os
import seaborn as sns
sns.set_theme(style="whitegrid")
import matplotlib.pyplot as plt
from _generate_histplot_supfig1and2and3_shout_score_distributions_allConditions_ import _generate_histplot_

# ========== Load and pre-process data required for generating plots: ==========
cell_results=pd.read_csv(os.path.join('../../results', 'cell_results.csv'))
p_values_cell_type=pd.read_csv(os.path.join('../../results', 'p_values_cell_type.csv'))
conditions_=['Eczema', 'T-Cell Lymphoma', 'Psoriasis']
conditions=['AD', 'CTCL', 'PSO']
conditions_abbreviations_dict=dict(zip(conditions_, conditions))
cell_results=cell_results.replace({"condition": conditions_abbreviations_dict})
p_values_cell_type=p_values_cell_type.replace({"condition_1": conditions_abbreviations_dict})
p_values_cell_type=p_values_cell_type.replace({"condition_2": conditions_abbreviations_dict})
celltypes=sorted(list(set(cell_results.cell_type)))
# celltypes=[i for i in celltypes if i != 'T-cells']
# celltypes.append('T-cells')

# ========== Create mosaic for plot: ==========
layout = [
    ["B-cells", "B-cells", "Basal keratinocytes", "Basal keratinocytes", "Endothelial cells", "Endothelial cells"],
    ["Fibroblasts", "Fibroblasts", "Langerhans cells", "Langerhans cells", "Macrophages", "Macrophages"],
    ["Melanocytes", "Melanocytes", "Smooth muscle cells", "Smooth muscle cells", "Suprabasal keratinocytes", "Suprabasal keratinocytes"],
    [".", "T-cells", "T-cells", "Unknown", "Unknown", "."]
]

scores=['local_entropy_5', 'local_homophily_5', 'egophily_5']
conditions=['AD', 'PSO', 'CTCL']
hue='condition'
histplot_count=0
for score in scores:
    histplot_count+=1
    fig, axes = plt.subplot_mosaic(layout, figsize=(25,15))
    celltype_count=0
    for celltype in celltypes:
        celltype_count+=1
        if celltype_count==11:
            legend=True
        else:
            legend=None
        data=cell_results[cell_results['cell_type']==celltype]
        title=celltype
        subplot_axis_id=axes[celltype]
        data_=[]
        for condition in conditions:
            data_.append(data[data['condition']==condition])
        data=pd.concat(data_, axis=0)
        legend_loc="center right"
        _generate_histplot_(score, data, subplot_axis_id, hue, legend_loc, title, legend=legend)
    
    # ========== Generate, save and show final plot (fig1): ==========
    # plt.subplots_adjust(wspace=0.02, hspace=0.5)
    plt.subplots_adjust(wspace=0.60, hspace=1.20)
    plt.savefig(f'supfig{histplot_count}.pdf', format='pdf', bbox_inches='tight')
    # fig.tight_layout(pad=100.0)
    plt.show()
    


