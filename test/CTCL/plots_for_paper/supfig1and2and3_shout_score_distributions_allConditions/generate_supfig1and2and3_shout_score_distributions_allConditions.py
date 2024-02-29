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
celltypes=list(set(cell_results.cell_type))
celltypes=[i for i in celltypes if i != 'T-cells']
celltypes.append('T-cells')
# ========== Create mosaic for plot: ==========


layout = [
    ["Basal keratinocytes", "Smooth muscle cells", "Melanocytes"],
    ["Unknown", "B-cells", "Fibroblasts"],
    [ "Suprabasal keratinocytes", "Endothelial cells", "Macrophages"],
    ["Langerhans cells", "T-cells", ""]
]


scores=['local_entropy_3', 'local_homophily_3', 'egophily_3']
conditions=['AD', 'CTCL', 'PSO']
hue='condition'
histplot_count=0
for score in scores:
    histplot_count+=1
    fig, axes = plt.subplot_mosaic(layout, layout='constrained', figsize=(20,10))
    subplot_str='A'
    for celltype in celltypes:
        data=cell_results[cell_results['cell_type']==celltype]
        title=subplot_str+f'. {celltype}'
        subplot_axis_id=axes[celltype]
        data_=[]
        for condition in conditions:
            data_.append(data[data['condition']==condition])
        data=pd.concat(data_, axis=0)
        if score=='local_entropy_3':
            legend_loc="upper left"
        elif score=='local_homophily_3':
            legend_loc="upper right"
        elif score=='egophily_3':
            legend_loc="upper right"
        _generate_histplot_(score, data, subplot_axis_id, hue, legend_loc, title)
        subplot_str=chr(ord(subplot_str)+1)
    suptitle=f'Distributions of {score} scores for AD, CTCL and PSO samples'
    plt.suptitle(suptitle)
    plt.savefig(f'supfig{histplot_count}.pdf', format='pdf', bbox_inches='tight')
    fig.tight_layout()
    plt.show()
    


