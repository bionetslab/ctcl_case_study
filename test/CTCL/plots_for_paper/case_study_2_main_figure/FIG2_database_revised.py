import pandas as pd
import seaborn as sns
sns.set_theme(style="whitegrid")
import matplotlib.pyplot as plt
import numpy as np
import os
import imageio as iio
import scanpy as sc
from _show_gray_image_ import _show_gray_image_
from _plot_pca_ import _plot_pca_, _scanpy_pca_, _scanpy_umap_
from _utilities_ import _substitute_value_in_dataframe_from_dict_

# ========== Load and pre-process data required for generating plots: ==========
pca_labeled = pd.read_csv('../../hpa_based_cell_type_assignment/pca_labeled.csv')
adata = sc.read_h5ad('../../hpa_based_cell_type_assignment/celltype_assigned_anndata.h5ad')

palette={"AD":(0.12156862745098039, 0.4666666666666667, 0.7058823529411765),
         "PSO": (0.17254901960784313, 0.6274509803921569, 0.17254901960784313),
         "CTCL": (1.0, 0.4980392156862745, 0.054901960784313725)}
celltypes_lookup={
        'B-cells': 'B cells',
        'Basal keratinocytes': 'Basal kers',
        'Endothelial cells': 'Endothelial',
        'Fibroblasts': 'Fibroblasts',
        'Langerhans cells': 'Langerhans',
        'Macrophages': 'Macrophages',
        'Melanocytes': 'Melanocytes',
        'Smooth muscle cells': 'SMCs',
        'Suprabasal keratinocytes': 'Supra kers',
        'T-cells': 'T cells',
        'Unknown': 'Unknown'
    }
# # ========== Create mosaic for "celltype_ratios": ==========
layout = [
    ["celltype_ratios", "celltype_ratios", "celltype_ratios", "celltype_ratios", "celltype_ratios", "celltype_ratios"],
    ["celltype_ratios", "celltype_ratios", "celltype_ratios", "celltype_ratios", "celltype_ratios", "celltype_ratios"],
    # ---
    ["PI", "PI", "PI", "PI&CD4", "PI&CD4", "PI&CD4"],
    ["PI", "PI", "PI", "PI&CD4", "PI&CD4", "PI&CD4"]
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
conditions_abbreviations_lookup={'Eczema': 'AD', 'Psoriasis': 'PSO', 'T-Cell Lymphoma': 'CTCL'}
cell_type_ratios_df=[]
for condition, ratios_dataframe in cell_type_ratios.items():
    cell_type_ratios_df.append(ratios_dataframe)
cell_type_ratios_df=pd.concat(cell_type_ratios_df, axis=0)
cell_type_ratios_df=_substitute_value_in_dataframe_from_dict_(cell_type_ratios_df, 'category', conditions_abbreviations_lookup)
# conditions=list(sorted(conditions_abbreviations_lookup.values()))
conditions=['AD', 'PSO', 'CTCL']
celltypes=sorted(list(set(cell_type_ratios_df.celltype)))

df=[]
for condition in conditions:
    df_=[]
    cell_type_ratios_df_filtered_by_condition=cell_type_ratios_df[cell_type_ratios_df.category==condition]
    for celltype in celltypes:
        cell_type_ratios_df_filtered_by_celltype=cell_type_ratios_df_filtered_by_condition[cell_type_ratios_df_filtered_by_condition.celltype==celltype]
        _df_=pd.DataFrame()
        _df_['celltype']=[celltype]
        _df_['category']=[condition]
        _df_['mean_values']=np.mean(cell_type_ratios_df_filtered_by_celltype.value)
        df_.append(_df_)
    df_=pd.concat(df_, axis=0)
    df.append(df_)
df=pd.concat(df, axis=0)


ax=sns.barplot(
    ax=axes["celltype_ratios"],
    data=df, # kind="bar",
    x="celltype", y="mean_values", hue="category",
    palette=palette
)

# plot_title='A. Mean relative expression'
plot_title='A'
title_loc='left'
ax.set_title(plot_title, loc=title_loc, fontsize=30, fontweight='bold') # fontsize=25, pad=20
ax.tick_params(axis='y', which='major', labelsize=20)
xticks=[]
for celltype, abbr in celltypes_lookup.items():
    xticks.append(abbr)
ax.set_xticklabels(xticks, size = 20)
# ylabel='$\overline{(\\frac{n_i}{N})}$'
ylabel='Mean cell type fraction'
ax.set_ylabel(ylabel, fontsize=35, labelpad=10)
ax.set_xlabel(None, fontsize=30, labelpad=5)

sns.move_legend(
    ax,
    "upper right",
    # bbox_to_anchor=[1.090, -0.375],
    ncol=3,
    title=None,
    frameon=True,
    markerscale=3,
    fontsize=25
)

# =========== Plot for "PI": ===========
img = iio.imread("Propidium iodide.jpg")
subplot_axis_id=axes["PI"]
title_loc='left'
# plot_title='B. Propidium iodide, CD4'
plot_title='B'
_show_gray_image_(img, subplot_axis_id, cmap='gray', title_loc=title_loc, plot_title=plot_title)
# img = iio.imread("CD4.jpg")
# subplot_axis_id=axes["PI"]
# _show_gray_image_(img, subplot_axis_id, alpha=0.75, legend=True, legend_labels=['Propidium iodide', 'CD4'])

# # # =========== Plot for "PI&CD4": ===========
img = iio.imread("Propidium iodide.jpg")
subplot_axis_id=axes["PI&CD4"]
title_loc='left'
# plot_title='B. Propidium iodide, CD4'
plot_title='C'
_show_gray_image_(img, subplot_axis_id, cmap='gray', title_loc=title_loc, plot_title=plot_title)
img = iio.imread("CD4.jpg")
subplot_axis_id=axes["PI&CD4"]
_show_gray_image_(img, subplot_axis_id, alpha=0.75, legend=True, legend_labels=['Propidium iodide', 'CD4'])

# ========== Generate, save and show final plot (fig1): ==========
# plt.subplots_adjust(wspace=0.02, hspace=0.5)
plt.subplots_adjust(wspace=0.5, hspace=0.4)
plt.savefig('fig1_database.pdf', format='pdf', bbox_inches='tight')
# fig.tight_layout(pad=100.0)
plt.show()

