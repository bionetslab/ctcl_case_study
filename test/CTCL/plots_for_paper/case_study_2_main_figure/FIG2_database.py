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
    ["PI", "PI", "PI", "PI", "PI", "PI"],
    ["PI", "PI", "PI", "PI", "PI", "PI"],
    # ---
    ["pca_celltype", "pca_celltype", "pca_celltype", "pca_condition", "pca_condition", "pca_condition"],
    ["pca_celltype", "pca_celltype", "pca_celltype", "pca_condition", "pca_condition", "pca_condition"],
    # ---
    ["pca_sex", "pca_sex", "pca_sex", "pca_age", "pca_age", "pca_age"],
    ["pca_sex", "pca_sex", "pca_sex", "pca_age", "pca_age", "pca_age"],
]

fig, axes = plt.subplot_mosaic(layout, figsize=(25,35))
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
    palette="dark", alpha=.6, # height=6
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

# # # =========== Plot for "CD45": ===========
# img = iio.imread("Propidium iodide.jpg")
# subplot_axis_id=axes["CD45"]
# title_loc='left'
# # plot_title='C. CD45, CD4'
# plot_title='C'
# _show_gray_image_(img, subplot_axis_id, cmap='gray', title_loc=title_loc, plot_title=plot_title)
# img = iio.imread("CD4.jpg")
# subplot_axis_id=axes["CD45"]
# _show_gray_image_(img, subplot_axis_id, alpha=0.75, legend=True, legend_labels=['CD45', 'CD4'])


# # =========== Plot for "pca_celltype": ===========
subplot_axis_id=axes["pca_celltype"]
hue='celltype'
title_loc='left'
# plot_title='D. PCA (colored by cell type)'
plot_title='C'
# Plot PCA using sns.scatterplot:
# _plot_pca_(pca_labeled, hue, subplot_axis_id, title_loc=title_loc, plot_title=plot_title)
# Plot PCA using sc.pl.pca:
_scanpy_umap_(adata, hue, subplot_axis_id, title_loc=title_loc, plot_title=plot_title) # sc.pl.pca(adata, color="celltype")

# # # =========== Plot for "pca_condition": ===========
subplot_axis_id=axes["pca_condition"]
hue='condition'
title_loc='left'
# plot_title='E. PCA (colored by condition)'
plot_title='D'
# Plot PCA using sns.scatterplot:
# _plot_pca_(pca_labeled, hue, subplot_axis_id, title_loc=title_loc, plot_title=plot_title)
# Plot PCA using sc.pl.pca:
_scanpy_umap_(adata, hue, subplot_axis_id, title_loc=title_loc, plot_title=plot_title) # sc.pl.pca(adata, color="condition")

# # # =========== Plot for "pca_sex": ===========
subplot_axis_id=axes["pca_sex"]
hue='sex'
title_loc='left'
# plot_title='F. PCA (colored by sex)'
plot_title='E'
palette = ["#F72585", "#4361EE"]
# Plot PCA using sns.scatterplot:
# _plot_pca_(pca_labeled, hue, subplot_axis_id, title_loc=title_loc, plot_title=plot_title, palette=palette)
# Plot PCA using sc.pl.pca:
_scanpy_umap_(adata, hue, subplot_axis_id, title_loc=title_loc, plot_title=plot_title) # sc.pl.pca(adata, color="sex")

# # # =========== Plot for "pca_age": ===========
# bin_dict={
#     1: '$20-29$',
#     2: '$30-39$',
#     3: '$40-49$',
#     4: '$50-59$',
#     5: '$60-69$',
#     6: '$70-79$',
#     7: '$80-89$'
# }
# subplot_axis_id=axes["pca_age"]
# hue='ages_binned_column'
# title_loc='left'
# plot_title='F. PCA (colored by age)'
# pca_labeled_=pca_labeled.copy()
# pca_labeled_=pca_labeled_.sort_values(by=['age'])
# ages_column=list(pca_labeled_.age)
# ages=sorted(list(set(pca_labeled_.age)))
# age_bins={}
# for age in ages:
#     bin=math.ceil((age-19)/10)
#     age_bins[age]=bin_dict[bin]
# ages_binned_column=[age_bins.get(item,item) for item in ages_column]
# pca_labeled_[hue]=ages_binned_column
# palette='Set1'
# alpha=0.25
# # Plot PCA using sns.scatterplot:
# # _plot_pca_(pca_labeled_, hue, subplot_axis_id, title_loc=title_loc, plot_title=plot_title, palette=palette, alpha=alpha)
# # Plot PCA using sc.pl.pca:
# _scanpy_pca_(adata, hue, subplot_axis_id, title_loc=title_loc, plot_title=plot_title) # sc.pl.pca(adata, color="age")

subplot_axis_id=axes["pca_age"]
hue='ages_binned_column'
title_loc='left'
# plot_title='G. PCA (colored by age)'
plot_title='F'
palette = 'Set1'
# Plot PCA using sns.scatterplot:
# _plot_pca_(pca_labeled, hue, subplot_axis_id, title_loc=title_loc, plot_title=plot_title, palette=palette)
# Plot PCA using sc.pl.pca:
_scanpy_umap_(adata, hue, subplot_axis_id, title_loc=title_loc, plot_title=plot_title) # sc.pl.pca(adata, color="sex")


# ========== Generate, save and show final plot (fig1): ==========
# plt.subplots_adjust(wspace=0.02, hspace=0.5)
plt.subplots_adjust(wspace=0.5, hspace=0.4)
plt.savefig('fig1_database.jpg', format='jpg', bbox_inches='tight')
# fig.tight_layout(pad=100.0)
plt.show()

