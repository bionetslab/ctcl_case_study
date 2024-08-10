import pandas as pd
import os
import seaborn as sns
sns.set_theme(style="whitegrid")
import matplotlib.pyplot as plt
from _plot_distributions_per_individual_celltypes_ import _plot_distributions_per_individual_celltypes_
from _plot_bk_neighborhood_enrichment_ import _plot_bk_neighborhood_enrichment_
from _plot_bk_centrality_measures_ import _plot_bk_centrality_measures_
from _generate_scatter_subplot_ import _generate_scatter_subplot_
from _plot_global_score_distributions_ import _plot_global_score_distributions_
import scanpy as sc

abbr={
              'B-cells': 'B',
              'Basal keratinocytes': 'BK',
              'Endothelial cells': 'En',
              'Fibroblasts': 'Fi',
              'Langerhans cells': 'La',
              'Macrophages': 'Ma',
              'Melanocytes': 'Me',
              'Smooth muscle cells': 'SMC',
              'Suprabasal keratinocytes': 'SK',
              'T-cells': 'T',
              'Unknown': 'Un'
      }

textstr = (r'$\bf{B}$: B-cells                                  '+
    r'$\bf{BK}$: Basal keratinocytes       '+
    r'$\bf{En}$: Endothelial cells              '+
    r'$\bf{Fi}$: Fibroblasts'+
    '\n'+
    r'$\bf{La}$: Langerhans cells                '+
    r'$\bf{Ma}$: Macrophages                '+
    r'$\bf{Me}$: Melanocytes                  '+
    r'$\bf{SMC}$: Smooth muscle cells'+
    '\n'+
    r'$\bf{SK}$: Suprabasal keratinocytes  '+
    r'$\bf{T}$: T cells                               '+
    r'$\bf{Un}$: Unknown')

props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

# ========== Load and pre-process data required for generating plots: ==========
cell_results=pd.read_csv(os.path.join('../../../results', 'cell_results.csv'))
p_values_cell_type=pd.read_csv(os.path.join('../../../results', 'p_values_cell_type.csv'))
squidpy_nhoodEnrichment_results=pd.read_csv(os.path.join('../../../results', 'squidpy_nhoodEnrichment_results.csv'))
p_values_cell_type_squidpy_nhoodEnrichment=pd.read_csv(os.path.join('../../../results', 'p_values_cell_type_squidpy_nhoodEnrichment.csv'))
squidpy_centralityScores_results=pd.read_csv(os.path.join('../../../results', 'squidpy_centralityScores_results.csv'))
p_values_cell_type_squidpy_centralityScores=pd.read_csv(os.path.join('../../../results', 'p_values_cell_type_squidpy_centralityScores.csv'))
conditions_=['Eczema', 'T-Cell Lymphoma', 'Psoriasis']
conditions=['AD', 'CTCL', 'PSO']
conditions_abbreviations_dict=dict(zip(conditions_, conditions))
cell_results=cell_results.replace({"condition": conditions_abbreviations_dict})
p_values_cell_type=p_values_cell_type.replace({"condition_1": conditions_abbreviations_dict})
p_values_cell_type=p_values_cell_type.replace({"condition_2": conditions_abbreviations_dict})
squidpy_nhoodEnrichment_results=squidpy_nhoodEnrichment_results.replace({"condition": conditions_abbreviations_dict})
p_values_cell_type_squidpy_nhoodEnrichment=p_values_cell_type_squidpy_nhoodEnrichment.replace({"condition_1": conditions_abbreviations_dict})
p_values_cell_type_squidpy_nhoodEnrichment=p_values_cell_type_squidpy_nhoodEnrichment.replace({"condition_2": conditions_abbreviations_dict})
squidpy_centralityScores_results=squidpy_centralityScores_results.replace({"condition": conditions_abbreviations_dict})
p_values_cell_type_squidpy_centralityScores=p_values_cell_type_squidpy_centralityScores.replace({"condition_1": conditions_abbreviations_dict})
p_values_cell_type_squidpy_centralityScores=p_values_cell_type_squidpy_centralityScores.replace({"condition_2": conditions_abbreviations_dict})

palette={"AD":(0.12156862745098039, 0.4666666666666667, 0.7058823529411765),
         "PSO": (0.17254901960784313, 0.6274509803921569, 0.17254901960784313),
         "CTCL": (1.0, 0.4980392156862745, 0.054901960784313725)}
fig = plt.figure(figsize=(25, 25))
# subfigs = fig.subfigures(3, 1, hspace=0.07)
subfigs = fig.subfigures(3, 1, hspace=0.0, height_ratios=[1, 2, 2])

# ========== Create mosaic for: I. T-cell plots ==========
layout = [
    ["A", "A", "B", "B", "C", "C"],
    ["A", "A", "B", "B", "C", "C"],
    ["D1", "D1", "D2", "D2", "D3", "D3"],
    ["D1", "D1", "D2", "D2", "D3", "D3"],
]
axes = subfigs[1].subplot_mosaic(layout)
# subfigs[1].set_facecolor('0.85')
subfigs[1].suptitle(r'$\bf{C}$ Local heterogeneity, centrality and neighborhood enrichment scores (T cells)', fontsize=30, x=0.3)

# ========== Data for plot A: ==========
heterogeneity_measure='entropy'
celltype='T-cells'
conditions=['AD', 'PSO', 'CTCL']
radii=[1, 5]
subplot_axis_id=axes["A"]
legend_loc=None
title_prefix=None # title_prefix='A '
title_loc='left'
plot_title=r'$\bf{C1}$'
# Generate plot-1:
_plot_distributions_per_individual_celltypes_(p_values_cell_type, cell_results, celltype, heterogeneity_measure, conditions, radii, subplot_axis_id, legend_loc=legend_loc, title_prefix=title_prefix, title_loc=title_loc, palette=palette, plot_title=plot_title)

# ========== Data for plot B: ==========
heterogeneity_measure='egophily'
celltype='T-cells'
conditions=['AD', 'PSO', 'CTCL']
radii=[1, 5]
subplot_axis_id=axes["B"]
legend_loc=None
title_prefix=None # title_prefix='B '
title_loc='left'
plot_title=r'$\bf{C2}$'
# Generate plot-B:
_plot_distributions_per_individual_celltypes_(p_values_cell_type, cell_results, celltype, heterogeneity_measure, conditions, radii, subplot_axis_id, legend_loc=legend_loc, title_prefix=title_prefix, title_loc=title_loc, palette=palette,plot_title=plot_title)

# ========== Data for plot C: ==========
celltype='T-cells'
conditions=['AD', 'PSO', 'CTCL']
title_prefix=None # title_prefix='C '
scores = ['degree_centrality', 'closeness_centrality']
data=squidpy_centralityScores_results[squidpy_centralityScores_results['celltypes']==celltype]
data_=[]
for condition in conditions:
    data_.append(data[data['condition']==condition])
data=pd.concat(data_, axis=0)
columns=scores+['condition']
data=data[columns]
subplot_axis_id=axes["C"]
title_loc="left"
plot_title=r'$\bf{C3}$'
legend_loc="upper right"
ylabel='Centrality scores'
bbox_to_anchor=(1.3, 1.00)
_plot_bk_centrality_measures_(p_values_cell_type_squidpy_centralityScores, data, celltype, conditions, scores, subplot_axis_id, legend_loc=legend_loc, title_prefix=title_prefix, title_loc=title_loc, ylabel=ylabel, palette=palette, bbox_to_anchor=bbox_to_anchor, plot_title=plot_title)

# ========== Data for plot D: ==========
reference_celltype='T-cells'
# celltypes_=['Smooth muscle cells', 'Macrophages', 'Endothelial cells', 'B-cells', 'Fibroblasts', 'Basal keratinocytes']
celltypes_=['Macrophages', 'Fibroblasts', 'Basal keratinocytes']
conditions=['AD', 'PSO', 'CTCL']
cnt=0
for celltype in celltypes_:
    cnt+=1
    title_prefix=None
    suptitle=None
    if cnt==1:
        # title_prefix='D '
        # suptitle='\t\t      '
        ylabel='Neighborhood\nenrichment\n(z-scores)'
        if celltype=='Suprabasal keratinocytes':
            plot_title=r'\bf{C4}\n                         T vs. SK'
        else:
            plot_title=str(r'$\bf{C4}$'+f'\n                         T vs. {abbr[celltype]}')
        title_loc="left"
    else:
        # title_prefix=None
        # suptitle=None
        ylabel=None
        if celltype=='Suprabasal keratinocytes':
            plot_title='T vs. suprabasal keratinocytes'
        else:
            plot_title=str(f'T vs. {abbr[celltype]}')
        title_loc="center"
    # if celltype=='Suprabasal keratinocytes':
    #     plot_title='T cells vs. suprabasal keratinocytes'
    # else:
    #     plot_title=str(f'T cells vs. {celltype.lower()}')
    data=pd.concat([squidpy_nhoodEnrichment_results[squidpy_nhoodEnrichment_results['celltype_1']==celltype][squidpy_nhoodEnrichment_results['celltype_2']=='T-cells'], squidpy_nhoodEnrichment_results[squidpy_nhoodEnrichment_results['celltype_2']==celltype][squidpy_nhoodEnrichment_results['celltype_1']=='T-cells']], axis=0)
    subplot_axis_id=axes[f"D{cnt}"]
    data_=[]
    for condition in conditions:
        data_.append(data[data['condition']==condition])
    data=pd.concat(data_, axis=0)
    _plot_bk_neighborhood_enrichment_(p_values_cell_type_squidpy_nhoodEnrichment, squidpy_nhoodEnrichment_results, data, celltype, reference_celltype, conditions, subplot_axis_id, title_prefix=title_prefix, plot_title=plot_title, suptitle=suptitle, ylabel=ylabel, palette=palette,title_loc=title_loc)
subfigs[1].subplots_adjust(wspace=0.45, hspace=1.3)

# ========== Create mosaic for: II. Basal keratinocyte plots ==========
layout = [
    ["A", "A", "A", "A", "A", "A", "A", "A", "A", "A", ".", ".", "B", "B", "B", "B", "B", "B", "B", "B", "B", "B", ".", ".", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C"],
    ["A", "A", "A", "A", "A", "A", "A", "A", "A", "A", ".", ".", "B", "B", "B", "B", "B", "B", "B", "B", "B", "B", ".", ".", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C"],
    ["D1", "D1", "D1", "D1", "D1", "D1", ".", "D2", "D2", "D2", "D2", "D2", "D2", ".", "D3", "D3", "D3", "D3", "D3", "D3", ".", "D4", "D4", "D4", "D4", "D4", "D4", ".", "D5", "D5", "D5", "D5", "D5", "D5"],
    ["D1", "D1", "D1", "D1", "D1", "D1", ".", "D2", "D2", "D2", "D2", "D2", "D2", ".", "D3", "D3", "D3", "D3", "D3", "D3", ".", "D4", "D4", "D4", "D4", "D4", "D4", ".", "D5", "D5", "D5", "D5", "D5", "D5"],
]
axes = subfigs[2].subplot_mosaic(layout)
subfigs[2].set_facecolor('1.00')
subfigs[2].suptitle(r'$\bf{D}$ Local heterogeneity, centrality and neighborhood enrichment scores (basal keratinocytes)', fontsize=30, x=0.35)

# ========== Data for plot A: ==========
heterogeneity_measure='entropy'
celltype='Basal keratinocytes'
conditions=['AD', 'PSO', 'CTCL']
radii=[1, 5]
subplot_axis_id=axes["A"]
legend_loc='lower right'
# title_prefix='A '
title_prefix=None
title_loc='left'
plot_title=r'$\bf{D1}$'
# Generate plot-1:
_plot_distributions_per_individual_celltypes_(p_values_cell_type, cell_results, celltype, heterogeneity_measure, conditions, radii, subplot_axis_id, legend_loc=legend_loc, title_prefix=title_prefix, title_loc=title_loc, palette=palette, plot_title=plot_title)

# ========== Data for plot B: ==========
heterogeneity_measure='homophily'
celltype='Basal keratinocytes'
conditions=['AD', 'PSO', 'CTCL']
radii=[1, 5]
subplot_axis_id=axes["B"]
legend_loc='upper left'
# title_prefix='B '
title_prefix=None
title_loc='left'
plot_title=r'$\bf{D2}$'
# Generate plot-B:
_plot_distributions_per_individual_celltypes_(p_values_cell_type, cell_results, celltype, heterogeneity_measure, conditions, radii, subplot_axis_id, legend_loc=legend_loc, title_prefix=title_prefix, title_loc=title_loc, palette=palette, plot_title=plot_title)

# ========== Data for plot C: ==========
celltype='Basal keratinocytes'
conditions=['AD', 'PSO', 'CTCL']
# title_prefix='C '
title_prefix=None
scores = ['degree_centrality', 'closeness_centrality']
data=squidpy_centralityScores_results[squidpy_centralityScores_results['celltypes']==celltype]
data_=[]
for condition in conditions:
    data_.append(data[data['condition']==condition])
data=pd.concat(data_, axis=0)
columns=scores+['condition']
data=data[columns]
subplot_axis_id=axes["C"]
title_loc="left"
plot_title=r'$\bf{D3}$'
legend_loc="lower right"
ylabel='Centrality scores'
bbox_to_anchor=(1.3, 0.51)
_plot_bk_centrality_measures_(p_values_cell_type_squidpy_centralityScores, data, celltype, conditions, scores, subplot_axis_id, legend_loc=legend_loc, title_prefix=title_prefix, title_loc=title_loc, ylabel=ylabel, palette=palette, bbox_to_anchor=bbox_to_anchor, plot_title=plot_title)

# ========== Data for plot D: ==========
reference_celltype='Basal keratinocytes'
celltypes_=['Macrophages', 'Melanocytes', 'Smooth muscle cells', 'Suprabasal keratinocytes', 'Langerhans cells']
conditions=['AD', 'PSO', 'CTCL']
radii=[1, 5]
cnt=0
for celltype in celltypes_:
    cnt+=1
    title_prefix=None
    suptitle=None
    if cnt==1:
        # title_prefix='D '
        # suptitle='\t\t      '
        ylabel='Neighborhood\nenrichment\n(z-scores)'
        title_loc="left"
        if celltype=='Suprabasal keratinocytes':
            plot_title=r'\bf{D4}\n         BK vs. SK'
        else:
            # plot_title=str(f'BK vs.\n{celltype.lower()}')
            plot_title=str(r'$\bf{D4}$'+f'\n         BK vs. {abbr[celltype]}')
    else:
        # title_prefix=None
        # suptitle=None
        ylabel=None
        title_loc="center"
        if celltype=='Suprabasal keratinocytes':
            plot_title='BK vs. SK'
        else:
            # plot_title=str(f'BK vs.\n{celltype.lower()}')
            plot_title=str(f'BK vs. {abbr[celltype]}')
    data=pd.concat([squidpy_nhoodEnrichment_results[squidpy_nhoodEnrichment_results['celltype_1']==celltype][squidpy_nhoodEnrichment_results['celltype_2']=='Basal keratinocytes'], squidpy_nhoodEnrichment_results[squidpy_nhoodEnrichment_results['celltype_2']==celltype][squidpy_nhoodEnrichment_results['celltype_1']=='Basal keratinocytes']], axis=0)
    subplot_axis_id=axes[f"D{cnt}"]
    data_=[]
    for condition in conditions:
        data_.append(data[data['condition']==condition])
    data=pd.concat(data_, axis=0)
    _plot_bk_neighborhood_enrichment_(p_values_cell_type_squidpy_nhoodEnrichment, squidpy_nhoodEnrichment_results, data, celltype, reference_celltype, conditions, subplot_axis_id, title_prefix=title_prefix, plot_title=plot_title, suptitle=suptitle, ylabel=ylabel, palette=palette, title_loc=title_loc)
subfigs[2].subplots_adjust(wspace=0.45, hspace=1.3)

subplot_axis_id.text(-15.10, -9.00, textstr, fontsize=27, verticalalignment='top', bbox=props)

subfigs_nested = subfigs[0].subfigures(1, 2, wspace=0.00, width_ratios=[0.5, 0.5])
# ========== Create mosaic for: IIIA. Global heterogeneity distributions ==========
layout = [
    ["global_heterogeneity_scores"]
]
axes = subfigs_nested[0].subplot_mosaic(layout)
# subfigs_nested[0].set_facecolor('0.85')
subfigs_nested[0].suptitle(r'$\bf{A}$ Global heterogeneity scores', fontsize=30, x=0.25)
sample_results=pd.read_csv(os.path.join('../../../results', 'sample_results.csv'))
p_values_sample_wise=pd.read_csv(os.path.join('../../../results', 'p_values_global.csv'))
sample_results=sample_results.replace({"condition": conditions_abbreviations_dict})
p_values_sample_wise=p_values_sample_wise.replace({"condition_1": conditions_abbreviations_dict})
p_values_sample_wise=p_values_sample_wise.replace({"condition_2": conditions_abbreviations_dict})
title_prefix='G '
data=sample_results[['global_entropy', 'global_homophily', 'condition']]
data=data.melt('condition', var_name='global_heterogeneity_measure', value_name='global_heterogeneity_score')
args = dict(x='global_heterogeneity_measure', y='global_heterogeneity_score', data=data, hue="condition", hue_order=list(conditions))
subplot_axis_id=axes["global_heterogeneity_scores"]
_plot_global_score_distributions_(p_values_sample_wise, sample_results, args, conditions, subplot_axis_id, palette=palette)

# ========== Create mosaic for: IIIB. Visual analyses ==========
layout = [
    ["60", "60", "71", "71"],
    ["60", "60", "71", "71"],
    ["295", "295", "307", "307"],
    ["295", "295", "307", "307"]
]
axes = subfigs_nested[1].subplot_mosaic(layout)
# subfigs_nested[1].set_facecolor('0.85')
subfigs_nested[1].suptitle(r'$\bf{B}$ Visual analyses (T cells, basal keratinocytes)', fontsize=30, x=0.37)
subfigs_nested[1].subplots_adjust(hspace=1.0)

# =========== Read in spatial data: ===========
# Case study 2: Samples 60, 82, 83, 70, 71, 72, 295, 296, 302, 303, 304, 305, 306 and 307 are of interest.
# axes.set_ylabel(ylabel, fontsize=15, labelpad=5)
list_of_files=['60', '71', '295', '307']
axes_idx=[(0,0), (0,1), (1,0), (1,1), (2,0), (2,1)]
list_of_filenames=[]
for i in list_of_files:
    list_of_filenames.append(i+'.h5ad')
file_count=-1
for filename in list_of_filenames:
    file_count+=1
    title_prefix=None
    title_loc=None
    if file_count==0:
        # title_prefix='E '
        # title_loc='left'
        ylabel='AD\nsamples'
        legend=True
    # elif file_count==3:
    elif file_count==2:
        # title_prefix='F '
        # title_loc='left'
        ylabel='CTCL\nsamples'
        legend=None
    else:
        # title_prefix=None
        # title_loc=None
        ylabel=None
        legend=None
    print(f'Reading file {filename}...')
    subplot_axis_id=axes[list_of_files[file_count]]
    adata = sc.read_h5ad('../../../data/'+filename)
    prop_iodide=adata.uns['spatial']['images']['Propidium iodide']
    # isns.imgplot(img, cmap='gray', ax=subplot_axis_id)
    # subplot_axis_id.imshow(prop_iodide, cmap='gray', aspect='auto')
    # subplot_axis_id.margins(x=0)
    # subplot_axis_id.subplots_adjust(wspace=0, hspace=0)
    # subplot_axis_id.grid(False)
    # subplot_axis_id.axis('off')
    
    # Prepare spatial data:
    spatial_data=pd.DataFrame(adata.obsm['spatial'])
    spatial_data=spatial_data.rename(columns={0:'x' , 1:'y'})
    spatial_data['celltype']=list(adata.obs['celltype'])
    # # ^^^ Plot for all celltypes in same scatter plot. (Comment out when generating plots separately per {celltype} vs 'Others'): ^^^
    # # Plot segmented data:
    # title=f'Segmented tissue (all cell types)'
    # savefig_name=f'{adata.uns["Group"][0]}_{filename}(patient_id {adata.uns["patient_id"][0]}).pdf'
    # _generate_scatter_subplot_(prop_iodide, spatial_data, 'y', 'x', 'celltype', subplot_axis_id, title, savefig_name)
    # # ^^^
    
    # celltype_of_interest='T-cells'
    spatial_data_filtered_by_Tcells=spatial_data[spatial_data['celltype']=='T-cells']
    spatial_data_filtered_by_BKs=spatial_data[spatial_data['celltype']=='Basal keratinocytes']
    spatial_data_other_celltypes=spatial_data[spatial_data['celltype']!='T-cells'][spatial_data['celltype']!='Basal keratinocytes']
    spatial_data_other_celltypes['celltype']='Other'
    spatial_data_filtered_by_cell_type=pd.concat([spatial_data_other_celltypes, spatial_data_filtered_by_Tcells, spatial_data_filtered_by_BKs], axis=0)
    palette = {"T-cells":"red", "Basal keratinocytes":"blue", "Other":"#e6f8d1"}
    ncol=2
    _generate_scatter_subplot_(prop_iodide, spatial_data_filtered_by_cell_type, 'y', 'x', 'celltype', subplot_axis_id, palette, ncol, title_prefix=title_prefix, ylabel=ylabel, add_empty_line_before_title=True, title_loc=title_loc, legend=legend)



# ========== Generate, save and show final plot (fig1): ==========
# plt.subplots_adjust(wspace=0.02, hspace=0.5)
# plt.subplots_adjust(wspace=0.45, hspace=3.3)
plt.savefig('fig3_subplot_mosaic.pdf', format='pdf', bbox_inches='tight')
# fig.tight_layout(pad=100.0)
plt.show()

