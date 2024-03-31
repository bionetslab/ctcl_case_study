# importing required library
import pandas as pd
import os
import seaborn as sns
sns.set_theme(style="whitegrid")
import matplotlib.pyplot as plt
from _plot_distributions_per_individual_celltypes_ import _plot_distributions_per_individual_celltypes_
from _plot_tcells_nhood_enrichment_ import _plot_tcells_nhood_enrichment_
from _plot_tcells_centrality_measures_ import _plot_tcells_centrality_measures_
from _utilities_ import _generate_scatter_subplot_
import scanpy as sc

# ========== Load and pre-process data required for generating plots: ==========
cell_results=pd.read_csv(os.path.join('../../results', 'cell_results.csv'))
p_values_cell_type=pd.read_csv(os.path.join('../../results', 'p_values_cell_type.csv'))
squidpy_nhoodEnrichment_results=pd.read_csv(os.path.join('../../results', 'squidpy_nhoodEnrichment_results.csv'))
p_values_cell_type_squidpy_nhoodEnrichment=pd.read_csv(os.path.join('../../results', 'p_values_cell_type_squidpy_nhoodEnrichment.csv'))
squidpy_centralityScores_results=pd.read_csv(os.path.join('../../results', 'squidpy_centralityScores_results.csv'))
p_values_cell_type_squidpy_centralityScores=pd.read_csv(os.path.join('../../results', 'p_values_cell_type_squidpy_centralityScores.csv'))
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

 
# creating grid for subplots
fig, axes = plt.subplots(6, 6, figsize=(10,10),gridspec_kw = {'wspace':0, 'hspace':0})
fig.set_figheight(20)
fig.set_figwidth(25)
 
ax_A = plt.subplot2grid(shape=(8, 6), loc=(0, 0), colspan=2, rowspan=2)
ax_B = plt.subplot2grid(shape=(8, 6), loc=(0, 2), colspan=2, rowspan=2)
ax_C = plt.subplot2grid(shape=(8, 6), loc=(0, 4), colspan=2, rowspan=2)
for count in range(1,7):
    exec_str=f'ax_D{count}=plt.subplot2grid(shape=(8, 6), loc=(2, {count-1}), colspan=1, rowspan=2)'
    print(exec_str)
    exec(exec_str)
list_of_files=['60', '83', '71']
count=0
for file in list_of_files:
    exec_str=f'ax_{file}=plt.subplot2grid(shape=(8, 6), loc=(4, {count}), colspan=2, rowspan=4)'
    print(exec_str)
    exec(exec_str)
    count+=2
 

# ========== Data for plot A: ==========
heterogeneity_measure='entropy'
celltype='T-cells'
conditions=['AD', 'CTCL']
radii=[1, 5]
subplot_axis_id=ax_A
legend_loc='lower right'
title_prefix='A. '
title_loc='left'
# Generate plot-1:
_plot_distributions_per_individual_celltypes_(p_values_cell_type, cell_results, celltype, heterogeneity_measure, conditions, radii, subplot_axis_id, legend_loc=legend_loc, title_prefix=title_prefix, title_loc=title_loc)


# ========== Data for plot B: ==========
heterogeneity_measure='egophily'
celltype='T-cells'
conditions=['AD', 'CTCL']
radii=[1, 5]
subplot_axis_id=ax_B
legend_loc='upper right'
title_prefix='B. '
title_loc='left'
# Generate plot-B:
_plot_distributions_per_individual_celltypes_(p_values_cell_type, cell_results, celltype, heterogeneity_measure, conditions, radii, subplot_axis_id, legend_loc=legend_loc, title_prefix=title_prefix, title_loc=title_loc)

# ========== Data for plot C: ==========
celltype='T-cells'
conditions=['AD', 'CTCL']
title_prefix='C. '
scores = ['degree_centrality', 'closeness_centrality']
data=squidpy_centralityScores_results[squidpy_centralityScores_results['celltypes']==celltype]
data_=[]
for condition in conditions:
    data_.append(data[data['condition']==condition])
data=pd.concat(data_, axis=0)
columns=scores+['condition']
data=data[columns]
subplot_axis_id=ax_C
title_loc="left"
legend_loc="lower right"
_plot_tcells_centrality_measures_(p_values_cell_type_squidpy_centralityScores, data, celltype, conditions, scores, subplot_axis_id, legend_loc=legend_loc, title_prefix=title_prefix, title_loc=title_loc)

# ========== Data for plot D: ==========
celltypes_=['Smooth muscle cells', 'Macrophages', 'Endothelial cells', 'B-cells', 'Fibroblasts', 'Basal keratinocytes']
conditions=['AD', 'CTCL']
radii=[1, 2, 3, 4, 5]
cnt=0
for celltype in celltypes_:
    cnt+=1
    if cnt==1:
        title_prefix='D. '
        suptitle='\t\t\t\t\t '
        ylabel='Neighborhood enrichment'
    else:
        title_prefix=None
        suptitle=None
        ylabel=None
    plot_title=str(f'{celltype}')
    data=pd.concat([squidpy_nhoodEnrichment_results[squidpy_nhoodEnrichment_results['celltype_1']==celltype][squidpy_nhoodEnrichment_results['celltype_2']=='T-cells'], squidpy_nhoodEnrichment_results[squidpy_nhoodEnrichment_results['celltype_2']==celltype][squidpy_nhoodEnrichment_results['celltype_1']=='T-cells']], axis=0)
    exec(f"subplot_axis_id=ax_D{cnt}")
    data_=[]
    for condition in conditions:
        data_.append(data[data['condition']==condition])
    data=pd.concat(data_, axis=0)
    _plot_tcells_nhood_enrichment_(p_values_cell_type_squidpy_nhoodEnrichment, squidpy_nhoodEnrichment_results, data, celltype, conditions, subplot_axis_id, title_prefix=title_prefix, plot_title=plot_title, suptitle=suptitle, ylabel=ylabel)



# # ========== Create mosaic for plot: ==========
# textstr_AD = 'Atopic dermatitis (AD)'
# # ax.hist(x, 50)
# # these are matplotlib.patch.Patch properties
# props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
# # place a text box in upper left in axes coords
# fig.text(0.23, 0.94, textstr_AD, fontsize=14,
#         verticalalignment='top', bbox=props)

# textstr_CTCL = 'Cutaneous T-cell lymphoma (CTCL)'
# # ax.hist(x, 50)
# # these are matplotlib.patch.Patch properties
# props = dict(boxstyle='round', facecolor='red', alpha=0.5)
# # place a text box in upper left in axes coords
# fig.text(0.55, 0.94, textstr_CTCL, fontsize=14,
#         verticalalignment='top', bbox=props)

# =========== Read in spatial data: ===========
# Case study 1: Samples 60, 82, 83, 70, 71, 72, 295, 296, 302, 303, 304, 305, 306 and 307 are of interest.
list_of_files=['60', '83', '71']
patient_ids=['Ecz\_01', 'Ecz\_12', 'Ecz\_06', 'TCL\_03', 'TCL\_07', 'TCL\_08']
axes_idx=[(0,0), (0,1), (1,0), (1,1), (2,0), (2,1)]
list_of_filenames=[]
for i in list_of_files:
    list_of_filenames.append(i+'.h5ad')
file_count=-1
for filename in list_of_filenames:
    file_count+=1
    print(f'Reading file {filename}...')
    exec(f'subplot_axis_id=ax_{list_of_files[file_count]}')
    adata = sc.read_h5ad('../../data/'+filename)
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
    
    celltype_of_interest='T-cells'
    spatial_data_filtered_by_cell_type=spatial_data[spatial_data['celltype']==celltype_of_interest]
    spatial_data_other_celltypes=spatial_data[spatial_data['celltype']!=celltype_of_interest]
    spatial_data_other_celltypes['celltype']='Other'
    spatial_data_filtered_by_cell_type=pd.concat([spatial_data_other_celltypes, spatial_data_filtered_by_cell_type], axis=0)
    palette = {f"{celltype_of_interest}":"red", "Other":"#e6f8d1"}
    ncol=2
    title_prefix=f'Patient\:id: {patient_ids[file_count]}'
    title=f'(sample {list_of_files[file_count]})'
    title=r'$\bf{' + title_prefix + '}$' + '\n' + title
    savefig_name=f'{adata.uns["Group"][0]}_{filename}(patient_id {adata.uns["patient_id"][0]})_{celltype_of_interest}.svg'
    _generate_scatter_subplot_(prop_iodide, spatial_data_filtered_by_cell_type, 'y', 'x', 'celltype', subplot_axis_id, palette, ncol)


# ========== Generate, save and show final plot (fig1): ==========
plt.subplots_adjust(wspace=0, hspace=0)
plt.savefig('fig1_gridspec.pdf', format='pdf', bbox_inches='tight')
# fig.tight_layout(pad=0.0)
plt.show()











