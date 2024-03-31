import pandas as pd
import os
import seaborn as sns
sns.set_theme(style="whitegrid")
import matplotlib.pyplot as plt
from _plot_distributions_per_individual_celltypes_ import _plot_distributions_per_individual_celltypes_
from _plot_tcells_nhood_enrichment_ import _plot_tcells_nhood_enrichment_
from _plot_tcells_centrality_measures_ import _plot_tcells_centrality_measures_
from _generate_scatterplot_ import _generate_scatterplot_
from _generate_histplot_ import _generate_histplot_
from _utilities_ import _generate_scatter_subplot_
import scanpy as sc
import seaborn_image as isns
import matplotlib as mpl
mpl.rcParams['savefig.pad_inches'] = 0
img = isns.load_image("polymer")

# =========== Read in spatial data: ===========
# Case study 1: Samples 60, 82, 83, 70, 71, 72, 295, 296, 302, 303, 304, 305, 306 and 307 are of interest.
list_of_files=['60', '83', '71', '295', '304', '307']
patient_ids=['Ecz\_01', 'Ecz\_12', 'Ecz\_06', 'TCL\_03', 'TCL\_07', 'TCL\_08']
axes_idx=[(0,0), (0,1), (1,0), (1,1), (2,0), (2,1)]
list_of_filenames=[]
for i in list_of_files:
    list_of_filenames.append(i+'.h5ad')
file_count=-1
for filename in list_of_filenames:
    file_count+=1
    print(f'Reading file {filename}...')
    adata = sc.read_h5ad('../../data/'+filename)
    fig, axes = plt.subplots(layout="constrained", figsize=(10,10))
    prop_iodide=adata.uns['spatial']['images']['Propidium iodide']
    # isns.imgplot(img, cmap='gray', ax=subplot_axis_id)
    ax=axes.imshow(prop_iodide, cmap='gray')
    # subplot_axis_id.margins(x=0)
    # subplot_axis_id.subplots_adjust(wspace=0, hspace=0)
    axes.grid(False)
    axes.axis('off')
    axes.margins(x=0)
    
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
    _generate_scatter_subplot_(prop_iodide, spatial_data_filtered_by_cell_type, 'y', 'x', 'celltype', axes, palette, ncol)
    # ========== Generate, save and show final plot (fig1): ==========
    plt.axis("off")
    plt.gca().set_axis_off()
    plt.subplots_adjust(wspace=0, hspace=0, bottom=0, top=0, left=0, right=0)
    plt.savefig(f'{list_of_files[file_count]}.jpg', format='jpg', bbox_inches='tight')
    # fig.tight_layout(pad=0.0)
    plt.show()














            
            
















