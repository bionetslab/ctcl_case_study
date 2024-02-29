import pandas as pd
import seaborn as sns
sns.set_theme(style="whitegrid")
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from utilities import _generate_scatter_subplot_in_mosaic_

if __name__ == '__main__':
    # ========== Create mosaic for plot: ==========
    # layout = [
    #     ["60", "295"],
    #     ["83", "304"],
    #     ["71", "307"]
    # ]
    # fig, axes = plt.subplot_mosaic(layout, figsize=(10,10))
    # fig.subplots_adjust(hspace=0.125, wspace=0.125)
    
    # fig = plt.figure(figsize=(20,20))
    # gs = gridspec.GridSpec(3, 2, figure=fig)
    # gs.update(wspace=0.0, hspace=0.05)
    
    fig, axes = plt.subplots(3,2, figsize=(5,8))
    

    # =========== Read in spatial data: ===========
    # Case study 1: Samples 60, 82, 83, 70, 71, 72, 295, 296, 302, 303, 304, 305, 306 and 307 are of interest.
    list_of_files=['60', '83', '71', '295', '304', '307']
    axes_idx=[(0,0), (0,1), (1,0), (1,1), (2,0), (2,1)]
    list_of_filenames=[]
    for i in list_of_files:
        list_of_filenames.append(i+'.h5ad')
    file_count=-1
    for filename in list_of_filenames:
        file_count+=1
        print(f'Reading file {filename}...')
        subplot_axis_id=axes[axes_idx[file_count][0], axes_idx[file_count][1]]
        adata = sc.read_h5ad('../../data/'+filename)
        prop_iodide=adata.uns['spatial']['images']['Propidium iodide']
        subplot_axis_id.imshow(prop_iodide, cmap='gray')
        subplot_axis_id.grid(False)
        subplot_axis_id.axis('off')
        
        # # Prepare spatial data:
        # spatial_data=pd.DataFrame(adata.obsm['spatial'])
        # spatial_data=spatial_data.rename(columns={0:'x' , 1:'y'})
        # spatial_data['celltype']=list(adata.obs['celltype'])
        # # # ^^^ Plot for all celltypes in same scatter plot. (Comment out when generating plots separately per {celltype} vs 'Others'): ^^^
        # # # Plot segmented data:
        # # title=f'Segmented tissue (all cell types)'
        # # savefig_name=f'{adata.uns["Group"][0]}_{filename}(patient_id {adata.uns["patient_id"][0]}).pdf'
        # # _generate_scatter_subplot_in_mosaic_(prop_iodide, spatial_data, 'y', 'x', 'celltype', subplot_axis_id, title, savefig_name)
        # # # ^^^
        
        # celltype_of_interest='T-cells'
        # spatial_data_filtered_by_cell_type=spatial_data[spatial_data['celltype']==celltype_of_interest]
        # spatial_data_other_celltypes=spatial_data[spatial_data['celltype']!=celltype_of_interest]
        # spatial_data_other_celltypes['celltype']='Other'
        # spatial_data_filtered_by_cell_type=pd.concat([spatial_data_filtered_by_cell_type, spatial_data_other_celltypes], axis=0)
        # palette = {f"{celltype_of_interest}":"red", "Other":"#e6f8d1"}
        # ncol=2
        # title=f'Segmented tissue ({celltype_of_interest})'
        # savefig_name=f'{adata.uns["Group"][0]}_{filename}(patient_id {adata.uns["patient_id"][0]})_{celltype_of_interest}.svg'
        # _generate_scatter_subplot_in_mosaic_(prop_iodide, spatial_data_filtered_by_cell_type, 'y', 'x', 'celltype', subplot_axis_id, title, savefig_name, palette, ncol)
    
    # ========== Generate, save and show final plot (fig1): ==========
    plt.savefig('fig2_caseStudy1_Tcells_ADvsCTCL.pdf', format='pdf', bbox_inches='tight')
    plt.show()
            
            
            

