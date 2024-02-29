import pandas as pd
import seaborn as sns
sns.set_theme(style="whitegrid")
import scanpy as sc
from utilities import _generate_scatter_plot_with_overlay_

if __name__ == '__main__':
    # =========== Read in spatial data: ===========
    # Case study 1: Samples 60, 82, 83, 70, 71, 72, 295, 296, 302, 303, 304, 305, 306 and 307 are of interest.
    list_of_filenames=['60', '82', '83', '70', '71', '72', '295', '296', '302', '303', '304', '305', '306', '307']
    for i in range(len(list_of_filenames)):
        list_of_filenames[i]+='.h5ad'
    
    for filename in list_of_filenames:
        print(f'Reading file {filename}...')
        adata = sc.read_h5ad('../../data/'+filename)
        prop_iodide=adata.uns['spatial']['images']['Propidium iodide']
        # Prepare spatial data:
        spatial_data=pd.DataFrame(adata.obsm['spatial'])
        spatial_data=spatial_data.rename(columns={0:'x' , 1:'y'})
        spatial_data['celltype']=list(adata.obs['celltype'])
        # # ^^^ Plot for all celltypes in same scatter plot. (Comment out when generating plots separately per {celltype} vs 'Others'): ^^^
        # # Plot segmented data:
        # title=f'Segmented tissue (all cell types)'
        # savefig_name=f'{adata.uns["Group"][0]}_{filename}(patient_id {adata.uns["patient_id"][0]}).pdf'
        # _generate_scatter_plot_with_overlay_(prop_iodide, spatial_data, 'y', 'x', 'celltype', title, savefig_name)
        # # ^^^
        
        # --- Next we plot celltype-wise: ---
        # Case study 1: 'T-cells' are of interest.
        # for all celltypes:
        # celltypes=list(set(adata.obs['celltype']))
        # for particular celltypes:
        celltypes=['T-cells']
        for celltype in celltypes:
            spatial_data_filtered_by_cell_type=spatial_data[spatial_data['celltype']==celltype]
            spatial_data_other_celltypes=spatial_data[spatial_data['celltype']!=celltype]
            spatial_data_other_celltypes['celltype']='Other'
            spatial_data_filtered_by_cell_type=pd.concat([spatial_data_filtered_by_cell_type, spatial_data_other_celltypes], axis=0)
            palette = {f"{celltype}":"red", "Other":"#e6f8d1"}
            ncol=2
            title=f'Segmented tissue ({celltype})'
            savefig_name=f'{adata.uns["Group"][0]}_{filename}(patient_id {adata.uns["patient_id"][0]})_{celltype}.svg'
            _generate_scatter_plot_with_overlay_(prop_iodide, spatial_data_filtered_by_cell_type, 'y', 'x', 'celltype', title, savefig_name, palette, ncol)
            
            
            

