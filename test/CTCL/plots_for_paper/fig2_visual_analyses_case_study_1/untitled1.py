import pandas as pd
import seaborn as sns
sns.set_theme(style="whitegrid")
import scanpy as sc
import matplotlib.pyplot as plt
from utilities import _generate_scatter_subplot_

if __name__ == '__main__':
    # ========== Create mosaic for plot: ==========
    # layout = [
    #     ["60", "295"],
    #     ["83", "304"],
    #     ["71", "307"]
    # ]
    layout = [
        ["295"]
    ]
    fig, axes = plt.subplot_mosaic(layout, figsize=(5,10))
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
    # list_of_files=['60', '83', '71', '295', '304', '307']
    # patient_ids=['Ecz\_01', 'Ecz\_12', 'Ecz\_06', 'TCL\_03', 'TCL\_07', 'TCL\_08']
    list_of_files=['295']
    patient_ids=['TCL\_03']
    axes_idx=[(0,0), (0,1), (1,0), (1,1), (2,0), (2,1)]
    list_of_filenames=[]
    for i in list_of_files:
        list_of_filenames.append(i+'.h5ad')
    file_count=-1
    for filename in list_of_filenames:
        file_count+=1
        print(f'Reading file {filename}...')
        subplot_axis_id=axes[list_of_files[file_count]]
        adata = sc.read_h5ad('../../data/'+filename)
        prop_iodide=adata.uns['spatial']['images']['Propidium iodide']
        subplot_axis_id.imshow(prop_iodide, cmap='gray')
        subplot_axis_id.grid(False)
        subplot_axis_id.axis('off')
        
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
        
        celltypes=spatial_data['celltype']
        spatial_data_filtered_by_cell_type=[]
        for celltype in celltypes:
            spatial_data_filtered_by_cell_type.append(spatial_data[spatial_data['celltype']==celltype])
        spatial_data_filtered_by_cell_type=pd.concat(spatial_data_filtered_by_cell_type, axis=0)
        # palette = {f"{celltype_of_interest}":"red", "Other":"#e6f8d1"}
        ncol=2
        title_prefix=f'Patient\:id: {patient_ids[file_count]}'
        title=f'(sample {list_of_files[file_count]})'
        title=r'$\bf{' + title_prefix + '}$' + '\n' + title
        title=None
        savefig_name=f'{adata.uns["Group"][0]}_{filename}(patient_id {adata.uns["patient_id"][0]}).jpg'
        _generate_scatter_subplot_(prop_iodide, spatial_data_filtered_by_cell_type, 'y', 'x', 'celltype', subplot_axis_id, title, savefig_name, ncol)
    
    # ========== Generate, save and show final plot (fig1): ==========
    plt.subplots_adjust(wspace=0)
    plt.savefig('295.jpg', format='jpg', bbox_inches='tight')
    plt.show()
















            
            

