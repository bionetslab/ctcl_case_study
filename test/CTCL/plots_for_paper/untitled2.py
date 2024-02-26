import pandas as pd
import scipy.stats as stats
import os
import itertools as itt
from IPython.display import set_matplotlib_formats
import seaborn as sns
sns.set_theme(style="whitegrid")
import matplotlib.pyplot as plt
import numpy as np
import itertools
from statannotations.Annotator import Annotator
from decimal import Decimal
from _plot_distributions_per_individual_celltypes_ import _plot_distributions_per_individual_celltypes_
from _plot_tcells_nhood_enrichment_ import _plot_tcells_nhood_enrichment_
import scanpy as sc
import random

def merge(list1, list2):
    merged_list = [(list1[i], list2[i]) for i in range(0, len(list1))]
    return merged_list

if __name__ == '__main__':
    cell_results=pd.read_csv('../results/cell_results.csv')
    cell_results=cell_results[cell_results['cell_type']=='T-cells']
    
    # ========== Prepare cell_results data by sample ids: ==========
    sample_ids=list(set(cell_results.sample_id))
    mean_local_entropy_3=[]
    mean_egophily_3=[]
    condition=[]
    df_by_sample_ids=pd.DataFrame()
    df_by_sample_ids['sample_ids']=sample_ids
    for sampleId in sample_ids:
        sample_res=cell_results[cell_results['sample_id']==sampleId]
        mean_local_entropy_3.append(np.mean(list(sample_res['local_entropy_3'])))
        mean_egophily_3.append(np.mean(list(sample_res['egophily_3'])))
        condition.append(list(sample_res.condition)[0])
    df_by_sample_ids['mean_local_entropy_3']=mean_local_entropy_3
    df_by_sample_ids['mean_egophily_3']=mean_egophily_3
    df_by_sample_ids['condition']=condition
    
    # ========== Prepare cell_results data by patient ids: ==========
    patient_ids=list(set(cell_results.patient_id))
    mean_local_entropy_3=[]
    mean_egophily_3=[]
    condition=[]
    df_by_patient_ids=pd.DataFrame()
    df_by_patient_ids['patient_ids']=patient_ids
    patientIds_sampleIds_lookup={}
    for patientId in patient_ids:
        patient_res=cell_results[cell_results['patient_id']==patientId]
        patientIds_sampleIds_lookup[patientId]=list(set(patient_res['sample_id']))
        mean_local_entropy_3.append(np.mean(list(patient_res['local_entropy_3'])))
        mean_egophily_3.append(np.mean(list(patient_res['egophily_3'])))
        condition.append(list(patient_res.condition)[0])
    df_by_patient_ids['mean_local_entropy_3']=mean_local_entropy_3
    df_by_patient_ids['mean_egophily_3']=mean_egophily_3
    df_by_patient_ids['condition']=condition
    
    # =========== Read in spatial data: ===========
    adatas={}
    for filename in os.listdir('../data'):
        filename = os.fsdecode(filename)
        print(f'Reading file {filename}...')
        if filename.endswith('.h5ad'):
            adata = sc.read_h5ad('../data/'+filename)
            adatas[filename]=adata
            prop_iodide=adata.uns['spatial']['images']['Propidium iodide']
            segmentation=adata.uns['spatial']['segmentation']
            x_coords=np.nonzero(adata.uns['spatial']['segmentation'])[0]
            y_coords=np.nonzero(adata.uns['spatial']['segmentation'])[1]
            coords=merge(x_coords, y_coords)
            res = list(map(segmentation.__getitem__, coords))
            fig, axes = plt.subplots(figsize=(10,10))
            axes.imshow(prop_iodide, cmap='gray')
            mapping=dict(zip(list(adata.obs['cell_id']), list(adata.obs['celltype'])))
            adata_df=pd.DataFrame()
            adata_df['x']=x_coords
            adata_df['y']=y_coords
            adata_df['celltype']=[mapping.get(item,item) for item in res]
            
            # ========== Subsample n cells per sample (i.e., section of tissue) for generation of scatter plot: ==========
            delta_sample=random.sample(coords, 10000)
            res = list(map(segmentation.__getitem__, delta_sample))
            mapping=dict(zip(list(adata.obs['cell_id']), list(adata.obs['celltype'])))
            df = pd.DataFrame(delta_sample, columns =['x', 'y'])
            df['celltype']=[mapping.get(item,item) for item in res]
            
            
            
            xlabel='x'
            ylabel='y'
            # ax = sns.scatterplot(data=adata_df, x="x", y="y", hue="celltype", palette='Set3', s=30)
            ax = sns.scatterplot(data=df, x="y", y="x", hue="celltype", palette='Set3', s=30)
            ax.set(xlabel=xlabel, ylabel=ylabel, title=f'Cell types')
            ax.invert_yaxis()
            ax.tick_params(left=False, bottom=False)
            sns.move_legend(
                ax, "upper right",
                bbox_to_anchor=(1.03, -0.08), ncol=3, title=None, frameon=False,
            )
            
            
            plt.savefig(f'{adata.uns["Group"][0]}_{filename}.pdf', format='pdf', bbox_inches='tight')
            plt.show()
            plt.close()
            
            # # ===== Scatter plot with overlay =======
            
            # local_heterogeneity_measure='local_entropy_2'
            # # =====
            # from IPython.display import set_matplotlib_formats
            # import seaborn as sns
            # sns.set_theme(style="whitegrid", palette="muted")
            # # ---
            # sns.set_style("whitegrid", {'axes.grid' : False})
            # # fig, (axes1, axes2) = plt.subplots(1,2)
            # layout = [["A"], ["B"]]
            # fig, axes = plt.subplot_mosaic(layout, figsize=(10,10))
            # # # plt.suptitle(f"\nPatient#{i} (condition: {_label_})", x=0.40, y=1.00)
            # # --> axes["A"].imshow(tissue_image_gray, cmap='gray')
            # # ---
            # xlabel='x'
            # ylabel='y'
            # # --> ax = sns.scatterplot(ax=axes["A"], data=DF_seqfish, x="x", y="y", hue=local_heterogeneity_measure, palette='RdBu', s=30)
            # # --> ax.set(xlabel=xlabel, ylabel=ylabel, title=f'Local heterogeneity score')
            # # ax.invert_yaxis()
            # axes["A"].tick_params(left=False, bottom=False)
            # # ---
            # # norm = plt.Normalize(df_coords['local_heterogeneity_score'].min(), df_coords['local_heterogeneity_score'].max())
            # norm = plt.Normalize(0, 1)
            # sm = plt.cm.ScalarMappable(cmap="RdBu", norm=norm)
            # sm.set_array([])
            # # Remove the legend and add a colorbar
            # # --> ax.get_legend().remove()
            # # ax.figure.colorbar(sm, location = 'right')
            # cax_=fig.add_axes([0.1, 0.6, 0.03, 0.3])
            # # --> ax.figure.colorbar(sm, cax=cax_)
            # # ========== X ==========
            # # ---
            # sns.set_style("whitegrid", {'axes.grid' : False})
            # # fig, axes = plt.subplots()
            # axes["B"].imshow(prop_iodide, cmap='gray')
            # # ---
            # xlabel='x'
            # ylabel='y'
            # ax = sns.scatterplot(ax=axes["B"], data=DF_seqfish, x="x", y="y", hue="Cluster", palette='Set3', s=30)
            # ax.set(xlabel=xlabel, ylabel=ylabel, title=f'Cell types')
            # # ax.invert_yaxis()
            # axes["B"].tick_params(left=False, bottom=False)
            # # ---
            # sns.move_legend(axes["B"], "upper left", bbox_to_anchor=(-1.0, .8), title='Cell types')
            #         # ---
            #         # # ---
            #         # # norm = plt.Normalize(df_coords['local_heterogeneity_score'].min(), df_coords['local_heterogeneity_score'].max())
            #         # norm = plt.Normalize(0, 1)
            #         # sm = plt.cm.ScalarMappable(cmap="RdBu", norm=norm)
            #         # sm.set_array([])
            #         # # Remove the legend and add a colorbar
            #         # ax.get_legend().remove()
            #         # ax.figure.colorbar(sm)
            #         # # ---
            # plt.xticks([])
            # plt.yticks([])
            # fig.tight_layout()
            # # # plt.savefig(f'{_label_}\Patient#{i} (condition: {_label_}).pdf', format='pdf', bbox_inches='tight')
            # plt.show()
            # plt.close()
            # # =====
            
            
            
            
            
            
            
            
            
            
            
            
            

