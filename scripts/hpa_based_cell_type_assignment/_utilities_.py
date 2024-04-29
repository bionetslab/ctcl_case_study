import os
import pandas as pd
import scanpy as sc
import math

def _prepare_spatial_data_for_celltype_assignment_(path_to_anndata_files):
    sex_lookup={'m': 'Male', 'f': 'Female'}
    anndata_allConditions_df=[]
    anndata_allConditions_metadata=[]
    for filename in os.listdir('../data'):
        filename = os.fsdecode(filename)
        if filename.endswith('.h5ad'):
            sample_id=filename.split('.h5ad')[0]
            adata = sc.read_h5ad('../data/'+filename)
            # Convert anndata object to pandas DataFrame:
            adata_df=adata.to_df()
            metadata_df=pd.DataFrame()
            metadata_df['cell_idx']=list(adata_df.index)
            metadata_df['cell_id']=list(adata.obs['cell_id'])
            metadata_df['patient_id']=adata.uns['patient_id'][0]
            metadata_df['patient_label']=list(adata.obsm['patient_label'])
            metadata_df['sample_id']=sample_id
            metadata_df['condition']=adata.uns['Group'][0]
            metadata_df['age']=int(adata.uns['Age'][0])
            metadata_df['sex']=sex_lookup[adata.uns['Sex'][0]]
            metadata_df['x']=list(pd.DataFrame(adata.obsm['spatial'])[0])
            metadata_df['y']=list(pd.DataFrame(adata.obsm['spatial'])[1])
            anndata_allConditions_df.append(adata_df)
            anndata_allConditions_metadata.append(metadata_df)
    anndata_allConditions_df=pd.concat(anndata_allConditions_df, axis=0)
    anndata_allConditions_df=anndata_allConditions_df.reset_index(drop=True)
    anndata_allConditions_metadata=pd.concat(anndata_allConditions_metadata, axis=0)
    anndata_allConditions_metadata=anndata_allConditions_metadata.reset_index(drop=True)
    return anndata_allConditions_df, anndata_allConditions_metadata

def _save_celltype_assignment_results_to_csv_(clustering_results, anndata_allConditions_metadata, adata_df):
    df=[]
    df_metadata=[]
    
    for celltype, cell_indices_list in clustering_results.items():
        celltypes=[celltype]*len(cell_indices_list)
        df_=adata_df.loc[cell_indices_list]
        df_metadata_=anndata_allConditions_metadata.loc[cell_indices_list]
        df_['celltype']=celltypes
        df_metadata_['celltype']=celltypes
        df.append(df_)
        df_metadata.append(df_metadata_)
    df=pd.concat(df, axis=0)
    df.sort_index(inplace=True)
    df_metadata=pd.concat(df_metadata, axis=0)
    df_metadata.sort_index(inplace=True)
    return df, df_metadata

def _bin_column_in_dataframe_(df, bin_dict, old_col_name, new_col_name):
    df_=df.copy()
    df_=df_.sort_values(by=[old_col_name])
    old_column=list(df_[old_col_name])
    old_col_unique=sorted(list(set(df_.age)))
    age_bins={}
    for age in old_col_unique:
        bin=math.ceil((age-19)/10)
        age_bins[age]=bin_dict[bin]
    binned_column=[age_bins.get(item,item) for item in old_column]
    df_[new_col_name]=binned_column
    return df_















