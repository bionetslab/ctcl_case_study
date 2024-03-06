import os
import pandas as pd
import scanpy as sc

def _prepare_spatial_data_for_celltype_assignment_(path_to_anndata_files):
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
            # _patientNumber_=list(anndata_allConditions.obsm['patient_number'])
            # patientNumber=list(np.unique(anndata_allConditions.obsm['patient_number']))
            # CellIndexNumber=list(anndata_allConditions.obsm['cell_index_number_within_patient'])
            anndata_allConditions_df.append(adata_df)
            anndata_allConditions_metadata.append(metadata_df)
    anndata_allConditions_df=pd.concat(anndata_allConditions_df, axis=0)
    anndata_allConditions_df=anndata_allConditions_df.reset_index(drop=True)
    anndata_allConditions_metadata=pd.concat(anndata_allConditions_metadata, axis=0)
    anndata_allConditions_metadata=anndata_allConditions_metadata.reset_index(drop=True)
    return anndata_allConditions_df, anndata_allConditions_metadata

def _save_celltype_assignment_results_to_csv_(clustering_results, anndata_allConditions_metadata):
    _celltype_=[]
    for i in clustering_results:
        l_=[i]*len(clustering_results[i])
        _celltype_.append(l_)
    _celltype_ = [x for xs in _celltype_ for x in xs]
    clustered_cells_df=anndata_allConditions_metadata
    clustered_cells_df['cell_type']=_celltype_
    return clustered_cells_df
    