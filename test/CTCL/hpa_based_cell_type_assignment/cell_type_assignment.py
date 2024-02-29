import pandas as pd
import math
from _plots_ import plot_heatmap
import os
import scanpy as sc
import sys
sys.path.insert(0, '../../')
from _fit_gaussian_mixture_model_ import _fit_gaussian_mixture_model_
from _calculate_spread_per_row_in_dataframe_ import _calculate_spread_per_row_in_dataframe_

if __name__ == '__main__':
    # Prepare spatial data for running the algorithm:
    anndata_allConditions_df=[]
    anndata_allConditions_metadata=[]
    for filename in os.listdir('../data'):
        filename = os.fsdecode(filename)
        if filename.endswith('.h5ad'):
            adata = sc.read_h5ad('../data/'+filename)
            # Convert anndata object to pandas DataFrame:
            adata_df=adata.to_df()
            metadata_df=pd.DataFrame()
            metadata_df['cell_idx']=list(adata_df.index)
            metadata_df['cell_id']=list(adata.obs['cell_id'])
            metadata_df['patient_id']=adata.uns['patient_id'][0]
            metadata_df['patient_label']=list(adata.obsm['patient_label'])
            # _patientNumber_=list(anndata_allConditions.obsm['patient_number'])
            # patientNumber=list(np.unique(anndata_allConditions.obsm['patient_number']))
            # CellIndexNumber=list(anndata_allConditions.obsm['cell_index_number_within_patient'])
            anndata_allConditions_df.append(adata_df)
    anndata_allConditions_df=pd.concat(anndata_allConditions_df, axis=0)
    anndata_allConditions_df=anndata_allConditions_df.reset_index(drop=True)
    anndata_allConditions_df_index=list(anndata_allConditions_df.index)    
    # Read in HPA data:
    expr_per_cell_type=pd.read_csv('cell_type_nTPM.csv')
    expr_per_cell_type_max_norm=pd.read_csv('cell_type_nTPM_max_norm.csv').rename(columns={'Unnamed: 0': "Celltypes"}).set_index('Celltypes')
    # Plot heatmap for CD4 and CD8 channels:
    plot_heatmap(expr_per_cell_type_max_norm, ['CD4','CD8A'], outfile='heatmap.jpg', title='Gene-wise max-normalized nTPM per cell type')
    list_of_channels_present_across_every_sample=['ADAM10', 'ITGAL', 'ITGAX', 'CD14', 'CD163', 'CD24', 'IL2RA', 'ITGB1',
    'CD3D', 'CD36', 'CD38', 'CD4', 'CD40', 'CD44', 'PTPRC', 'ICAM1', 'CD55',
    'NCAM1', 'CD6', 'SELP', 'CD63', 'CD68', 'CD69', 'CD8A', 'CD9', 'FAS',
    'KRT14', 'HLA-A', 'HLA-DRA', 'NOTCH1', 'NOTCH3', 'VIM']
    # Plot heatmap for all channels present across every sample:
    plot_heatmap(expr_per_cell_type_max_norm, list_of_channels_present_across_every_sample, outfile='heatmap.jpg', title='Gene-wise max-normalized nTPM per cell type')
    expr_per_cell_type_max_norm=expr_per_cell_type_max_norm[list_of_channels_present_across_every_sample]
    
    clustering_tree={}
    clustering_results={}
    clustered_cells=[]
    errorStatement="No_errors_all_celltypes_assigned"
    suddenStopFlag=0
    
    while len(list(expr_per_cell_type_max_norm.index))>0: # All columns have been eliminated from the HPA data.
        if len(anndata_allConditions_df)==0: # All cell types have been assigned
            errorStatement="No_more_cells_left_to_assign_in_database"
            break
        anndata_allConditions_colNames=list(anndata_allConditions_df.columns) # Retrieve all the channel names from the spatial data
        anndata_allConditions_colValues=[] # Retrieve all the column values from the spatial data
        for i in anndata_allConditions_colNames:
            anndata_allConditions_colValues.append(list(anndata_allConditions_df[i]))
        anndata_allConditions_dict=dict(zip(anndata_allConditions_colNames, anndata_allConditions_colValues))
        
        good_split_genes, geneName_upperThreshold_dict=_fit_gaussian_mixture_model_(anndata_allConditions_dict)
        
        # ===========================
        # Steps pending after fitting GMM model, and finding list of good split genes:
        # 1. Calculate spread per celltype per gene in set {good_split_genes}.
        # 2. Arrange all genes in set {good_split_genes} in descending order, per celltype.
        # 3. Pick {gene-g, celltype-C} that maximizes spread.
        # 4. Assign g+: C.
        # 5. Repeat step-1.
        # ===========================
        
        celltypes=list(expr_per_cell_type_max_norm.index) # List of celltypes present in HPA data.
        len_spread_scores_allColumns_for_index_desc_list, list_of_spread_scores_allColumns_per_row, column_values_desc_for_all_rows=_calculate_spread_per_row_in_dataframe_(expr_per_cell_type_max_norm)
        _celltypes_impGenes_=dict(zip(celltypes, list_of_spread_scores_allColumns_per_row))
        _celltypes_maxImpGeneValues_=dict(zip(celltypes, column_values_desc_for_all_rows))
        
        
        # Pick gene per celltype that maximizes spread:
        _pick_column_maximizing_spread()
        
        error_flag=1
        for k in range(len_spread_scores_allColumns_for_index_desc_list):
        
            _celltypes_=[]
            _impGenes_=[]
            _maxImpGeneValues_=[]
            
            for l in _celltypes_impGenes_:
                _celltypes_.append(l)
                _impGenes_.append(_celltypes_impGenes_[l][k])
                _maxImpGeneValues_.append(_celltypes_maxImpGeneValues_[l][k])
            
            
            celltypes_good=[]
            impGenes_good=[]
            maxImpGeneValues_good=[]
            
            celltypes_bad=[]
            impGenes_bad=[]
            maxImpGeneValues_bad=[]
            
            celltypes_counter=-1
            for i in _celltypes_:
                celltypes_counter+=1
                if _impGenes_[celltypes_counter] in good_split_genes:
                    celltypes_good.append(i)
                    impGenes_good.append(_impGenes_[celltypes_counter])
                    maxImpGeneValues_good.append(_maxImpGeneValues_[celltypes_counter])
                else:
                    celltypes_bad.append(i)
                    impGenes_bad.append(_impGenes_[celltypes_counter])
                    maxImpGeneValues_bad.append(_maxImpGeneValues_[celltypes_counter])
            
            
            # # ---
            # if len(celltypes_good)==0:
            #     errorStatement="No_good_split_genes_found"
            #     break
            # # ---
            
            if len(celltypes_good)>0:
                error_flag=0
                _pos_=pd.Series(maxImpGeneValues_good).idxmax()
                if math.isnan(_pos_):
                    errorStatement="No_good_split_genes_found"
                    suddenStopFlag=1
                    break
                print("Maximum Index position from the list: ", _pos_)
                anndata_clustered=anndata_allConditions_df[anndata_allConditions_df[impGenes_good[_pos_]]>geneName_upperThreshold_dict[impGenes_good[_pos_]][0]]
                anndata_clustered_indices=list(anndata_clustered.index)
                clustered_cells.append(anndata_clustered_indices)
                # ---
                clustering_tree[impGenes_good[_pos_]]=celltypes_good[_pos_]
                clustering_results[celltypes_good[_pos_]]=anndata_clustered_indices
                # ---
                anndata_allConditions_df=anndata_allConditions_df.drop(anndata_clustered_indices)
                expr_per_cell_type_max_norm=expr_per_cell_type_max_norm.drop(celltypes_good[_pos_])
                break
            # ---
        if suddenStopFlag==1:
            break
        if error_flag==1:
            errorStatement="No_good_split_genes_found"
            break
    
    clustered_cells = [
        x
        for xs in clustered_cells
        for x in xs
    ]
    
    non_clustered_cells=list(set(anndata_allConditions_df_index).difference(set(clustered_cells)))
    
    clustering_results['Unknown']=non_clustered_cells
    
    _celltype_=[]
    _patientIndex_=[]
    for i in clustering_results:
        l_=[i]*len(clustering_results[i])
        _celltype_.append(l_)
        _patientIndex_.append(clustering_results[i])
    
    _celltype_ = [x for xs in _celltype_ for x in xs]
    
    _patientIndex_ = [x for xs in _patientIndex_ for x in xs]
    
    Clustered_Cells_DF=pd.DataFrame()
    Clustered_Cells_DF.index=anndata_allConditions_df_index
    Clustered_Cells_DF['patient_number']=_patientNumber_
    Clustered_Cells_DF['cell_index_number_within_patient']=CellIndexNumber
    Clustered_Cells_DF=Clustered_Cells_DF.drop(non_clustered_cells, axis=0)
    
    Index_Celltype_DF=pd.DataFrame()
    Index_Celltype_DF.index=_patientIndex_
    Index_Celltype_DF['celltype']=_celltype_
    
    Clustered_Cells_DF=Clustered_Cells_DF.join(Index_Celltype_DF)
    
    # ----------
    
    