import pandas as pd
import math


def _pick_column_maximizing_spread_(no_of_spread_scores_per_row, _celltypes_impGenes_, _celltypes_spreadScores_, good_split_genes, adata_df, geneName_upperThreshold_dict, clustering_tree, clustering_results, clustered_cells, hpa_data, errorStatement, suddenStopFlag):
    error_flag=1
    for column in range(no_of_spread_scores_per_row):
    
        celltypes_=[]
        celltypes_good=[]
        celltypes_bad=[]
        impGenes_=[]
        impGenes_good=[]
        impGenes_bad=[]
        spreadScores_=[]
        spreadScores_good=[]
        spreadScores_bad=[]
        
        for l in _celltypes_impGenes_:
            celltypes_.append(l)
            impGenes_.append(_celltypes_impGenes_[l][column])
            spreadScores_.append(_celltypes_spreadScores_[l][column])
        celltypes_counter=-1
        for i in celltypes_:
            celltypes_counter+=1
            if impGenes_[celltypes_counter] in good_split_genes:
                celltypes_good.append(i)
                impGenes_good.append(impGenes_[celltypes_counter])
                spreadScores_good.append(spreadScores_[celltypes_counter])
            else:
                celltypes_bad.append(i)
                impGenes_bad.append(impGenes_[celltypes_counter])
                spreadScores_bad.append(spreadScores_[celltypes_counter])
        
        if len(celltypes_good)>0:
            error_flag=0
            _pos_=pd.Series(spreadScores_good).idxmax()
            if math.isnan(_pos_):
                errorStatement="No_good_split_genes_found"
                suddenStopFlag=1
                break
            print("Maximum index position from the list: ", _pos_)
            anndata_clustered=adata_df[adata_df[impGenes_good[_pos_]]>geneName_upperThreshold_dict[impGenes_good[_pos_]][0]]
            anndata_clustered_indices=list(anndata_clustered.index)
            clustered_cells.append(anndata_clustered_indices)
            # ---
            clustering_tree[impGenes_good[_pos_]]=celltypes_good[_pos_]
            anndata_clustered_indices=[int(i) for i in anndata_clustered_indices]
            clustering_results[celltypes_good[_pos_]]=anndata_clustered_indices
            # ---
            adata_df=adata_df.drop(anndata_clustered_indices)
            hpa_data=hpa_data.drop(celltypes_good[_pos_])
            break
        # ---
    return error_flag, errorStatement, suddenStopFlag, clustered_cells, clustering_tree, clustering_results, adata_df, hpa_data
    
    