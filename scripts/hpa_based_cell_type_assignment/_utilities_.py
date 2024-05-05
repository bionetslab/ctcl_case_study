import pandas as pd
import math
import numpy as np
import statistics
from sklearn.mixture import GaussianMixture

def _save_celltype_assignment_results_(clustering_results, anndata_allConditions_metadata, adata_df):
    df=[]
    df_metadata=[]
    for celltype, cell_indices_list in clustering_results.items():
        celltypes=[celltype]*len(cell_indices_list)
        cell_indices_list_str = [str(x) for x in cell_indices_list]
        # cell_indices_list_int = [int(x) for x in cell_indices_list]
        df_=adata_df.loc[cell_indices_list_str]
        df_metadata_=anndata_allConditions_metadata.loc[cell_indices_list_str]
        df_['celltype']=celltypes
        df_metadata_['celltype']=celltypes
        df.append(df_)
        df_metadata.append(df_metadata_)
    df=pd.concat(df, axis=0)
    df_metadata=pd.concat(df_metadata, axis=0)
    df_combined=pd.concat([df, df_metadata], axis=1)
    return df, df_metadata, df_combined

def _generate_metadata_df_(anndata_allConditions_df, anndata_allConditions):
    metadata_df=pd.DataFrame()
    metadata_df['cell_idx']=list(anndata_allConditions_df.index)
    metadata_df['cell_id']=list(anndata_allConditions.obsm['cell_index_number_within_patient'])
    metadata_df['_index_']=list(anndata_allConditions.obsm['_index_'])
    # metadata_df['patient_label']=list(anndata_allConditions.obsm['patient_number'])
    metadata_df['condition']=list(anndata_allConditions.obsm['label'])
    metadata_df['x']=list(anndata_allConditions.obsm['x'])
    metadata_df['y']=list(anndata_allConditions.obsm['y'])
    metadata_df.index=[str(x) for x in list(metadata_df.index)]
    return metadata_df

def declare_celltype_annotated_df_from_anndata(anndata_df, anndata):
    df=anndata_df.copy()
    # cell_coords_df=pd.DataFrame([l.tolist() for l in anndata.obsm['_cell_coords_'].tolist()])
    # cell_coords_df=cell_coords_df.rename(columns={0:'x', 1:'y'})
    _sample_id_=list(anndata.obsm['patient_number'])
    Labels=list(anndata.obsm['label'])
    df['x']=list(anndata.obsm['x'])
    df['y']=list(anndata.obsm['y'])
    df['sample_id']=_sample_id_
    df['condition']=list(Labels)
    return df, _sample_id_, Labels

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
            anndata_clustered_indices=[str(x) for x in anndata_clustered_indices]
            adata_df=adata_df.drop(anndata_clustered_indices)
            hpa_data=hpa_data.drop(celltypes_good[_pos_])
            break
        # ---
    return error_flag, errorStatement, suddenStopFlag, clustered_cells, clustering_tree, clustering_results, adata_df, hpa_data   

def _calculate_spread_per_row_in_dataframe_(dataframe):
    list_of_genes=[]
    spread_scores=[]
    for index in dataframe.index:
        row=dataframe.loc[[index]]
        update_df = dataframe.drop(index)
        row_complement=pd.DataFrame(update_df.max()).T.rename(index={0: index})
        spread_scores_allColumns_for_index=row.sub(row_complement) # Spread score (x-max(-x)) for index=index, across all columns present in HPA data.
        index_with_highest_spread=spread_scores_allColumns_for_index.index[0]
        spread_scores_allColumns_for_index_desc_dict=spread_scores_allColumns_for_index[spread_scores_allColumns_for_index.iloc[-1,:].sort_values(ascending=False).index].T[index_with_highest_spread].to_dict()
        spread_scores_allColumns_for_index_desc_list=list(spread_scores_allColumns_for_index_desc_dict.keys())
        no_of_spread_scores_per_row=len(spread_scores_allColumns_for_index_desc_list)
        list_of_genes.append(spread_scores_allColumns_for_index_desc_list)
        spread_scores.append(list(spread_scores_allColumns_for_index_desc_dict.values()))
    return no_of_spread_scores_per_row, list_of_genes, spread_scores

def _fit_gaussian_mixture_model_(anndata_allConditions_dict):
    anndata_allConditions_colNames=list(anndata_allConditions_dict.keys())
    anndata_allConditions_gmm=[]
    anndata_allConditions_gmm_labels=[]
    anndata_allConditions_gmm_upper_threshold=[]
    anndata_allConditions_gmm_lower_threshold=[]
    # -----
    good_split_genes=[]
    bad_split_genes=[]
    # -----
    for i in anndata_allConditions_dict:
        X_list=list(anndata_allConditions_dict[i])
        X=np.array(X_list)
        X=X.reshape(-1, 1)
        gm = GaussianMixture(n_components=2, random_state=0).fit(X)
        anndata_allConditions_gmm.append(gm)
        labels = gm.predict(X)
        anndata_allConditions_gmm_labels.append(labels)
        X_0=[]
        X_1=[]
        cnt=-1
        for j in labels:
            cnt+=1
            if j==0:
                X_0.append(X_list[cnt])
                X_1.append(X_list[cnt])
        X=[X_0, X_1]
        means_=gm.means_
        stddevs_0= statistics.pstdev(X_0)
        stddevs_1= statistics.pstdev(X_1)
        stddevs_=[stddevs_0, stddevs_1]
        # -----
        upper_threshold=means_[1]-1.96*stddevs_[1]
        anndata_allConditions_gmm_upper_threshold.append(upper_threshold)
        lower_threshold=means_[0]+1.96*stddevs_[0]
        anndata_allConditions_gmm_lower_threshold.append(lower_threshold)
        # -----
        if upper_threshold>lower_threshold:
            good_split_genes.append(i)
        else:
            bad_split_genes.append(i)
    geneName_upperThreshold_dict=dict(zip(anndata_allConditions_colNames, anndata_allConditions_gmm_upper_threshold))
    return good_split_genes, geneName_upperThreshold_dict

def read_hpa_data(hpa_data_path, channels):
    expr_per_cell_type_max_norm=pd.read_csv(hpa_data_path)
    expr_per_cell_type_max_norm=expr_per_cell_type_max_norm.rename(columns={'Unnamed: 0': "Celltypes"})
    expr_per_cell_type_max_norm=expr_per_cell_type_max_norm.set_index('Celltypes')
    expr_per_cell_type_max_norm=expr_per_cell_type_max_norm[channels]
    return expr_per_cell_type_max_norm

def preprocess_anndata_dataframe(anndata_allConditions_df, anndata_allConditions):
    anndata_allConditions_df_index=list(anndata_allConditions_df.index)
    _sample_id_=list(anndata_allConditions.obsm['patient_number'])
    CellIndexNumber=list(anndata_allConditions.obsm['cell_index_number_within_patient'])
    return anndata_allConditions_df_index, _sample_id_, CellIndexNumber

def _save_anndata_as_h5ad_(anndata, filename, df, list_of_obsms):
    for col in list_of_obsms:
        anndata.obsm[col]=np.array(df[col])
    return anndata
