from _utilities_ import _save_celltype_assignment_results_, _generate_metadata_df_, declare_celltype_annotated_df_from_anndata, _pick_column_maximizing_spread_, _calculate_spread_per_row_in_dataframe_, _fit_gaussian_mixture_model_, read_hpa_data, preprocess_anndata_dataframe, _save_anndata_as_h5ad_
import scanpy as sc
import hdf5plugin

hpa_data_path='data/cell_type_nTPM_max_norm.csv'
channels=['ITGAL', 'ITGAX', 'CD14', 'CD163', 'LY75', 'MRC1', 'CD24', 'IL2RA', 'ITGB1', 'CD3D', 
    'CD36', 'CD38', 'CD4', 'CD40', 'CD44', 'PTPRC', 'CD52', 'ICAM1', 'CD55', 'NCAM1', 'CD6', 'SELP', 
    'CD63', 'CD68', 'CD69', 'CD8A', 'CD9', 'FAS', 'KRT14', 'HLA-A', 'HLA-DRA', 'NOTCH1', 'NOTCH3', 'PPARG', 'CTNNB1']
expr_per_cell_type_max_norm=read_hpa_data(hpa_data_path, channels)
clustering_tree={}
clustering_results={}
anndata_allConditions=sc.read_h5ad('anndata_allConditions.h5ad')
anndata_allConditions_df=anndata_allConditions.to_df()
metadata_df=_generate_metadata_df_(anndata_allConditions_df, anndata_allConditions)
celltype_annotated_df, _sample_id_, Labels=declare_celltype_annotated_df_from_anndata(anndata_allConditions_df, anndata_allConditions)
anndata_allConditions_df_index, _sample_id_, CellIndexNumber=preprocess_anndata_dataframe(anndata_allConditions_df, anndata_allConditions)


# =================== Celltype assignment algorithm: ===================
# Step-1: Fit GMM model.
# Step-2: Find list of good split genes.
# Step-3: Calculate spread per celltype per gene in set {good_split_genes}.
# Step-4: Arrange all genes in set {good_split_genes} in descending order, per celltype.
# Step-5: Pick {gene-g, celltype-C} that maximizes spread.
# Step-6: Assign g+: C.
# ** Step-7: Repeat step-1. ** (Already taken care of by the outer while loop.)
# Step-8: Assign unassigned cells to 'Unknown' type.
# Step-9: Map back assigned celltypes to the original cells.
# Step-10: Save celltype assignment results to csv file.
# ======================================================================
clustering_tree={}
clustering_results={}
all_columns_list=list(anndata_allConditions_df.columns)
anndata_allConditions_df_index=list(anndata_allConditions_df.index)
clustered_cells=[]
errorStatement="No_errors_all_celltypes_assigned"
suddenStopFlag=0
while len(list(expr_per_cell_type_max_norm.index))>0:
    if len(anndata_allConditions_df)==0:
        errorStatement="No_more_cells_left_to_assign_in_database"
        break
    anndata_allConditions_colNames=list(anndata_allConditions_df.columns)
    anndata_allConditions_colValues=[]
    for i in anndata_allConditions_colNames:
        anndata_allConditions_colValues.append(list(anndata_allConditions_df[i]))
    anndata_allConditions_dict=dict(zip(anndata_allConditions_colNames, anndata_allConditions_colValues))
    # Steps 1 (Fit GMM model) and 2 (Find list of good split genes):
    good_split_genes, geneName_upperThreshold_dict=_fit_gaussian_mixture_model_(anndata_allConditions_dict)
    celltypes=list(expr_per_cell_type_max_norm.index)
    # Step-3 (Calculate spread per celltype per gene in set {good_split_genes}) and 4 (Arrange all genes in set {good_split_genes} in descending order, per celltype):
    no_of_spread_scores_per_row, list_of_genes, spread_scores=_calculate_spread_per_row_in_dataframe_(expr_per_cell_type_max_norm)
    _celltypes_impGenes_=dict(zip(celltypes, list_of_genes))
    _celltypes_spreadScores_=dict(zip(celltypes, spread_scores))
    # Step-5 (Pick {gene-g, celltype-C} that maximizes spread) and 6 (Assign g+: C):
    error_flag, errorStatement, suddenStopFlag, clustered_cells, clustering_tree, clustering_results, anndata_allConditions_df, expr_per_cell_type_max_norm=_pick_column_maximizing_spread_(no_of_spread_scores_per_row, _celltypes_impGenes_, _celltypes_spreadScores_, good_split_genes, anndata_allConditions_df, geneName_upperThreshold_dict, clustering_tree, clustering_results, clustered_cells, expr_per_cell_type_max_norm, errorStatement, suddenStopFlag)
    if suddenStopFlag==1:
        break
    if error_flag==1:
        errorStatement="No_good_split_genes_found"
        break
clustered_cells = [x for xs in clustered_cells for x in xs]
non_clustered_cells=list(set(anndata_allConditions_df_index).difference(set(clustered_cells)))
# Step-8 (Assign unassigned cells to 'Unknown' type):
clustering_results['Unknown']=non_clustered_cells
# Step-9 (Map back assigned celltypes to the original cells):
clustered_cells_df, clustered_cells_metadata, clustered_cells_combined=_save_celltype_assignment_results_(clustering_results, metadata_df, celltype_annotated_df)
# Save celltype assignment results:
clustered_cells_combined.to_csv('result/sample-wise celltypes (HPA-based clustering).csv')
# Save df as anndata object containing all samples:
adata=sc.AnnData(clustered_cells_combined[list(set(channels).intersection(set(clustered_cells_combined.columns)))])
# Save anndata object as .h5ad file:
filename='result/celltype_assigned_anndata.h5ad'
list_of_obsms=['x', 'y', 'sample_id', 'condition', 'celltype', 'cell_idx', 'cell_id', '_index_']
adata=_save_anndata_as_h5ad_(adata, filename, clustered_cells_combined, list_of_obsms)
adata.write_h5ad(filename, compression=hdf5plugin.FILTERS["zstd"], compression_opts=hdf5plugin.Zstd(clevel=5).filter_options)

# # ============= Generate separate .h5ad files: =============
# sample_ids=list(np.unique(clustered_cells_combined.sample_id))
# for i in sample_ids:
#     filename=f'result/{str(int(i))}.h5ad'
#     clustered_cells_combined_filtered=clustered_cells_combined[clustered_cells_combined['sample_id']==i]
#     adata_filtered=sc.AnnData(clustered_cells_combined_filtered[list(set(channels).intersection(set(clustered_cells_combined.columns)))])
#     list_of_obsms=['x', 'y', 'sample_id', 'condition', 'celltype', 'cell_idx', 'cell_id', '_index_']
#     adata_filtered=_save_anndata_as_h5ad_(adata_filtered, filename, clustered_cells_combined_filtered, list_of_obsms)
#     adata_filtered.write_h5ad(filename, compression=hdf5plugin.FILTERS["zstd"], compression_opts=hdf5plugin.Zstd(clevel=5).filter_options)
# # ===========================================================
