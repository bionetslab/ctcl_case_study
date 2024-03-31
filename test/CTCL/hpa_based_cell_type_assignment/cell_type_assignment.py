import pandas as pd
from _plots_ import plot_heatmap
from _fit_gaussian_mixture_model_ import _fit_gaussian_mixture_model_
from _calculate_spread_per_row_in_dataframe_ import _calculate_spread_per_row_in_dataframe_
from _pick_column_maximizing_spread_ import _pick_column_maximizing_spread_
from _utilities_ import _prepare_spatial_data_for_celltype_assignment_, _save_celltype_assignment_results_to_csv_, _bin_column_in_dataframe_
from _pca_ import _perform_pca_, _plot_pca_
from _generate_anndata_ import _generate_anndata_
import scanpy as sc
import hdf5plugin

if __name__ == '__main__':
    path_to_anndata_files='../data'
    # Prepare spatial data for running the algorithm:
    anndata_allConditions_df, anndata_allConditions_metadata=_prepare_spatial_data_for_celltype_assignment_(path_to_anndata_files) 
    anndata_allConditions_df.to_csv('anndata_allConditions_proteinExpressions.csv')
    anndata_allConditions_metadata.to_csv('anndata_allConditions_metadata.csv')
    adata_df=anndata_allConditions_df.copy()
    # Read in HPA data:
    expr_per_cell_type_max_norm=pd.read_csv('cell_type_nTPM_max_norm.csv').rename(columns={'Unnamed: 0': "Celltypes"}).set_index('Celltypes')
    # Plot heatmap for CD4 and CD8 channels:
    plot_heatmap(expr_per_cell_type_max_norm, ['CD4','CD8A'], outfile='heatmap.jpg', title='Gene-wise max-normalized nTPM per cell type')
    # Enter list of column names as gene symbols. The HPA dataset will not be able to map other formats, including CD codes which most of the columns are in:
    # list_of_channels_present_across_every_sample=['ADAM10', 'ITGAL', 'ITGAX', 'CD14', 'CD163', 'CD24', 'IL2RA', 'ITGB1',
    # 'CD3D', 'CD36', 'CD38', 'CD4', 'CD40', 'CD44', 'PTPRC', 'ICAM1', 'CD55',
    # 'NCAM1', 'CD6', 'SELP', 'CD63', 'CD68', 'CD69', 'CD8A', 'CD9', 'FAS',
    # 'KRT14', 'HLA-A', 'HLA-DRA', 'NOTCH1', 'NOTCH3', 'VIM']
    list_of_channels_present_across_every_sample_=list(anndata_allConditions_df.columns)
    channels_present_in_hpa_data=list(expr_per_cell_type_max_norm)
    unidentified_protein_channels=list(set(list_of_channels_present_across_every_sample_).difference(set(channels_present_in_hpa_data)))
    unidentified_protein_lookup_channels=['CD3', 'CD8', 'HLA-ABC', 'HLA-DR']
    hpa_protein_lookup_channels=[['CD3D', 'CD3E', 'CD3G'][0], ['CD8A', 'CD8B', 'CD8B2'][0], ['HLA-A', 'HLA-B', 'HLA-C'][0], ['HLA-DRA', 'HLA-DRB1', 'HLA-DRB5'][0]]
    unidentified_protein_lookup=dict(zip(unidentified_protein_lookup_channels, hpa_protein_lookup_channels))
    hpa_protein_lookup=dict(zip(hpa_protein_lookup_channels, unidentified_protein_lookup_channels))
    anndata_allConditions_df=anndata_allConditions_df.rename(columns=unidentified_protein_lookup)
    list_of_channels_present_across_every_sample=[unidentified_protein_lookup.get(item,item) for item in list_of_channels_present_across_every_sample_]
    # Plot heatmap for all channels present across every sample:
    plot_heatmap(expr_per_cell_type_max_norm, list_of_channels_present_across_every_sample, outfile='heatmap.jpg', title='Gene-wise max-normalized nTPM per cell type')
    expr_per_cell_type_max_norm=expr_per_cell_type_max_norm[list_of_channels_present_across_every_sample]
    # =================== Celltype assignment algorithm: ===================
    # Steps pending after fitting GMM model, and finding list of good split genes:
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
    clustered_cells_df, clustered_cells_metadata=_save_celltype_assignment_results_to_csv_(clustering_results, anndata_allConditions_metadata, adata_df)
    # Save celltype assigned cells as csv:
    clustered_cells_df.to_csv('clustered_cells_df.csv')
    clustered_cells_metadata.to_csv('clustered_cells_metadata.csv')
    # Bin column in dataframe:
    bin_dict={
        1: '$20-29$',
        2: '$30-39$',
        3: '$40-49$',
        4: '$50-59$',
        5: '$60-69$',
        6: '$70-79$',
        7: '$80-89$'
    }
    clustered_cells_metadata=_bin_column_in_dataframe_(clustered_cells_metadata, bin_dict, 'age', 'ages_binned_column')
    # Save as anndata object:
    adata=_generate_anndata_(clustered_cells_df, clustered_cells_metadata, label='celltype', obs=['celltype', 'condition', 'sex', 'age', 'ages_binned_column'])
    # Save as h5ad file:
    filename='celltype_assigned_anndata.h5ad'
    # Perform PCA using scanpy:
    sc.tl.pca(adata)
    sc.pl.pca(adata, color="celltype")
    # Perform UMAP using scanpy:
    sc.tl.umap(adata)
    sc.pl.umap(adata, color="celltype")
    # Write adata file:
    adata.write_h5ad(
        filename,
        compression=hdf5plugin.FILTERS["zstd"],
        compression_opts=hdf5plugin.Zstd(clevel=5).filter_options
    )
    # Perform PCA using sklearn.decomposition.PCA:
    pca, pca_labeled=_perform_pca_(clustered_cells_df, clustered_cells_metadata, label='celltype')
    pca.to_csv('pca.csv')
    pca_labeled.to_csv('pca_labeled.csv')
    _plot_pca_(pca_labeled)
    
