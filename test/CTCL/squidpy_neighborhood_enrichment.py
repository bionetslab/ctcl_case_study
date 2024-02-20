import sys
import os
import time
import scanpy as sc
import scanpy as sc
import squidpy as sq
import numpy as np
import itertools as itt


if __name__ == '__main__':
    cnt=-1
    for filename in os.listdir('data'):
        filename = os.fsdecode(filename)
        if filename.endswith('.h5ad'):
            cnt+=1
            one_or_zero_flag=-1
            adata = sc.read_h5ad('data/'+filename)
            # ---
            # patient_ids.append(i)
            # rec_lab=np.unique(pickle_[i].obsm[dependent_variable_name])[0]
            # recurrence_labels.append(rec_lab)
            # ---
            celltypes=list(sorted(set(adata.obs['celltype'])))
            cluster_pairs_list=list(itt.combinations(celltypes, 2))
            # ---
            
            # ---------------------*************************---------------------
            print(i)
            # squidpy.gr.spatial_neighbors(adata)
            squidpy.gr.spatial_neighbors(adata, coord_type='generic', delaunay=True)
            # squidpy.gr.spatial_neighbors(adata, coord_type='generic', delaunay=False, n_neighs=6)
            nhood_enrichment=squidpy.gr.nhood_enrichment(adata, cluster_key='celltype', copy=True, backend="multiprocessing", show_progress_bar=False)
            # squidpy.gr.nhood_enrichment(adata, cluster_key='celltype', backend="multiprocessing", show_progress_bar=False)
            # sq.pl.nhood_enrichment(adata, cluster_key="celltype")
            nhood_enrichment_zscore=nhood_enrichment[0]
            # nhood_enrichment_zscores.append(nhood_enrichment_zscore)
            # ---
            nhood_enrichment_count=nhood_enrichment[1]
            upper_zscore_matrix = np.triu(np.array(nhood_enrichment_zscore), 1)
            zscore_list=[]
            for j in upper_zscore_matrix:
                for k in j:
                    if k!=0 and k!=None:
                        zscore_list.append(k)
            upper_count_matrix = np.triu(np.array(nhood_enrichment_count), 1)
            count_list=[]
            for j in upper_count_matrix:
                for k in j:
                    if k!=0 and k!=None:
                        count_list.append(k)
            dict_nhoodEnrichment_zscore=dict(zip(cluster_pairs_list, zscore_list))
            dict_nhoodEnrichment_count=dict(zip(cluster_pairs_list, count_list))
            # ---
            Dict_NhoodEnrichment_Zscore.append(dict_nhoodEnrichment_zscore)
            Dict_NhoodEnrichment_Count.append(dict_nhoodEnrichment_count)
            # ---
            
            
            
            
            
            
            
            
            
            
            
            
            if count>0:
                break
            adata = sc.read_h5ad('data/'+filename)
            
            sq.gr.spatial_neighbors(adata, coord_type='generic', delaunay=True)
            nhood_enrichment=sq.gr.nhood_enrichment(adata, cluster_key='celltype', copy=True, backend="multiprocessing", show_progress_bar=False)
            nhood_enrichment_zscore=nhood_enrichment[0]
            nhood_enrichment_count=nhood_enrichment[1]
            
            
            
            
            # start = time.time()
            # all_scores(adata, cluster_key='celltype', radii=[1, 2, 3, 4, 5], num_cell_types=11)
            # end = time.time()
            # adata.uns['SHouT_execution_time'] = end - start
            # adata.write_h5ad(os.path.join('results', filename))