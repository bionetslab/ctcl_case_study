import os
import scanpy as sc
import squidpy as sq
import numpy as np
import itertools as itt
import pandas as pd

if __name__ == '__main__':
    cnt=-1
    sq_ne_results=[]
    for filename in os.listdir('../../data'):
        filename = os.fsdecode(filename)
        if filename.endswith('.h5ad'):
            cnt+=1
            one_or_zero_flag=-1
            adata = sc.read_h5ad('../../data/'+filename)
            celltypes=list(sorted(set(adata.obs['celltype'])))
            cluster_pairs_list=list(itt.combinations(celltypes, 2))
            nshape_interactions=len(cluster_pairs_list)
            sq.gr.spatial_neighbors(adata, coord_type='generic', delaunay=True)
            nhood_enrichment=sq.gr.nhood_enrichment(adata, cluster_key='celltype', copy=True, backend="multiprocessing", show_progress_bar=False)
            nhood_enrichment_zscore=nhood_enrichment[0]
            nhood_enrichment_count=nhood_enrichment[1]
            n_shape=np.shape(nhood_enrichment_zscore)[0]
            upper_zscore_matrix = np.triu(np.array(nhood_enrichment_zscore), 1)
            zscore_list=[]
            for j in range(1, n_shape):
                for k in range(j, n_shape):
                    zscore_list.append(upper_zscore_matrix[j,k])
            upper_count_matrix = np.triu(np.array(nhood_enrichment_count), 1)
            count_list=[]
            for j in range(1, n_shape):
                for k in range(j, n_shape):
                    count_list.append(upper_count_matrix[j,k])
            df=pd.DataFrame(cluster_pairs_list, columns =['celltype_1', 'celltype_2'])
            df['sq_ne_zscore']=zscore_list
            df['sq_ne_count']=count_list
            df['patient_id']=[adata.uns['patient_id'][0]]*nshape_interactions
            df['condition']=[adata.uns['Group'][0]]*nshape_interactions
            sq_ne_results.append(df)       
    sq_ne_results=pd.concat(sq_ne_results, axis=0)
    sq_ne_results.to_csv(os.path.join('../../results', 'squidpy_nhoodEnrichment_results.csv'))
    
    
    
    
    