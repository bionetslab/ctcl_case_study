import os
import scanpy as sc
import squidpy as sq
import numpy as np
import itertools as itt
import pandas as pd

if __name__ == '__main__':
    cnt=-1
    sq_centrality_scores_results=[]
    for filename in os.listdir('data'):
        filename = os.fsdecode(filename)
        if filename.endswith('.h5ad'):
            cnt+=1
            one_or_zero_flag=-1
            adata = sc.read_h5ad('data/'+filename)
            celltypes=list(sorted(set(adata.obs['celltype'])))
            nshape=len(celltypes)
            sq.gr.spatial_neighbors(adata, coord_type='generic', delaunay=True)
            centrality_scores=sq.gr.centrality_scores(adata, cluster_key='celltype', copy=True, backend="multiprocessing", show_progress_bar=False)
            centrality_scores = centrality_scores.reset_index().rename(columns={"index":"celltypes"})
            centrality_scores['patient_id']=[adata.uns['patient_id'][0]]*nshape
            centrality_scores['condition']=[adata.uns['Group'][0]]*nshape
            sq_centrality_scores_results.append(centrality_scores)
    centrality_scores=pd.concat(sq_centrality_scores_results, axis=0)
    centrality_scores.to_csv(os.path.join('results', 'squidpy_centralityScores_results.csv'))
    
    
    
    