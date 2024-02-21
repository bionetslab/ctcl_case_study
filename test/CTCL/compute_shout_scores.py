import sys
sys.path.insert(0, '../../')
from shout import all_scores
import os
import time
import scanpy as sc



if __name__ == '__main__':
    count=-1
    for filename in os.listdir('data'):
        filename = os.fsdecode(filename)
        if filename.endswith('.h5ad'):
            count+=1
            if count>0:
                break
            adata = sc.read_h5ad('data/'+filename)
            start = time.time()
            all_scores=all_scores(adata, cluster_key='celltype', radii=[1, 2, 3, 4, 5], num_cell_types=11, copy=True)
            end = time.time()
            adata.uns['SHouT_execution_time'] = end - start
            adata.write_h5ad(os.path.join('results', filename))
