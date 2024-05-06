import pandas as pd
import os
import seaborn as sns
sns.set_theme(style="whitegrid")
import scanpy as sc
import time
import sys
sys.path.insert(0, '../../../')
from SHouT.shout import all_scores
sns.set_theme(style="whitegrid")
radii=[1, 5, 10, 20, 40, 80, 100]

times=[]
for filename in os.listdir('../../../data'):
    filename = os.fsdecode(filename)
    if filename.endswith('.h5ad'):
        print(filename)
        adata=sc.read_h5ad('../../../data/'+filename)
        times_=pd.DataFrame()
        for radius in radii:
            start = time.time()
            all_scores(adata, cluster_key='celltype', radii=[radius])
            end = time.time()
            times_[radius]=[end-start]
        times.append(times_)
times=pd.concat(times, axis=0)
data_all_cells=pd.melt(times)
data_all_cells.to_csv('../../../results/data_all_cells.csv')