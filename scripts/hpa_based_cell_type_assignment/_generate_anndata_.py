import scanpy as sc
import numpy as np

def _generate_anndata_(df, metadata, label=None, obs=None):
    if label:
        x=df.drop(columns=[label])
    else:
        x=df
    adata=sc.AnnData(x)
    adata_obs=metadata[obs]
    for i in obs:
        adata.obs[i]=np.array(adata_obs[i])
    return adata