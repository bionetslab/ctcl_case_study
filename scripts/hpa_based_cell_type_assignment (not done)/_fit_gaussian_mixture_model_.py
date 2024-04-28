import numpy as np
from sklearn.mixture import GaussianMixture
import statistics
import pandas as pd

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
