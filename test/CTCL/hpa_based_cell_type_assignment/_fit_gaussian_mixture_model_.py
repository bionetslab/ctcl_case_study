import numpy as np
from sklearn.mixture import GaussianMixture
import statistics
import pandas as pd

def _fit_gaussian_mixture_model_(anndata_allConditions_dict):
    anndata_allConditions_gmm=[] # List to store GMM models.
    anndata_allConditions_gmm_upper_threshold=[] # List of lower threshold for higher 'mode' of bimodal distribution, one per channel.
    anndata_allConditions_gmm_lower_threshold=[] # List of upper threshold for lower 'mode' of bimodal distribution, one per channel.
    good_split_genes=[] # Columns where upper_threshold of lower distribution is lesser than lower threshold of higher distribution.
    bad_split_genes=[] # Columns where above condition is not satisfied.
    
    for column in anndata_allConditions_dict:
        X_list=list(anndata_allConditions_dict[column]) # Column value used to fit GMM model.
        X=np.array(X_list).reshape(-1,1) # Series to fit GMM model on.
        gm = GaussianMixture(n_components=2, random_state=0).fit(X)
        anndata_allConditions_gmm.append(gm)
        labels = gm.predict(X)
        X_df=pd.DataFrame() # Create DataFrame of series with corresponding label (i.e., whether data point belongs to lower or higher distribution).
        X_df['X']=X_list
        X_df['labels']=labels
        X_0=list(X_df[X_df['labels']==0]['X']) # Data points in series belonging to lower distribution.
        X_1=list(X_df[X_df['labels']==1]['X']) # Data points in series belonging to higher distribution.
        
        means_=gm.means_
        stddevs_0= statistics.pstdev(X_0)
        stddevs_1= statistics.pstdev(X_1)
        stddevs_=[stddevs_0, stddevs_1]
        
        upper_threshold=means_[1]-1.96*stddevs_[1]
        anndata_allConditions_gmm_upper_threshold.append(upper_threshold)
        lower_threshold=means_[0]+1.96*stddevs_[0]
        anndata_allConditions_gmm_lower_threshold.append(lower_threshold)
        
        if upper_threshold>lower_threshold:
            good_split_genes.append(column)
        else:
            bad_split_genes.append(column)
    
    geneName_upperThreshold_dict=dict(zip(list(anndata_allConditions_dict.keys()), anndata_allConditions_gmm_upper_threshold))
    
    return good_split_genes, geneName_upperThreshold_dict
