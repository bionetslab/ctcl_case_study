import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import statistics
from sklearn.decomposition import PCA
sns.set_theme(style="whitegrid", palette="muted")

if __name__ == '__main__':
    anndata_allConditions_df=pd.read_csv('../data/anndata_allConditions_proteinExpressions.csv')
    anndata_allConditions_metadata=pd.read_csv('../hpa_based_cell_type_assignment/HPA_based_celltype_assignment_results.csv')
    pca = PCA()
    Xt=pca.fit_transform(anndata_allConditions_df)
    pca_df=pd.DataFrame(Xt)
    pca_df_labeled=pca_df.copy()
    pca_df_labeled['celltype']=list(anndata_allConditions_metadata['cell_type'])
    pca_df_labeled['Labels']=list(anndata_allConditions_metadata['condition'])
    
    layout = [
        ["A"],
        ["B"]
    ]
    fig, axes = plt.subplot_mosaic(layout, layout='constrained', figsize=(20,10))
    
    subplot_axis_id=axes["A"]
    _x_=0
    _y_=1
    xlabel='PC_'+str(_x_+1)
    ylabel='PC_'+str(_y_+1)
    ax = sns.scatterplot(data=pca_df_labeled, x=_x_, y=_y_, hue="celltype", s=2, ax=subplot_axis_id)
    ax.set(xlabel=xlabel, ylabel=ylabel, title='PCA')
    fig.tight_layout()
    subplot_axis_id=axes["B"]
    _x_=0
    _y_=1
    xlabel='PC_'+str(_x_+1)
    ylabel='PC_'+str(_y_+1)
    subplot_axis_id=axes["B"]
    ax = sns.scatterplot(data=pca_df_labeled, x=_x_, y=_y_, hue="Labels", s=2, ax=subplot_axis_id)
    ax.set(xlabel=xlabel, ylabel=ylabel, title='PCA')
    fig.tight_layout()
    plt.savefig(f'PCA_celltypes_and_conditions.pdf', format='pdf', bbox_inches='tight')
    plt.show()
    pca_Labels_dict_perFeature={}
    for j in list(pca_df.columns):
        pca_Labels_dict={}
        _labels_=list(np.unique(pca_df_labeled.Labels))
        pca_df_Labels_stdDev=statistics.pstdev(list(pca_df_labeled[j]))
        for i in _labels_:
            pca_Labels_dict[i]=pca_df_Labels_stdDev / statistics.pstdev(list(pca_df_labeled[pca_df_labeled['Labels']==i][0]))
        pca_Labels_dict_perFeature[j]=pca_Labels_dict
    pca_celltype_dict_perFeature={}
    for j in list(pca_df.columns):
        pca_celltype_dict={}
        _labels_=list(np.unique(pca_df_labeled.celltype))
        pca_df_celltype_stdDev=statistics.pstdev(list(pca_df_labeled[j]))
        for i in _labels_:
            pca_celltype_dict[i]=pca_df_celltype_stdDev / statistics.pstdev(list(pca_df_labeled[pca_df_labeled['celltype']==i][0]))
        pca_celltype_dict_perFeature[j]=pca_celltype_dict
    DF_pca_celltype_dict_perFeature=pd.DataFrame(pca_celltype_dict_perFeature).T
    DF_pca_Labels_dict_perFeature=pd.DataFrame(pca_Labels_dict_perFeature).T
    DF_PCA=pd.concat([DF_pca_celltype_dict_perFeature, DF_pca_Labels_dict_perFeature], axis=1)
    DF_PCA.index.names = ['PC']
    sns.set_theme()
    sns.set_style("whitegrid")
    plt.figure(figsize=(10,6), tight_layout=True)
    for i in list(DF_PCA.columns):
        print(i)
        if i in list(DF_pca_Labels_dict_perFeature.columns):
            plt.plot(DF_PCA[[i]], '--', linewidth=2)
        else:
            plt.plot(DF_PCA[[i]], 'o-', linewidth=2)
    plt.xlabel('Principal component #')
    plt.ylabel('Relative homogeneity score\n(Cluster quality score)')
    plt.title('Cluster quality across principal components')
    plt.legend(title='Celltypes/ conditions', title_fontsize = 13, labels=list(DF_PCA.columns))
    plt.savefig('lineplot_ClusterQuality_CelltypesAndConditions(Labels).pdf', format='pdf', bbox_inches='tight')
    plt.show()
    # -----------