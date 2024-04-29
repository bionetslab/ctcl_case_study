from sklearn.decomposition import PCA
import pandas as pd
import matplotlib.pyplot as plt

def _plot_pca_(pca_labeled):
    import seaborn as sns
    sns.set_theme(style="whitegrid", palette="muted")
    
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
    ax = sns.scatterplot(data=pca_labeled, x=_x_, y=_y_, hue="celltype", s=2, ax=subplot_axis_id)
    ax.set(xlabel=xlabel, ylabel=ylabel, title='PCA')
    fig.tight_layout()
    # plt.savefig(f'PCA_{xlabel}vs{ylabel}_(hue-celltype).pdf', format='pdf', bbox_inches='tight')
    # plt.show()
    # ---
    # -----------------------
    # ---
    subplot_axis_id=axes["B"]
    # ---
    _x_=0
    _y_=1
    xlabel='PC_'+str(_x_+1)
    ylabel='PC_'+str(_y_+1)
    subplot_axis_id=axes["B"]
    ax = sns.scatterplot(data=pca_labeled, x=_x_, y=_y_, hue="condition", s=2, ax=subplot_axis_id)
    ax.set(xlabel=xlabel, ylabel=ylabel, title='PCA')
    fig.tight_layout()
    plt.savefig('PCA_celltypes_and_conditions.pdf', format='pdf', bbox_inches='tight')
    plt.show()

def _perform_pca_(df, metadata, label=None):
    if label:
        x=df.drop(columns=[label])
    else:
        x=df
    pca_ = PCA()
    x_pca=pca_.fit_transform(x)
    pca=pd.DataFrame(x_pca)
    pca_labeled=pca.copy()
    pca_labeled['x']=list(metadata['x'])
    pca_labeled['y']=list(metadata['y'])
    pca_labeled['sample_id']=list(metadata['sample_id'])
    pca_labeled['celltype']=list(metadata['celltype'])
    pca_labeled['condition']=list(metadata['condition'])
    pca_labeled['age']=list(metadata['age'])
    pca_labeled['sex']=list(metadata['sex'])
    return pca, pca_labeled








