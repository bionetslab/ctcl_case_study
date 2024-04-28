import seaborn as sns
import scanpy as sc

def _plot_pca_(pca_labeled, hue, axes, legend_loc=None, title_prefix=None, plot_title=None, title_loc=None, xlabel=None, palette=None, alpha=None, s=0.25):
    sns.set_theme(style="whitegrid", palette="muted")
    _x_=0
    _y_=1
    xlabel='PC '+str(_x_+1)
    ylabel='PC '+str(_y_+1)
    ax = sns.scatterplot(data=pca_labeled, x=str(_x_), y=str(_y_), hue=hue, s=s, ax=axes, palette=palette, alpha=alpha)
    # ax.set(xlabel=xlabel, ylabel=ylabel, title='PCA')
    
    sns.set(font_scale = 1.2)
    sns.set_style("white")
    sns.set(font_scale = 1.2)
    sns.set_style("white")
    if title_prefix:
        if plot_title:
            plot_title=title_prefix + plot_title
        else:
            plot_title=title_prefix
    ax.set_title(plot_title, loc=title_loc, fontsize=30, fontweight='bold') # fontsize=25, pad=20
    ax.set_xlabel(xlabel, fontsize=30) # fontsize=25, labelpad=20
    ax.set_ylabel(ylabel, fontsize=30) # fontsize=25, labelpad=20
    # ax.get_legend().set_visible(False)
    # ax.legend.remove()
    yticks=[]
    for j in ax.get_yticks():
        yticks.append(round(j,1))
    ax.set_yticklabels(yticks, size = 25) # size = 20
    ax.tick_params(axis='x', which='major', labelsize=25)
    ax.grid(False)
    sns.move_legend(ax, loc='upper right', title=None, fontsize=25, markerscale=4)
    
def _scanpy_pca_(pca_labeled, hue, axes, legend_loc=None, title_prefix=None, plot_title=None, title_loc=None, xlabel=None, palette=None, alpha=None, s=5):
    sns.set_theme(style="whitegrid", palette="muted")
    _x_=0
    _y_=1
    xlabel='PC '+str(_x_+1)
    ylabel='PC '+str(_y_+1)
    ax = sc.pl.pca(pca_labeled, color=hue, show=False, size=s, ax=axes, palette=palette, alpha=alpha, title='')
    sns.set(font_scale = 1.2)
    sns.set_style("white")
    if title_prefix:
        if plot_title:
            plot_title=title_prefix + plot_title
        else:
            plot_title=title_prefix
    ax.set_title(plot_title, loc=title_loc, fontsize=30, fontweight='bold') # fontsize=25, pad=20
    ax.set_xlabel(xlabel, fontsize=30) # fontsize=25, labelpad=20
    ax.set_ylabel(ylabel, fontsize=30) # fontsize=25, labelpad=20
    yticks=[]
    for j in ax.get_yticks():
        yticks.append(round(j,1))
    ax.set_yticklabels(yticks, size = 25) # size = 20
    ax.tick_params(axis='x', which='major', labelsize=25)
    ax.grid(False)
    sns.move_legend(ax, loc='upper right', title=None, fontsize=25, markerscale=4)

def _scanpy_umap_(pca_labeled, hue, axes, legend_loc=None, title_prefix=None, plot_title=None, title_loc=None, xlabel=None, palette=None, alpha=None, s=5):
    sns.set_theme(style="whitegrid", palette="muted")
    _x_=0
    _y_=1
    xlabel='UMAP '+str(_x_+1)
    ylabel='UMAP '+str(_y_+1)
    ax = sc.pl.umap(pca_labeled, color=hue, show=False, size=s, ax=axes, palette=palette, alpha=alpha, title='')
    sns.set(font_scale = 1.2)
    sns.set_style("white")
    if title_prefix:
        if plot_title:
            plot_title=title_prefix + plot_title
        else:
            plot_title=title_prefix
    ax.set_title(plot_title, loc=title_loc, fontsize=30, fontweight='bold') # fontsize=25, pad=20
    ax.set_xlabel(xlabel, fontsize=30) # fontsize=25, labelpad=20
    ax.set_ylabel(ylabel, fontsize=30) # fontsize=25, labelpad=20
    yticks=[]
    for j in ax.get_yticks():
        yticks.append(round(j,1))
    ax.set_yticklabels(yticks, size = 25) # size = 20
    ax.tick_params(axis='x', which='major', labelsize=25)
    ax.grid(False)
    sns.move_legend(ax, loc='upper right', title=None, fontsize=25, markerscale=4)
    
    
    
    
    
    